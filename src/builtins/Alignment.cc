#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
#include "computation/computation.H"
#include "dp/2way.H"
#include "imodel/imodel.H"
#include "computation/expression/expression.H"
#include <boost/dynamic_bitset.hpp>
#include "alignment/alignment.H"
#include <boost/optional.hpp>
#include <tuple>

using boost::optional;
using std::tuple;

extern "C" closure builtin_function_pairwise_alignment_probability_from_counts(OperationArgs& Args)
{
    const matrix<int>& counts = Args.evaluate(0).as_<Box<matrix<int>>>();
    const indel::PairHMM& Q = Args.evaluate(1).as_<indel::PairHMM>();

    using namespace A2;

    log_double_t P=1;

    // Account for S-? start probability
    for(int i=0;i<Q.size2();i++)
	if (counts(states::S,i))
	    P *= Q.start(i);

    // Account for the mass of transitions
    for(int i=0;i<3;i++)
	for(int j=0;j<3;j++) {
	    log_double_t Qij = Q(i,j);
	    // FIXME - if we propose really bad indel parameters, we can get log(Q_ij) where Qij == 0
	    if (counts(i,j))
		P *= pow(Qij,counts(i,j));
	}
  
    // Account for ?-E end probability
    if (not counts(states::S,states::E))
	for(int i=0;i<Q.size1();i++)
	    if (counts(i,states::E))
		P *= Q(i,states::E);

    return {P};
}

extern "C" closure builtin_function_pairwise_alignment_length1(OperationArgs& Args)
{
    return {Args.evaluate(0).as_<pairwise_alignment_t>().length1()};
}

extern "C" closure builtin_function_pairwise_alignment_length2(OperationArgs& Args)
{
    return {Args.evaluate(0).as_<pairwise_alignment_t>().length2()};
}

extern "C" closure builtin_function_transition_counts(OperationArgs& Args)
{
    const pairwise_alignment_t& A = Args.evaluate(0).as_<pairwise_alignment_t>();

    Box<matrix<int>> counts(5,5,0);

    using namespace A2;

    int prev = states::S;
    for(int column=0;column<A.size();column++)
    {
	counts(prev,A.get_state(column))++;
	prev = A.get_state(column);
    }
    counts(prev, states::E)++;

    return counts;
}

using std::vector;

extern "C" closure builtin_function_rs07_lengthp(OperationArgs& Args)
{
    double e = Args.evaluate(0).as_double();
    if (e < 0.0)
	throw myexception()<<"Error: mean indel length cannot be < 1, but was set to "<<1.0/(1.0-e)<<"!";

    int l = Args.evaluate(1).as_int();

    if (l < 0)
	return {0.0};
    else if (l==0)
	return {1.0};
    else
	return {1.0-e};
}

extern "C" closure builtin_function_rs07_branch_HMM(OperationArgs& Args)
{
    double e = Args.evaluate(0).as_double();
    if (e < 0.0)
	throw myexception()<<"Error: mean indel length cannot be < 1, but was set to "<<1.0/(1.0-e)<<"!";

    double D = Args.evaluate(1).as_double();
    double heat = Args.evaluate(2).as_double();
    constructor in_training_c = Args.evaluate(3).head().as_<constructor>();
    bool in_training = true;
    if (in_training_c.f_name == "Prelude.False")
	in_training = false;

    using namespace A2::states;

    // Here D = rate * t

    // Return a model with all probabilities zero if e==1.
    // Scaling time by 1/(1.0-e) doesn't work if e==1.
    if (e >= 1)
	return indel::PairHMM();

    // (1-e) * delta / (1-delta) = P(indel)
    // But move the (1-e) into the RATE to make things work
    double mu = D/(1.0-e);
    double P_indel = 1.0 - exp(-mu);
    double A = P_indel;

    if (in_training) A = std::min(A,0.005);

    double delta = A/(1+A);

    // Note: If the branch is disconnected, then t < -0.5
    //  if (t < -0.5) delta = 0.5;

    double f = 0.1; //unaligned fraction
    delta = pow(delta, heat) * pow(f/(1+f),1-heat);
    e = 1.0 - pow(1.0 - e, heat);

    if (1 - 2*delta <0)
	throw myexception()<<"indel model: we need (delta <= 0.5), but delta = "<<delta;

    if (e > 1)
	throw myexception()<<"indel model: we need (epsilon <= 1), but epsilon = "<<e;
    
    assert(delta >= 0 and delta <= 1);
    assert(e >= 0 and e < 1);

    // transition probabilities default to *zero*
    indel::PairHMM Q;

    Q(S ,S ) = 0;
    Q(S ,M ) = 1 - 2*delta;
    Q(S ,G1) = delta;
    Q(S ,G2) = delta;
    Q(S ,E ) = 1 - delta;

    Q(M ,S ) = 1;
    Q(G1,S ) = 1;
    Q(G2,S ) = 1;

    // turn the model into a fragment model
    fragmentize(Q,e);

    remove_one_state(Q,S);

    Q.start_pi(S)  = 0;
    Q.start_pi(M)  = 1;
    Q.start_pi(G1) = 0;
    Q.start_pi(G2) = 0;
    Q.start_pi(E)  = 0;

    return Q;
}

extern "C" closure builtin_function_rs05_branch_HMM(OperationArgs& Args)
{
    double e = Args.evaluate(0).as_double();
    double delta = Args.evaluate(1).as_double();
    double t = Args.evaluate(2).as_double();
    double heat = Args.evaluate(3).as_double();
    constructor in_training_c = Args.evaluate(4).head().as_<constructor>();
    bool in_training = true;
    if (in_training_c.f_name == "Prelude.False")
	in_training = false;

    if (in_training) delta = std::min(delta,0.005);

  // Return a model with all probabilities zero if e==1.
    if (e >= 1)
	return indel::PairHMM();

    double f = 0.1; //unaligned fraction
    delta = pow(delta, heat) * pow(f/(1+f),1-heat);
    e = 1.0 - pow(1.0 - e, heat);

    if (delta > 0.5)
	throw myexception()<<"RS05_branch_HMM: we need (delta <= 0.5), but delta = "<<delta;

    if (e > 1.0)
	throw myexception()<<"RS05_branch_HMM: we need (epsilon <= 1), but epsilon = "<<e;

    assert(delta >= 0.0 and delta <= 1.0);
    assert(e >= 0.0 and e < 1.0);

    indel::PairHMM Q;
    using namespace A2::states;

    Q(S ,S ) = 0;
    Q(S ,M ) = 1 - 2*delta;
    Q(S ,G1) = delta;
    Q(S ,G2) = delta;
    Q(S ,E ) = 0;

    Q(M ,S ) = 1;
    Q(G1,S ) = 1;
    Q(G2,S ) = 1;

    // For the states G1, G2 fragment lengths are Geometric(e)
    fragmentize(Q,e,G1);
    fragmentize(Q,e,G2);

    // For the states M, G1, G2 we might exit with probability t
    exitize(Q,t,M ,E);
    exitize(Q,t,G1,E);
    exitize(Q,t,G2,E);

    // When moving from another state, continue until we are not in S
    remove_one_state(Q,S);

    Q.start_pi(S)  = 0;
    Q.start_pi(M)  = 1;
    Q.start_pi(G1) = 0;
    Q.start_pi(G2) = 0;
    Q.start_pi(E)  = 0;

    return Q;
}

// f_M(s) = [ ME  + s(MGxGE - MExGG) ] / [ 1 - s(GG + MM) + s^2(MMxGG - MGxGM) ]

extern "C" closure builtin_function_rs05_lengthp(OperationArgs& Args)
{
    indel::PairHMM QE = Args.evaluate(0).as_<indel::PairHMM>();
    int l = Args.evaluate(1).as_int();

    using namespace A2::states;

    remove_one_state(QE,G2);

    //--------------- Remove the 'G2' State ----------------------//
    double MM = QE(M,M);
    double MG = QE(M,G1);
    double ME = QE(M,E);

    double GM = QE(G1,M);
    double GG = QE(G1,G1);
    double GE = QE(G1,E);

    //----- Calculate roots of q(s); we assume its quadratic -----//
    double C = 1;
    double B = -(GG + MM);
    double A = MM*GG - MG*GM;
    if (A == 0) return {0.0};
    
    double sqr_det = sqrt(B*B-4.0*A*C);
    double r1 = (-B - sqr_det)/(2*A);
    double r2 = (-B + sqr_det)/(2*A);

    //------------ Calculate the coefficients f_M[l] ------------//
    double P;
    if (l==0)
	P = ME;
    else {
	double P1 = pow(r1,-l-1);
	double P2 = pow(r2,-l-1);

	// Calculate q[l] and q[l-i] (I've proved that all q[i]>=0)
	double q_l   = 1.0/ (A*(r2-r1)) * (P1 - P2);
	double q_lm1 = 1.0/ (A*(r2-r1)) * (P1*r1 - P2*r2);

	// Calculate f_M[l] from the q[i] (*IS* this always positive?)
	P = ME*q_l + (MG*GE - ME*GG)*q_lm1;
    }
    return {P};
}

extern "C" closure builtin_function_bitmask_from_alignment(OperationArgs& Args)
{
    using boost::dynamic_bitset;

    auto arg0 = Args.evaluate(0);
    const auto& A = arg0. as_<alignment>();

    int seq = Args.evaluate(1).as_int();

    int L = A.length();
    
    object_ptr<Box<dynamic_bitset<>>> mask_(new Box<dynamic_bitset<>>(L));
    auto& mask = *mask_;

    for(int i=0;i<L;i++)
	if (A.character(i,seq))
	    mask.flip(i);

    return mask_;
}

extern "C" closure builtin_function_load_alignment(OperationArgs& Args)
{
    std::string filename = Args.evaluate(0).as_<String>();

    object_ptr<alignment> A(new alignment(DNA(),filename));

    return A;
}

#include "alignment/alignment-constraint.H"

extern "C" closure builtin_function_leaf_alignment_constraint(OperationArgs& Args)
{
    // 1. Read the arguments

    // 1a. Matrix of all columns to constrain
    auto M_ = Args.evaluate(0);
    auto& M = M_.as_<Box<matrix<int>>>();

    // 1b: width of constraint
    int delta = Args.evaluate(1).as_int();

    // 1c: which sequence to get the constraint for
    int s = Args.evaluate(2).as_int();

    // 1d: sequence length
    int L = Args.evaluate(3).as_<Vector<int>>().size();

    // 2. Construct the constraint object to return
    object_ptr<alignment_constraints> leaf_con( new alignment_constraints );

    for(int i=0;i<M.size1();i++)
    {
	int j = M(i,s);
	if (j < 0) continue;

	optional<int> j_minus_delta = j - delta;
	optional<int> j_plus_delta = j + delta;

	if (*j_minus_delta < 0) j_minus_delta = boost::none;
	if (*j_plus_delta >= L) j_plus_delta  = boost::none;

	leaf_con->push_back(alignment_constraint{i, j_minus_delta, j_plus_delta, 1});
    }

    return leaf_con;
}

template<typename T>
boost::optional<T> max(const boost::optional<T>& x, const boost::optional<T>& y)
{
    if (not x) return y;
    if (not y) return x;
    return std::max(*x,*y);
}

template<typename T>
boost::optional<T> min(const boost::optional<T>& x, const boost::optional<T>& y)
{
    if (not x) return y;
    if (not y) return x;
    return std::min(*x,*y);
}

extern "C" closure builtin_function_merge_alignment_constraints(OperationArgs& Args)
{
    using boost::get;

    // 1. Read the arguments
    auto con_x_ = Args.evaluate(0);
    auto& con_x = con_x_.as_<alignment_constraints>();

    auto a_xz_ = Args.evaluate(1);
    auto& a_xz = a_xz_.as_<pairwise_alignment_t>();

    auto con_y_ = Args.evaluate(2);
    auto& con_y = con_y_.as_<alignment_constraints>();

    auto a_yz_ = Args.evaluate(3);
    auto& a_yz = a_yz_.as_<pairwise_alignment_t>();

    // 2. Construct the object to return
    object_ptr<alignment_constraints> con_z( new alignment_constraints );

    auto max_z_le_x = get_max_y_le_x(a_xz);
    auto max_z_le_y = get_max_y_le_x(a_yz);

    auto min_z_ge_x = get_min_y_ge_x(a_xz);
    auto min_z_ge_y = get_min_y_ge_x(a_yz);

    con_z->reserve(con_x.size()+con_y.size());

    for(int i=0,j=0;i<con_x.size() or j<con_y.size();)
    {
	int id = std::min(get<0>(con_x[i]), get<0>(con_y[j]) );

	optional<int> zmax_x;
	optional<int> zmax_y;

	optional<int> zmin_x;
	optional<int> zmin_y;

	bool have_x_con = i < con_x.size() and get<0>(con_x[i]) == id;
	bool have_y_con = j < con_y.size() and get<0>(con_y[j]) == id;

	int xnum = 0;
	int ynum = 0;

        // 3a. Get the zmax(X-Delta) and zmax(X+Delta)
	if (have_x_con)
	{
	    auto zmax_x = lookup(max_z_le_x, get<1>(con_x[i]));
	    auto zmin_x = lookup(min_z_ge_x, get<2>(con_x[i]));
	    xnum = get<3>(con_x[i]);
	    i++;
	}

	// 3b. Get the zmax(X-Delta) and zmax(X+Delta)
	if (have_y_con)
	{
	    auto zmax_y = lookup(max_z_le_y, get<1>(con_y[j]));
	    auto zmin_y = lookup(min_z_ge_y, get<2>(con_y[j]));
	    ynum = get<3>(con_y[i]);
	    j++;
	}

	// 3c. Check that the constraint is met on the X and Y sequences.
	if (zmax_x and zmin_y and not (*zmax_x < *zmin_y) ) throw myexception()<<"X-D !>= Y+D constraints failed!";
	if (zmax_y and zmin_x and not (*zmax_y < *zmin_x) ) throw myexception()<<"Y-D !>= X+D constraints failed!";

	// max_{xy in (X U Y)-Delta} max {z <= xy}
	auto zmax = std::max(zmax_x, zmax_y);
	// min_{xy in (X U Y)+Delta} min {z >= xy}
	auto zmin = std::min(zmin_x, zmin_y);
	int znum = xnum + ynum;

	con_z->push_back(alignment_constraint{id, zmax, zmin, znum});
    }

    return con_z;
}
