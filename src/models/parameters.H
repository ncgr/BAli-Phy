/*
  Copyright (C) 2004-2012 Benjamin Redelings

  This file is part of BAli-Phy.

  BAli-Phy is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free
  Software Foundation; either version 2, or (at your option) any later
  version.

  BAli-Phy is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with BAli-Phy; see the file COPYING.  If not see
  <http://www.gnu.org/licenses/>.  */

///
/// \file parameters.H
///
/// \brief This file defines the Parameters class which holds the model and state
///        for the MCMC.
///

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <valarray>
#include <boost/optional.hpp>

#include "sequence/alphabet.H"
#include "tree/sequencetree.H"
#include "models/model.H"
#include "imodel/imodel.H"
#include "alignment/alignment.H"

#include "tools/partition.H"
#include "cow-ptr.H"
#include "cached_value.H"
#include "dp/2way.H"
#include "TreeInterface.H"

//------------------------------- parameter-containing class --------------------------//

class Parameters;
class Likelihood_Cache_Branch;
class model_t;

struct data_partition_constants
{
    /// The pairwise alignment for each branch
    std::vector<int> pairwise_alignment_for_branch;

    int likelihood_index;

    int cl_index;

    std::vector<int> conditional_likelihoods_for_branch;

    std::vector<int> leaf_sequence_indices;

    std::vector<int> sequence_length_indices;

    std::vector<int> sequence_length_pr_indices;

    std::vector<int> transition_p_method_indices;

    std::vector<int> branch_HMM_indices;

    int alignment_prior_index = -1;

    /// The names and comments for the sequences
    std::vector< sequence > seqs;

    /// The leaf sequence data of this partition
    std::vector< std::vector<int> > sequences;

    boost::shared_ptr<alphabet> a;

    std::vector<int> branch_HMM_type;

    /// Alignment constraint
    matrix<int> alignment_constraint;

    int likelihood_calculator;

    data_partition_constants(Parameters* P, int i, const alignment&, const std::vector<int>& counts, int like_calc);
};

/// Each data_partition is a model with one parameter: mu (the branch mean)
class data_partition
{
    friend class Parameters;
  
public:
    const Parameters* P;

protected:
    int partition_index;

    const data_partition_constants& DPC() const;

public:

    /// Is the alignment allowed to vary?
    bool variable_alignment() const;
  
    int likelihood_calculator() const;

    alignment A() const;

    TreeInterface t() const;

    /// Cached Transition Probabilities
    const EVector& transition_P(int b) const;
  
    const std::vector<int>& get_sequence(int i) const;

    int n_base_models() const;

    // Perhaps this should really be per-mixture-component.
    int n_states() const;

    // Perhaps this should also be per-mixture-component.
    std::vector<unsigned> state_letters() const;
  
    std::vector<double> distribution() const;

    Matrix WeightedFrequencyMatrix() const;

    Matrix FrequencyMatrix() const;

    const Likelihood_Cache_Branch& cache(int b) const;

    int subst_root() const;

    /// cached branch HMMs
    const indel::PairHMM& get_branch_HMM(int b) const;
    std::vector<indel::PairHMM> get_branch_HMMs(const std::vector<int>&) const;

    double sequence_length_pr(int i) const;

    double get_beta() const;

    // get the alphabet for partition i
    const alphabet& get_alphabet() const;

    bool has_IModel() const;

    int seqlength(int n) const;

    expression_ref get_pairwise_alignment_(int b) const;

    const pairwise_alignment_t& get_pairwise_alignment(int b) const;

    bool pairwise_alignment_is_unset(int b) const;

    const matrix<int>& alignment_constraint() const;

    // The prior, likelihood, and probability
    log_double_t prior_alignment() const;
    log_double_t likelihood() const;
    log_double_t heated_likelihood() const;

    data_partition(const Parameters* P, int i);
};

struct mutable_data_partition: public data_partition
{
    void set_alignment(const alignment&);

    void unalign_sequences();

    void set_pairwise_alignment(int b, const pairwise_alignment_t&);

    void unset_pairwise_alignment(int b);

    mutable_data_partition(const Parameters* P, int i);
};

struct smodel_methods
{
    int main;
    int n_base_models;
    int n_states;
    int state_letters;
    int distribution;
    int weighted_frequency_matrix;
    int frequency_matrix;
    int get_alphabet;
    int rate;

    // but these ones take an argument!
    int transition_p; // int branch
    int frequencies; // int 
  
    // also:
    // Q1: can we use evaluate_expression( ) to pass in arguments?  what are the pitfalls?
    // Q2: how could we control the peeling from haskell?
    // 
    smodel_methods(const expression_ref&, context&);
};

class maybe_parameter
{
    int param = -1;
    int value = -1;
public:
    int get_value(const context* C) const
	{
	    if (is_const())
		return value;
	    else
		return C->get_parameter_value(param).as_int();
	}
    void set_value(context* C, int v) { assert(not is_const()); C->set_parameter_value(param,v);}
    bool is_const() const {return param == -1;}
    maybe_parameter(int p):param(p) {}
    maybe_parameter(int,int v):value(v) {}
};

struct tree_constants: public Object
{
    // The tree
    int tree_head = -1;
    std::vector<std::vector<maybe_parameter>> parameters_for_tree_node;
    std::vector<std::tuple<maybe_parameter,maybe_parameter,maybe_parameter>> parameters_for_tree_branch;

    std::vector<int> branch_duration_index;
    std::vector<boost::optional<int>> branch_duration_regs;

    int n_leaves = -1;

    std::vector<std::string> node_labels;
  
    tree_constants* clone() const {return new tree_constants(*this);}

    tree_constants(Parameters*, const SequenceTree&, const model_t&);
};

struct parameters_constants: public Object
{
    /// The substitution::Model
    std::vector<smodel_methods> SModels;
    std::vector<boost::optional<int>> smodel_for_partition;

    /// The IndelModel
    std::vector<boost::optional<int>> imodel_for_partition;
    int n_imodels;
    int imodels_index = -1;

    // Scale for partition
    std::vector<boost::optional<int>> scale_for_partition;
    int n_scales;
    std::vector<int> scale_parameter_indices;

    std::vector< std::vector<int> > branch_length_indices;

    int substitution_branch_lengths_index = -1;

    matrix<int> branch_transition_p_indices;

    /// Do we have an Exponential (0) or Gamma-0.5 (1) prior on branch lengths?
    int branch_prior_type = 0;

    /// Topology and Branch-Length Constraint
    SequenceTree TC;

    /// Alignment-constrained branches
    std::vector<int> AC;

    std::vector<int> branch_HMM_type;

    /// Tree partitions to weight
    std::vector<Partition> partitions;
    std::vector<log_double_t> partition_weights;

    /// Is pinning enabled or not?
    std::vector<double> constants;
    /// Is pinning enabled or not?
    int features = 0;

    std::vector<double> beta_series;

    std::vector<double> all_betas;

    std::vector<data_partition_constants> DPC;

    parameters_constants* clone() const {return new parameters_constants(*this);}

    parameters_constants(const std::vector<alignment>& A, const SequenceTree&, 
			 const std::vector<model_t>&,
			 const std::vector<boost::optional<int>>&,
			 const std::vector<model_t>&,
			 const std::vector<boost::optional<int>>&,
			 const std::vector<boost::optional<int>>&);
};

/// A class to contain all the MCMC state except the alignment
class Parameters: public Model
{
    friend class data_partition;
    friend class TreeInterface;
public:
    object_ptr<parameters_constants> PC;
    object_ptr<tree_constants> TC;

    int beta_index;

    int subst_root_index;

    bool variable_alignment_;
    int variable_alignment_param;

private:

    std::vector<int> affected_nodes;
    std::vector<std::vector<int>*> branches_from_affected_node;

    bool variable_alignment_from_param() const;

public:

    bool variable_alignment() const;

    void variable_alignment(bool b);

            data_partition get_data_partition(int i) const;

    mutable_data_partition get_data_partition(int i);

    int n_data_partitions() const {return PC->DPC.size();}

            data_partition operator[](int i) const {return get_data_partition(i);}

    mutable_data_partition operator[](int i)       {return get_data_partition(i);}

    boost::optional<int> smodel_index_for_partition(int i) const {return PC->smodel_for_partition[i];}

    boost::optional<int> imodel_index_for_partition(int i) const {return PC->imodel_for_partition[i];}

    std::vector<int> partitions_for_imodel(int i) const;

    std::vector<int> partitions_for_scale(int i) const;

    boost::optional<int> scale_index_for_partition(int i) const {return PC->scale_for_partition[i];}
  
    void set_beta(double b);
    double get_beta() const;

    TreeInterface t() const;

    std::vector<std::string> get_labels() const;
  
    expression_ref my_tree() const;
    expression_ref my_imodels() const;
    expression_ref my_substitution_branch_lengths() const;

    void reconnect_branch(int s1, int t1, int t2);
    void prune_subtree(const tree_edge& b_subtree);
    void regraft_subtree(const tree_edge& b_subtree, const tree_edge& b_target);
    void begin_modify_topology();
    void end_modify_topology();

private:
    void exchange_subtrees(int b1, int b2);
public:
    void NNI(const tree_edge& b1, const tree_edge& b2, bool disconnect_subtree=false);
    void NNI(int b1, int b2, bool disconnect_subtree=false);
    void show_h_tree() const;

    /// Did we most recently hit beta==1 (1) or beta=0 (0)
    int updown;
    int branch_prior_type() const {return PC->branch_prior_type;}
    int branch_HMM_type(int b) const {
	assert(0 <= b and b < PC->branch_HMM_type.size());
	return PC->branch_HMM_type[b];
    }
  
    // The prior, likelihood, and probability
    log_double_t prior_no_alignment() const;
    log_double_t prior_alignment() const;

    /// How many substitution models?
    int n_smodels() const {return PC->SModels.size();}
  
    /// How many substitution models?
    int n_imodels() const {return PC->n_imodels;}

    /// How many substitution models?
    int n_branch_scales() const {return PC->n_scales;}

    object_ptr<const alphabet> get_alphabet_for_smodel(int i) const;

    /// Set branch 'b' to have length 'l', and compute the transition matrices
    void setlength_unsafe(int b, double l);
    void setlength(int b,double l); 
  
    void select_root(int b) const;
    void set_root(int n) const;
    int subst_root() const;

    double branch_mean() const;

    int branch_scale_index(int i) const;

    boost::optional<int> branch_scale_modifiable_reg(int i) const;

    double branch_scale(int i) const;

    void branch_scale(int i,double); 

    double get_branch_subst_rate(int p, int b) const;

    Parameters* clone() const {return new Parameters(*this);}

    std::string name() const {return "Full Model";}

    Parameters(const std::shared_ptr<module_loader>& L,
	       const std::vector<alignment>& A, const SequenceTree&, 
	       const std::vector<model_t>&,
	       const std::vector<boost::optional<int>>&,
	       const std::vector<model_t>&,
	       const std::vector<boost::optional<int>>&,
	       const model_t&,
	       const std::vector<int>&,
	       const key_map_t&);

    Parameters(const std::shared_ptr<module_loader>& L,
	       const std::vector<alignment>& A, const SequenceTree&, 
	       const std::vector<model_t>&,
	       const std::vector<boost::optional<int>>&,
	       const std::vector<model_t>&,
	       const std::vector<boost::optional<int>>&,
	       const std::vector<model_t>&,
	       const std::vector<boost::optional<int>>&,
	       const model_t&,
	       const std::vector<int>&,
	       const key_map_t&);
};

bool accept_MH(const Model& P1, const Model& P2, log_double_t rho);
std::string write_newick(const Parameters&, bool=true);

#endif
