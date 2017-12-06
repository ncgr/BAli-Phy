#include "translate.H"
#include "parse.H"
#include <vector>
#include <set>
#include "rules.H"
#include "myexception.H"
#include "util.H"
#include <iostream>
#include "setup.H"

using std::vector;
using std::set;
using std::pair;
using std::map;
using std::string;
using boost::property_tree::ptree;
using boost::optional;

// Translate pass M+FM -> RCTMC[M,FM]
void pass1(const Rules& R, ptree& p)
{
    // 1. Handle children.
    for(auto& child: p)
	pass1(R, child.second);
    
    // 2. Convert e.g. TN+F -> RCTMC[TN,F]
    if (unify(R.get_result_type(p),parse_type("FrequencyModel[_]")) and p.count("submodel"))
    {
	ptree q = p.get_child("submodel");
	auto& r = p;
	r.erase("submodel");
	
	ptree result = {};
	result.put_value("RCTMC");
	result.push_back({"Q",q});
	result.push_back({"R",r});
	p = result;
    }
}

/// True if some conversion function can be applied to the expression of type t1, so that it is of type t2
equations convertible_to(ptree& model, const type_t& t1, type_t t2)
{
    auto E = unify(t1, t2);
    if (E)
	return E;

    if (t2.get_value<string>() == "Double")
    {
	t2.put_value("Int");
	E = convertible_to(model, t1, t2);
	if (E)
	{
	    ptree result;
	    result.put_value("intToDouble");
	    result.push_back({"x",model});
	    model = result;
	}
    }
    else if (t2.get_value<string>() == "MultiMixtureModel")
    {
	t2.put_value("MixtureModel");
	E = convertible_to(model,t1,t2);
	if (E)
	{
	    ptree result;
	    result.put_value("MultiMixtureModel");
	    result.push_back({"submodel",model});
	    model = result;
	}
    }
    else if (t2.get_value<string>() == "MixtureModel")
    {
	t2.put_value("RA");
	E = convertible_to(model,t1,t2);
	if (E)
	{
	    ptree result;
	    result.put_value("UnitMixture");
	    result.push_back({"submodel",model});
	    model = result;
	}
    }
    else if (t2.get_value<string>() == "RA")
    {
	t2.put_value("ExchangeModel");
	E = convertible_to(model,t1,t2);
	if (E)
	{
	    ptree result;
	    result.put_value("RCTMC");
	    result.push_back({"Q",model});
	    result.push_back({"R",ptree("F")});
	    model = result;
	}
    }

    return E;
}

set<string> find_rule_type_vars(const ptree& rule)
{
    set<string> vars = find_variables_in_type(rule.get_child("result_type"));
    for(const auto& x: rule.get_child("args"))
	add(vars, find_variables_in_type( x.second.get_child("arg_type") ) );
    return vars;
}

Rule substitute_in_rule_types(const map<string,term_t>& renaming, Rule rule)
{
    substitute(renaming, rule.get_child("result_type") );
    substitute(renaming, rule.get_child("constraints") );
    for(auto& x: rule.get_child("args"))
    {
	ptree& arg_type = x.second.get_child("arg_type");
	substitute( renaming, arg_type );
    }
    return rule;
}

Rule substitute_in_rule_types(const equations& renaming, Rule rule)
{
    substitute(renaming, rule.get_child("result_type") );
    substitute(renaming, rule.get_child("constraints") );
    for(auto& x: rule.get_child("args"))
    {
	ptree& arg_type = x.second.get_child("arg_type");
	substitute( renaming, arg_type );
    }
    return rule;
}

Rule freshen_type_vars(Rule rule, const set<string>& bound_vars)
{
    // 1. Find variables in rule type
    set<string> rule_type_variables = find_rule_type_vars(rule);

    // 2. Make substitutions in rule type
    auto renaming = alpha_rename(rule_type_variables, bound_vars);
    return substitute_in_rule_types(renaming, rule);
}


set<string> find_type_variables_from_scope(const vector<pair<string, ptree>>& scope)
{
    set<string> vars;
    for(auto& x: scope)
	add(vars, find_variables_in_type(x.second));
    return vars;
}

optional<ptree> type_for_var_in_scope(const string& name, const vector<pair<string,ptree>>& in_scope)
{
    for(int i=int(in_scope.size())-1;i>=0;i--)
	if (in_scope[i].first == name)
	    return in_scope[i].second;
    return boost::none;
}


vector<pair<string, ptree>> extend_scope(const ptree& rule, int skip, const vector<pair<string, ptree>>& scope)
{
    auto scope2 = scope;
    int i=0;
    for(const auto& arg: rule.get_child("args"))
    {
	i++;
	if (i < skip) continue;

	const auto& argument = arg.second;
	string arg_name = argument.get<string>("arg_name");
	type_t arg_required_type = get_type_for_arg(rule, arg_name);

	scope2.push_back({arg_name, arg_required_type});
    }

    return scope2;
}


typedef pair<vector<type_t>, type_t> instance_def;
typedef pair<vector<type_t>, type_t> class_def;
typedef pair<vector<type_t>, type_t> data_def;

vector<instance_def> get_instances()
{
    vector<instance_def> instances;
    instances.push_back( {{}, parse_type("Alphabet[Numeric]")}     );
    instances.push_back( {{}, parse_type("Alphabet[DNA]")}         );
    instances.push_back( {{}, parse_type("Alphabet[RNA]")}         );
    instances.push_back( {{}, parse_type("Alphabet[AA]")}          );
    instances.push_back( {{}, parse_type("Alphabet[Triplets[a]]")} );
    instances.push_back( {{}, parse_type("Alphabet[Codons[a]]")}   );

    instances.push_back( {{}, parse_type("Nucleotides[DNA]")} );
    instances.push_back( {{}, parse_type("Nucleotides[RNA]")} );

    instances.push_back( {{}, parse_type("Triplets[Codons[a],a]")} );
    instances.push_back( {{}, parse_type("Triplets[Triplets[a],a]")} );

    instances.push_back( {{}, parse_type("Num[Int]")}    );
    instances.push_back( {{}, parse_type("Num[Double]")} );

    return instances;
}

vector<class_def> get_classes()
{
    vector<class_def> classes;
    classes.push_back({{parse_type("Alphabet[a]"),parse_type("Alphabet[b]")},  parse_type("Triplets[a,b]")});
    classes.push_back({{},  parse_type("Alphabet[a]")});
    classes.push_back({{parse_type("Alphabet[a]")},  parse_type("Nucleotides[a]")});
    classes.push_back({{},  parse_type("Num[a]")});
    return classes;
}

vector<class_def> get_data()
{
    vector<class_def> data;
    data.push_back({{parse_type("Nucleotides[a]")},parse_type("Codons[a]")});
    data.push_back({{parse_type("Alphabet[a]")},parse_type("Triplets[a]")});
    data.push_back({{parse_type("Alphabet[a]")},parse_type("ExchangeModel[a]")});
    data.push_back({{parse_type("Alphabet[a]")},parse_type("RA[a]")});
    data.push_back({{parse_type("Alphabet[a]")},parse_type("FrequencyModel[a]")});
    data.push_back({{parse_type("Alphabet[a]")},parse_type("MixtureModel[a]")});
    data.push_back({{parse_type("Alphabet[a]")},parse_type("MultiMixtureModel[a]")});
    data.push_back({{},parse_type("IndelModel")});
    data.push_back({{},parse_type("AA")});
    data.push_back({{},parse_type("DNA")});
    data.push_back({{},parse_type("RNA")});
    return data;
}

auto instances = get_instances();
auto classes = get_classes();
auto data = get_data();

pair<vector<type_t>,type_t> freshen_vars(const pair<vector<type_t>,type_t>& x, const set<string>& bound_vars)
{
    set<string> type_vars;
    for(auto& y: x.first)
	add(type_vars, find_variables_in_type(y));
    add(type_vars, find_variables_in_type(x.second));
    auto renaming = alpha_rename(type_vars, bound_vars);

    auto x2 = x;
    for(auto& y: x2.first)
	substitute(renaming, y);
    substitute(renaming, x2.second);

    return x2;
}


bool is_pure_constraint(const type_t& t)
{
    for(auto& arg: t)
	if (not is_variable(arg.second)) return false;

    return true;
}

equations simplify_constraint(type_t constraint, const Rules& R)
{
//    std::cerr<<"Simplifying constraint: "<<unparse_type(constraint)<<std::endl;

    if (is_pure_constraint(constraint))
    {
	equations E;
	E.add_constraint(constraint);
	return E;
    }

    for(auto& instance: instances)
    {
	equations E;
	instance = freshen_vars(instance, find_variables_in_type(constraint));
//	std::cerr<<"   instance = "<<unparse_type(instance.second)<<std::endl;
	for(auto& sub_constraint: instance.first)
	    E.add_constraint(sub_constraint);
//	std::cerr<<"   proposed E: "<<show(E)<<std::endl;
	E = E && unify(constraint, instance.second);
	if (E)
	{
	    E.eliminate_except(find_variables_in_type(constraint));
//	    std::cerr<<"   successful E: "<<show(E)<<std::endl;
	    return E;
	}
    }
    throw myexception()<<"No instance for constraint '"<<unparse_type(constraint)<<"'!";
    return equations(false);
}

equations simplify_constraints(const equations& E, const Rules& R)
{
    equations E2 = E;
    E2.get_constraints().clear();

//    std::cerr<<"Simplifying constraints: "<<show(E2)<<std::endl;

    for(auto& constraint: E.get_constraints())
    {
	auto constraint2 = constraint;
	substitute(E, constraint2);

	E2 = E2 && simplify_constraint(constraint2, R);
    }
    return E2;
}

// OK, so 'model' is going to have arg=value pairs set, but not necessarily in the right order.
equations pass2(const Rules& R, const ptree& required_type, ptree& model, set<string> bound_vars, const vector<pair<string, ptree>>& scope)
{
    // 0a. Any variables in required_type must be listed as bound
    assert(includes(bound_vars, find_variables_in_type(required_type)));
    // 0b. Any type variables in scope must also be listed as bound
    assert(includes(bound_vars, find_type_variables_from_scope(scope)));

    auto name = model.get_value<string>();

    // 1. Get result type and the rule, if there is one.
    type_t result_type;
    optional<Rule> rule;
    if (can_be_converted_to<int>(name))
	result_type=type_t("Int");
    else if (can_be_converted_to<double>(name))
	result_type=type_t("Double");
    else if (name.size()>=2 and name[0] == '"' and name.back() == '"' and required_type.get_value<string>() == "String") 
	result_type=type_t("String");
    else if (auto type = type_for_var_in_scope(name, scope))
	result_type = *type;
    else
    {
	rule = R.require_rule_for_func(name);
	rule = freshen_type_vars(*rule, bound_vars);

	//	std::cout<<"name = "<<name<<" required_type = "<<unparse_type(required_type)<<"  result_type = "<<unparse_type(result_type)<<std::endl;

	result_type = rule->get_child("result_type");
    }

    // 2. Unify required type with rule result type
    auto E = unify(result_type, required_type);
    if (rule)
	for(const auto& constraint: rule->get_child("constraints"))
	    E.add_constraint(constraint.second);

//    std::cerr<<show(E)<<std::endl;
    
    // 3. Attempt a conversion if the result_type and the required_type don't match.
    if (not E)
    {
	if (convertible_to(model, result_type, required_type))
	    return pass2(R, required_type, model, bound_vars, scope);
	else
	    throw myexception()<<"Term '"<<model.get_value<string>()<<"' of type '"<<unparse_type(result_type)
			       <<"' cannot be converted to type '"<<unparse_type(required_type)<<"'";
    }

    // 4. If this is a constant or variable, then we are done here.
    if (not rule)
    {
	if (not model.empty())
	    throw myexception()<<"Term '"<<model.get_value<string>()<<"' of type '"<<unparse_type(result_type)
			       <<"' should not have arguments!";
	return E;
    }

    // We need to be able to check constraints here, but we might have Num[a] => a=ExchangeModel[a]
    try
    {
	E = E && simplify_constraints(E, R);
    }
    catch (myexception& e)
    {
	std::ostringstream o;
	o<<"Parsing: '"<<unparse(model)<<"':\n  ";
	e.prepend(o.str());
	throw e;
    }

    // 5.1 Update required type and rules with discovered constraints
    rule = substitute_in_rule_types(E, *rule);

    // 5.2 Record any new variables that we are using as bound variables
    //     (I think that only rules can introduce new variables)
    add(bound_vars, find_rule_type_vars(*rule));
    
    // 5.3. Check type of arguments if pass_arguments
    if (rule->get("pass_arguments",false) or rule->get("list_arguments",false))
    {
	for(auto& child: model)
	{
	    type_t arg_required_type = get_type_for_arg(*rule, "*");
	    E = E && pass2(R, arg_required_type, child.second, bound_vars, scope);
	    if (not E)
	    {
		auto x = arg_required_type;
		substitute(E,x);
		throw myexception()<<"Expression '"<<unparse(child.second)<<"' is not of required type "<<unparse_type(x)<<"!";
	    }
	    E = E && simplify_constraints(E, R);
	    add(bound_vars, E.referenced_vars());
	}
	E.eliminate_except(find_variables_in_type(required_type));
	return E;
    }

    // 6. Complain if undescribed arguments are supplied
    for(auto& child: model)
	if (get_arg(*rule, child.first).get("no_apply",false))
	    throw myexception()<<"Rule for function '"<<rule->get<string>("name")<<"' doesn't allow specifying a value for '"<<child.first<<"'.";

    // Create the new model tree with args in correct order
    ptree model2(name);

    // 7. Handle arguments in rule order
    int skip=0;
    for(const auto& arg: rule->get_child("args"))
    {
	skip++;
	const auto& argument = arg.second;

	string arg_name = argument.get<string>("arg_name");
	bool arg_is_required = (arg_name != "*") and not argument.get("no_apply",false);

	// If this is a supplied argument
	if (model.count(arg_name))
	{
	    type_t arg_required_type = get_type_for_arg(*rule, arg_name);
	    substitute(E, arg_required_type);
	    ptree arg_value = model.get_child(arg_name);
	    E = E && pass2(R, arg_required_type, arg_value, bound_vars, extend_scope(*rule,skip,scope));
	    if (not E)
		throw myexception()<<"Expression '"<<unparse(arg_value)<<"' is not of required type "<<unparse_type(arg_required_type)<<"!";
	    E = E && simplify_constraints(E, R);
	    model2.push_back({arg_name, arg_value});
	    add(bound_vars, E.referenced_vars());
	}

	// If there's no default arg 
	else if (not argument.count("default_value"))
	{
	    if (arg_is_required)
		throw myexception()<<"Command '"<<name<<"' missing required argument '"<<arg_name<<"'";
	}
	else
	{
	    auto arg_required_type = argument.get_child("arg_type");
	    substitute(E, arg_required_type);
	    auto default_arg = argument.get_child("default_value");
	    E = E && pass2(R, arg_required_type, default_arg, bound_vars, extend_scope(*rule, skip, scope));
	    if (not E)
		throw myexception()<<"Expression '"<<unparse(default_arg)<<"' is not of required type "<<unparse_type(arg_required_type)<<"!";
	    E = E && simplify_constraints(E, R);
	    add(bound_vars, E.referenced_vars());

	    model2.push_back({arg_name, default_arg});
	}
    }

    model = model2;

    auto keep = find_variables_in_type(required_type);
    add(keep, find_type_variables_from_scope(scope));
    E.eliminate_except(keep);

    return E;
}

std::pair<ptree,equations> translate_model(const Rules& R, const ptree& required_type, ptree model, const vector<pair<string,ptree>>& scope)
{
    pass1(R, model);
    auto E = pass2(R, required_type, model, find_variables_in_type(required_type), scope);
    return {model,E};
}

