#ifndef SIMPLIFIER_H
#define SIMPLIFIER_H

#include <set>
#include <map>
#include <utility>
#include "computation/expression/expression_ref.H"
#include "computation/expression/let.H"
#include "computation/expression/var.H"

struct simplifier_options
{
    bool pre_inline_unconditionally = true;
    bool post_inline_unconditionally = true;
    bool let_float_from_case = true;
    bool let_float_from_apply = true;
    bool let_float_from_let = true;
    bool case_of_constant = true;
    bool case_of_variable = true;
    bool case_of_case = true;
    bool beta_reduction = true;
    int max_iterations = 4;
    int inline_threshhold = 8;
    double keenness = 1.5;

    bool fully_lazy = true;

    bool dump_parsed = false;
    bool dump_renamed = false;
    bool dump_desugared = false;
};

class module_loader;

typedef std::pair<expression_ref,occurrence_info> bound_variable_info;

typedef std::map<var, bound_variable_info> in_scope_set;

int get_n_lambdas1(const expression_ref& E);

expression_ref peel_n_lambdas1(const expression_ref& E, int n);

bool is_trivial(const expression_ref& E);

std::vector<CDecls> simplify_module(const simplifier_options&, const std::map<var, expression_ref>&, const std::set<var>&, const std::vector<CDecls>& decl_groups);

#endif
