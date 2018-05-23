#ifndef ALIGNMENT_CONSTRAINT_H
#define ALIGNMENT_CONSTRAINT_H
/*
   Copyright (C) 2004-2007,2009 Benjamin Redelings

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

#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>
#include <tuple>
#include <utility>
#include "alignment/alignment.H"
#include "tree/sequencetree.H"
#include "dp/2way.H"
#include "io.H"

typedef std::tuple<int, boost::optional<int>, boost::optional<int>, int> alignment_constraint;
typedef Vector<alignment_constraint> alignment_constraints;

std::vector<boost::optional<int>> get_max_y_le_x(const pairwise_alignment_t& a_xy);
std::vector<boost::optional<int>> get_min_y_ge_x(const pairwise_alignment_t& a_xy);
boost::optional<int> lookup(const std::vector<boost::optional<int>>& array, const boost::optional<int>& index);
alignment_constraints merge_alignment_constraints(const alignment_constraints& con_x, const pairwise_alignment_t& a_xz,
						  const alignment_constraints& con_y, const pairwise_alignment_t& a_yz,
						  const std::pair<matrix<int>,std::vector<int>>& totals);

std::vector<int> count_constrained_characters(const matrix<int>& m);
std::pair<matrix<int>,std::vector<int>> constraint_matrix_from_alignment(const alignment& A, int n);

std::vector<std::pair<int,int>> get_yboundaries_from_cons(int I, int J, const alignment_constraints& x_con, const alignment_constraints& y_con);

void check_constraints(const alignment_constraints& con_x, const alignment_constraints& con_y,
		       const pairwise_alignment_t& a_xy, const std::pair<matrix<int>, std::vector<int>>& con_matrix);

//---

matrix<int> load_alignment_constraint(const std::string& filename,SequenceTree& T);

std::vector< std::vector<int> > get_pins(const matrix<int>& constraint,const alignment& A,
					 const boost::dynamic_bitset<>& group1,const boost::dynamic_bitset<>& group2,
					 const std::vector<int>& seq1,const std::vector<int>& seq2);

boost::dynamic_bitset<> constraint_satisfied(const matrix<int>& constraint,const alignment& A);

void report_constraints(const boost::dynamic_bitset<>& s1, const boost::dynamic_bitset<>& s2, int p);

bool any_branches_constrained(const std::vector<int>&,const SequenceTree& T,const SequenceTree& TC, const std::vector<int>&);

std::vector< std::pair<int,int> > get_x_ranges_for_band(int D, 
							const std::vector<int>& seq1, 
							const std::vector<int>& seq2, 
							const std::vector<int>& seq12);
std::vector< std::pair<int,int> > get_y_ranges_for_band(int D, 
							const std::vector<int>& seq1, 
							const std::vector<int>& seq2, 
							const std::vector<int>& seq12);

std::vector< std::pair<int,int> > get_yboundaries_from_pins(int I, int J, const std::vector< std::vector<int> >& pins);

std::vector< std::pair<int,int> > boundaries_intersection(const std::vector< std::pair<int,int> >&,const std::vector< std::pair<int,int> >&);

#endif
