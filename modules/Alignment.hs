module Alignment where
{
  import Tree;
  builtin pairwise_alignment_length1 1 "pairwise_alignment_length1" "Alignment";
  builtin pairwise_alignment_length2 1 "pairwise_alignment_length2" "Alignment";
  builtin transition_counts 1 "transition_counts" "Alignment";
  builtin pairwise_alignment_probability_from_counts 2 "pairwise_alignment_probability_from_counts" "Alignment";
  builtin builtin_load_alignment 1 "load_alignment" "Alignment";
  
  builtin merge_alignment_constraints  4 "merge_alignment_constraints" "Alignment";
  builtin leaf_alignment_constraint    4 "leaf_alignment_constraint"   "Alignment";

  branch_hmms (model,_) distances n_branches = listArray' $ map (model distances) [0..2*n_branches-1];
  
  alignment_branch_pr a hmms b = pairwise_alignment_probability_from_counts (transition_counts (a!b)) (hmms!b);
  
  seqlength a tree node = pairwise_alignment_length1 (a!b) where { b = head $ edgesOutOfNode tree node};
  
  product' = foldl' (*) (doubleToLogDouble 1.0);
  alignment_pr_top a tree hmm = product' $ map (alignment_branch_pr a hmm) [0..numBranches tree - 1];
  alignment_pr_bot a tree (_,lengthp) = (product' $ map (lengthp . seqlength a tree) (internal_nodes tree))^2;
  alignment_pr a tree hmm model = (alignment_pr_top a tree hmm)/(alignment_pr_bot a tree model);
  alignment_pr1 seq (_,lengthp) = lengthp (sizeOfVectorInt seq);
  load_alignment filename = builtin_load_alignment (listToString filename);

  alignment_constraints t m delta as seqs= let {constraints = mkArray (2*numBranches t) con_func;
                                                con_func b = case edgesBeforeEdge t b of
                                                               {
                                                                 [] -> let {n=sourceNode t b} in
                                                                       leaf_alignment_constraint m delta n (seqs!n);
                                                                 [b1,b2] -> merge_alignment_constraints (constraints!b1) (as!b1) (constraints!b2) (as!b2)
                                                               }
                                               }
                                           in constraints;
}  
