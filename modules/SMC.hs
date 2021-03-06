module SMC where

import Distributions
builtin builtin_smc_density 4 "smc_density" "SMC"
smc_density rho_over_theta rates level_boundaries sequences = builtin_smc_density rho_over_theta rates' level_boundaries' sequences
                                                                where rates' = list_to_vector rates
                                                                      level_boundaries' = list_to_vector level_boundaries
smc rho_over_theta rates level_boundaries = ProbDensity (smc_density rho_over_theta rates level_boundaries) (error "SMC has no quantile") (return 0) (listToString("AlignmentRangeString"))

