/*
   Copyright (C) 2005-2007, 2012-13, 2017 Benjamin Redelings

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

#ifndef POW2_H
#define POW2_H

#include <cmath>
#include <cstdlib>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef HAVE_LOG2
#ifndef log2
#ifdef M_LOG2E
inline double log2(double x)
{
  return log(x)*M_LOG2E;
}
#elif M_LOG2_E
inline double log2(double x)
{
  return log(x)*M_LOG2_E;
}
#else
inline double log2(double x)
{
  constexpr double temp = 1.0/log(2);
  return log(x)*temp;
}
#endif
#endif
#endif

namespace fp_scale {

  constexpr int max = 50;
  constexpr int shift = max;

  extern double table[max*2+1];

  inline double pow2(int i) 
  {
    // NOTE: we can't use abs(i), because abs(INT_MIN) == INT_MIN < 0
    if (i < -max or i > max)
      return exp2(i);
    else
      return table[i+shift];
  }

  void initialize();

  // 11 bits for exponent in double precision
  // 10 bits for negative exponents: exponent range in [0,1024)
  // 9 bits usable if we keep at half the exponent range: [0,512) 
  constexpr double hi_cutoff = 1.0e154;  // 2**-512 == 10**-154
  constexpr double lo_cutoff = 1.0/hi_cutoff;  // 2**-512 == 10**-154
}

using fp_scale::pow2;

#endif
