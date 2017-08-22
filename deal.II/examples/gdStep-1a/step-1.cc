/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

#include "step-1.h"

void step1_grid ()
{
  Triangulation<2> divideIt;
  GridGenerator::hyper_cube (divideIt);
  divideIt.refine_global (5);
  std::ofstream out ("step-1.eps");
  GridOut grid_out;
  grid_out.write_eps (divideIt, out);
  std::cout << "Grid written to step-1.eps" << std::endl;
}

int main (void)
{
  step1_grid ();
  return 0;
}
