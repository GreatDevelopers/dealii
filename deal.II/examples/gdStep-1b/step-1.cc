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

unsigned int No_of_divisionC = 0; // Divisions along Circumference; 0
  // mean code will try to keep aspect ratio of elements to minimum.
  
unsigned int No_of_divisionR = 2; // Divisions along Radial
//unsigned int No_of_divisionR = No_of_divisionC;

void step1_grid ()
{
  const Point<2> center (1,0);
  const double inner_radius = 0.5,
               outer_radius = 1.0;

  const HyperShellBoundary<2> boundary_description(center);

  Triangulation<2> divideIt;

  divideIt.set_boundary (0, boundary_description);

  GridGenerator::hyper_shell (divideIt,
                              center, inner_radius, outer_radius,
                              No_of_divisionC, false);

  for (unsigned int step=0; step < No_of_divisionR; ++step)
    {
      Triangulation<2>::active_cell_iterator
      cell = divideIt.begin_active(),
      endc = divideIt.end();

      for (; cell!=endc; ++cell)
 
         for (unsigned int v=0;
             v < GeometryInfo<2>::vertices_per_cell;
             ++v)
            {
             const double distance_from_center
              = center.distance (cell->vertex(v));

            if (std::fabs(distance_from_center - inner_radius) 
               < TOLRANCE)
              {
                cell->set_refine_flag ();
                break;
              }
          }
      divideIt.execute_coarsening_and_refinement ();
    }

  std::ofstream out ("step-1.eps");
  GridOut grid_out;
  grid_out.write_eps (divideIt, out);

//  divideIt.set_boundary (0);
  std::cout << "Grid written to step-1.eps" << std::endl;
 }

int main (void)
{
  step1_grid ();
  return 0;
}
