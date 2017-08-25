/* ---------------------------------------------------------------------
 * @f$Id: @ref step_2 "step-2".cc 30526 2013-08-29 20:06:27Z felix.gruber @f$
 *
 * Copyright (C) 1999 - 2013 by the deal.II authors
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
 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */

#include "step-2.h"

unsigned int No_of_divisionC = 0; // Divisions along Circumference; 0
  // mean code will try to keep aspect ratio of elements to minimum.
  
unsigned int No_of_divisionR = 5; // Divisions along Radial
//unsigned int No_of_divisionR = No_of_divisionC;

void make_grid (Triangulation<2> &divideIt)
{
  const Point<2> center (1,0);
  const double inner_radius = 0.5,
               outer_radius = 1.0;
  GridGenerator::hyper_shell (divideIt,
                              center, inner_radius, outer_radius,
                              No_of_divisionC, false);
  static const HyperShellBoundary<2> boundary_description(center);
  divideIt.set_boundary (0, boundary_description);
  for (unsigned int step=0; step< No_of_divisionR; ++step)
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
               < TOLERANCE)
              {
                cell->set_refine_flag ();
                break;
              }
          }
      divideIt.execute_coarsening_and_refinement ();
    }
}
void distribute_dofs (DoFHandler<2> &dof_handler)
{
  static const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs (finite_element);
  CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs(),
                                                        dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from (compressed_sparsity_pattern);
  std::ofstream out ("sparsity_pattern.1");
  sparsity_pattern.print_gnuplot (out);
}
void renumber_dofs (DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee (dof_handler);
  CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs(),
                                                        dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from (compressed_sparsity_pattern);
  std::ofstream out ("sparsity_pattern.2");
  sparsity_pattern.print_gnuplot (out);
}
int main ()
{
  Triangulation<2> divideIt;
  make_grid (divideIt);
  DoFHandler<2> dof_handler (divideIt);
  distribute_dofs (dof_handler);
  renumber_dofs (dof_handler);
}
