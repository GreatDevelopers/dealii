//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/trilinos_precondition_amg.h>

#ifdef DEAL_II_USE_TRILINOS

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <ml_include.h>
#include <ml_MultiLevelPreconditioner.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  PreconditionAMG::PreconditionAMG () 
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
             :
             communicator (MPI_COMM_WORLD)
#endif
  {}



  void
  PreconditionAMG::
  initialize (const SparseMatrix                    &matrix,
	      const bool                             elliptic,
	      const bool                             higher_order_elements,
	      const double                           aggregation_threshold,
	      const std::vector<std::vector<bool> > &null_space,
	      const bool                             output_details)
  {
    const unsigned int n_rows = matrix.m();
    const unsigned int null_space_dimension = null_space.size();

				        // Build the AMG preconditioner.
    Teuchos::ParameterList parameter_list;
  
    if (elliptic)
      {
	ML_Epetra::SetDefaults("SA",parameter_list);
	parameter_list.set("smoother: type", "Chebyshev");
	parameter_list.set("smoother: sweeps", 4);
      }
    else
      {
	ML_Epetra::SetDefaults("NSSA",parameter_list);
	parameter_list.set("aggregation: type", "Uncoupled");
	parameter_list.set("aggregation: block scaling", true);
      }
  
    parameter_list.set("aggregation: threshold", aggregation_threshold);
    
    if (output_details)
      parameter_list.set("ML output", 10);
    else
      parameter_list.set("ML output", 0);
  
    if (higher_order_elements)
      parameter_list.set("aggregation: type", "MIS");

    const Epetra_Map * domain_map = &(matrix.matrix->DomainMap());
    
    Epetra_MultiVector null_space_modes (*domain_map, null_space_dimension);
  
    if (null_space_dimension > 1)
      {
	Assert (n_rows == null_space[0].size(),
		ExcDimensionMismatch(n_rows,
				     null_space[0].size()));
	Assert (n_rows == (unsigned int)null_space_modes.GlobalLength(),
		ExcDimensionMismatch(n_rows,
				     null_space_modes.GlobalLength()));

	const unsigned int my_size = domain_map->NumMyElements();
	Assert (my_size == (unsigned int)domain_map->MaxLID()+1,
		ExcDimensionMismatch (my_size, domain_map->MaxLID()+1));
	
				        // Reshape null space as a
				        // contiguous vector of
				        // doubles so that Trilinos
				        // can read from it.
	for (unsigned int d=0; d<null_space_dimension; ++d)
	  for (unsigned int row=0; row<my_size; ++row)
	    {
	      int global_row_id = domain_map->GID(row);
	      null_space_modes.ReplaceMyValue(row, d, 
				 (double)null_space[d][global_row_id]);
	    }
  
	parameter_list.set("null space: type", "pre-computed");
	parameter_list.set("null space: dimension", null_space_modes.NumVectors());
	parameter_list.set("null space: vectors", null_space_modes.Values());
      }

    multilevel_operator = Teuchos::rcp (new ML_Epetra::MultiLevelPreconditioner(
				      *matrix.matrix, parameter_list, true));

    if (output_details)
      multilevel_operator->PrintUnused(0);
  }



  void
  PreconditionAMG::
  initialize (const ::dealii::SparseMatrix<double>  &deal_ii_sparse_matrix,
	      const bool                             elliptic,
	      const bool                             higher_order_elements,
	      const double                           aggregation_threshold,
	      const std::vector<std::vector<bool> > &null_space,
	      const bool                             output_details)
  {
    const unsigned int n_rows = deal_ii_sparse_matrix.m();
  
				        // Init Epetra Matrix, avoid
				        // storing the nonzero
				        // elements.

    Map.reset (new Epetra_Map(n_rows, 0, communicator));

    Matrix.reset();
    Matrix = boost::shared_ptr<SparseMatrix> (new SparseMatrix());

    Matrix->reinit (*Map, deal_ii_sparse_matrix);
    Matrix->compress();

    initialize (*Matrix, elliptic, higher_order_elements, 
		aggregation_threshold, null_space, output_details);
  }
  
  
  void PreconditionAMG::
  reinit ()
  {
    multilevel_operator->ReComputePreconditioner();
  }



  void PreconditionAMG::vmult (VectorBase       &dst,
			       const VectorBase &src) const
  {
    Assert (dst.vector->Map().SameAs(multilevel_operator->OperatorRangeMap()),
	    ExcNonMatchingMaps("dst"));
    Assert (src.vector->Map().SameAs(multilevel_operator->OperatorDomainMap()),
	    ExcNonMatchingMaps("src"));
    
    const int ierr = multilevel_operator->ApplyInverse (*src.vector, *dst.vector);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



				        // For the implementation of
				        // the <code>vmult</code>
				        // function we note that
				        // invoking a call of the
				        // Trilinos preconditioner
				        // requires us to use Epetra
				        // vectors as well. It is
				        // faster to provide a view,
				        // i.e., feed Trilinos with a
				        // pointer to the data, so we
				        // avoid copying the content
				        // of the vectors during the
				        // iteration. In the
				        // declaration of the right
				        // hand side, we need to cast
				        // the source vector (that is
				        // <code>const</code> in all
				        // deal.II calls) to
				        // non-constant value, as this
				        // is the way Trilinos wants
				        // to have them.
  void PreconditionAMG::vmult (dealii::Vector<double>       &dst,
			       const dealii::Vector<double> &src) const
  {
    Assert (Map->SameAs(Matrix->matrix->RowMap()),
	    ExcNonMatchingMaps("dst"));
    Assert (Map->SameAs(Matrix->matrix->RowMap()),
	    ExcNonMatchingMaps("src"));
    
    Epetra_Vector LHS (View, multilevel_operator->OperatorDomainMap(),
		       dst.begin());
    Epetra_Vector RHS (View, multilevel_operator->OperatorRangeMap(),
		       const_cast<double*>(src.begin()));
  
    const int res = multilevel_operator->ApplyInverse (RHS, LHS);
  
    Assert (res == 0,
	    ExcMessage ("Trilinos AMG MultiLevel preconditioner returned "
			"with an error!"));
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
