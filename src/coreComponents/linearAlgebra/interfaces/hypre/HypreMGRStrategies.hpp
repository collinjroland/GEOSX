/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreMGRStrategies.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_

#include "common/DataTypes.hpp"

#include "HyprePreconditioner.hpp"

//#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
//#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
//#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
//
#include <_hypre_utilities.h>
//#include <_hypre_parcsr_ls.h>
//#include <_hypre_IJ_mv.h>
//#include <krylov.h>

namespace geosx
{

/**
 * @todo Add a detailed description with an example
 *
 * @brief Compute an array of unique component labels.
 *
 * @details
 *
 * @param numComponentsPerField array of number of components per field
 * @param numLocalDofsPerField array of local number of dofs per field
 * @return array1d of HYPRE_Int labels
 */
inline array1d< HYPRE_Int > computeLocalDofComponentLabels( arraySlice1d< localIndex const > const & numComponentsPerField,
                                                            arraySlice1d< localIndex const > const & numLocalDofsPerField )
{
  array1d< HYPRE_Int > ret;
  HYPRE_Int numFields = LvArray::integerConversion< HYPRE_Int >( numLocalDofsPerField.size() );

  if( numFields > 0 )
  {
    HYPRE_Int numTotalLocalDof = 0;
    for( HYPRE_Int i = 0; i < numFields; ++i )
    {
      numTotalLocalDof += LvArray::integerConversion< HYPRE_Int >( numLocalDofsPerField[i] );
    }

    ret.resize( numTotalLocalDof );

    HYPRE_Int firstLabel = 0;
    HYPRE_Int istr= 0;
    HYPRE_Int iend;
    HYPRE_Int numComp;
    for( HYPRE_Int iFld = 0; iFld < numFields; ++iFld )
    {
      numComp = LvArray::integerConversion< HYPRE_Int > ( numComponentsPerField[iFld] );
      array1d< HYPRE_Int > vectorLabels( numComp );
      for ( HYPRE_Int k = 0; k < numComp; ++k )
      {
        vectorLabels[k] = k + firstLabel;
      }
      iend = istr + LvArray::integerConversion< HYPRE_Int >( numLocalDofsPerField[iFld] );;
      for( localIndex i = istr; i < iend; i += numComp )
      {
        for ( integer k = 0; k < numComp; ++k )
        {
          ret[i+k] = vectorLabels[k];
        }
      }
      istr += iend;
      firstLabel += numComp;
    }
  }
  return ret;
}
//
///**
// * @brief Convert GEOSX global index value to hypre bigint
// * @param index the input value
// * @return the converted value
// */
//inline HYPRE_BigInt toHYPRE_BigInt( globalIndex const index )
//{
//  return LvArray::integerConversion< HYPRE_BigInt >( index );
//}
//
///**
// * @brief Converts a non-const array from GEOSX globalIndex type to HYPRE_BigInt
// * @param[in] index the input array
// * @return the converted array
// */
//inline HYPRE_BigInt * toHYPRE_BigInt( globalIndex * const index )
//{
//  return reinterpret_cast< HYPRE_BigInt * >(index);
//}
//
///**
// * @brief Converts a const array from GEOSX globalIndex type to HYPRE_BigInt
// * @param[in] index the input array
// * @return the converted array
// */
//inline HYPRE_BigInt const * toHYPRE_BigInt( globalIndex const * const index )
//{
//  return reinterpret_cast< HYPRE_BigInt const * >(index);
//}
//
///**
// * @brief Container for hypre preconditioner function pointers.
// */
//struct HyprePrecFuncs
//{
//  /// Alias for setup function type
//  using SetupFunc = HYPRE_Int (*)( HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector );
//
//  /// Alias for apply function type
//  using ApplyFunc = HYPRE_Int (*)( HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector );
//
//  /// Alias for destroy function type
//  using DestroyFunc = HYPRE_Int (*)( HYPRE_Solver );
//
//  SetupFunc setup{};     ///< pointer to setup function
//  ApplyFunc apply{};     ///< pointer to apply function
//  DestroyFunc destroy{}; ///< pointer to destroy function
//};
//
///**
// * @brief Container for hypre Krylov solver function pointers.
// */
//struct HypreSolverFuncs
//{
//  /// Alias for set preconditioner function type
//  using SetPrecondFunc = HYPRE_Int ( * )( HYPRE_Solver,
//                                          HYPRE_PtrToParSolverFcn,
//                                          HYPRE_PtrToParSolverFcn,
//                                          HYPRE_Solver );
//
//  /// Alias for setup function type
//  using SetupFunc = HYPRE_Int ( * )( HYPRE_Solver,
//                                     HYPRE_ParCSRMatrix,
//                                     HYPRE_ParVector,
//                                     HYPRE_ParVector );
//
//  /// Alias for solve function type
//  using SolveFunc = HYPRE_Int ( * )( HYPRE_Solver,
//                                     HYPRE_ParCSRMatrix,
//                                     HYPRE_ParVector,
//                                     HYPRE_ParVector );
//
//  /// Alias for get number of iterations function type
//  using GetNumIter = HYPRE_Int ( * )( HYPRE_Solver solver,
//                                      HYPRE_Int * num_iterations );
//
//  /// Alias for get final residual norm function type
//  using GetFinalNorm = HYPRE_Int ( * )( HYPRE_Solver solver,
//                                        HYPRE_Real * norm );
//
//  /// Alias for destroy function type
//  using DestroyFunc = HYPRE_Int ( * )( HYPRE_Solver );
//
//  SetPrecondFunc setPrecond{}; ///< pointer to set preconditioner function
//  SetupFunc setup{};           ///< pointer to setup function
//  SolveFunc solve{};           ///< pointer to solve function
//  GetNumIter getNumIter{};     ///< pointer to get number of iterations function
//  GetFinalNorm getFinalNorm{}; ///< pointer to get final residual norm function
//  DestroyFunc destroy{};       ///< pointer to destroy function
//};

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_*/
