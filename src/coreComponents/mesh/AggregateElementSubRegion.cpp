/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include "AggregateElementSubRegion.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include <mpi.h>

namespace geosx
{
AggregateElementSubRegion::AggregateElementSubRegion( string const & name,
                                                      dataRepository::ManagedGroup * const parent):
  ElementSubRegionBase( name, parent )
{
  RegisterViewWrapper(viewKeyStruct::elementCenterString, &m_elementCenter, 0 );
  RegisterViewWrapper(viewKeyStruct::elementVolumeString, &m_elementVolume, 0 );
}

AggregateElementSubRegion::~AggregateElementSubRegion()
{
}

void AggregateElementSubRegion::CreateFromFineToCoarseMap( localIndex nbAggregates,
                                                           array1d< localIndex > const & fineToCoarse,
                                                           array1d< R1Tensor > const & barycenters,
                                                           array1d< real64 > const & volumes )
{
  this->resize( nbAggregates );
  m_elementCenter = barycenters;
  m_elementVolume = volumes;
  m_nbFineCellsPerCoarseCell.resize( nbAggregates + 1 );
  m_fineToCoarse = fineToCoarse;
  m_fineByAggregates.resize( fineToCoarse.size() ) ;
  
  /// First loop to count the number of fine cells per coarse cell
  for( localIndex fineCell = 0; fineCell < fineToCoarse.size(); fineCell++ )
  {
    m_nbFineCellsPerCoarseCell[1 + fineToCoarse[fineCell]]++;
  }

  /// Second loop to cumulate the number of fine cells
  for( localIndex coarseCell = 1; coarseCell <= nbAggregates; coarseCell++ )
  {
    m_nbFineCellsPerCoarseCell[coarseCell] += m_nbFineCellsPerCoarseCell[coarseCell-1]; 
  }

  /// Third loop to order the index of the fine cells
  array1d< localIndex > offset( nbAggregates );
  for( localIndex fineCell = 0; fineCell < fineToCoarse.size(); fineCell++ )
  {
    localIndex coarseCell = fineToCoarse[fineCell];
    m_fineByAggregates[m_nbFineCellsPerCoarseCell[coarseCell] + offset[coarseCell]++] = fineCell;
  }

  int mpiSize;
  int mpiRank;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );

  array1d< localIndex > nbAggregatesPerRank(mpiSize);
  MPI_Allgather( &nbAggregates, 1, MPI_LONG_LONG,nbAggregatesPerRank.data(), 1, MPI_LONG_LONG, MPI_COMM_GEOSX);
  globalIndex offsetForGlobalIndex = 0;
  for( localIndex i = 0; i < mpiRank; i++)
  {
    offsetForGlobalIndex += nbAggregatesPerRank[i];
  }
  GEOS_LOG_RANK("offset of the subregion " << offsetForGlobalIndex);
  GEOS_LOG_RANK("size local2global " << this->m_localToGlobalMap.size());
  GEOS_LOG_RANK("size global2local " << this->m_globalToLocalMap.size());
  for( localIndex i = 0; i < nbAggregates; i++)
  {
    this->m_localToGlobalMap[i] = i + offsetForGlobalIndex;
  }

}

void AggregateElementSubRegion::ComputeGhosts()
{
   ElementRegion const * elementRegion = this->getParent()->getParent()->group_cast<ElementRegion const *>();
   std::cout << "Compute ghost "<< std::endl;
       GEOS_LOG_RANK( "in compute ghost size "<< this->size());
   elementRegion->forElementSubRegions([&]( auto const * const subRegion )-> void
   {
     std::cout << "sub region size " << subRegion->size() << std::endl;
     std::cout << "sub region ghost element " << subRegion->GetNumberOfGhosts() << std::endl;
     localIndex oldSize = m_fineToCoarse.size();
     m_fineToCoarse.resize(  oldSize  + subRegion->GetNumberOfGhosts());
      auto & aggregateIndex =
        subRegion->template getWrapper< array1d< localIndex > > ("aggregateIndex")->reference();
       GEOS_LOG_RANK(subRegion->m_globalToLocalMap.size());
       GEOS_LOG_RANK(subRegion->m_localToGlobalMap.size());
       GEOS_LOG_RANK("aggregateIndex ssize " << aggregateIndex.size());
       GEOS_LOG_RANK("oldsize " <<oldSize);
     for(localIndex i = oldSize ; i < m_fineToCoarse.size();i++)
     {
       m_fineToCoarse[i] = aggregateIndex[i];
     }
   });
   GEOS_LOG_RANK("new fine to coarse size" << m_fineToCoarse.size());
   GEOS_LOG_RANK("new local to glocal size " << m_localToGlobalMap.size());
   GEOS_LOG_RANK("new global to local size " << m_globalToLocalMap.size());
   GEOS_ERROR_IF(true,"staph");
}
}
