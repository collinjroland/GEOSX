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
  m_fineCellCenters.resize( fineToCoarse.size() );
  m_fineCellVolumes.resize( fineToCoarse.size() );
  
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
     GEOS_LOG_RANK("gloval index "<< i + offsetForGlobalIndex);
    this->m_globalToLocalMap[i + offsetForGlobalIndex] = i;
  }

}

struct newAggregate
{
  real64 volume;
  R1Tensor center;
  localIndex rankFrom;
  localIndex indexFrom;
  localIndex newIndex;
  array1d< localIndex> fineGhostCells;
  bool operator==(const newAggregate& rhs) const
  {
    return rhs.indexFrom == indexFrom;
  }
  bool operator<(const newAggregate& rhs) const
  {
    return indexFrom < rhs.indexFrom;
  }

};
void AggregateElementSubRegion::ComputeGhosts()
{
   ElementRegion const * elementRegion = this->getParent()->getParent()->group_cast<ElementRegion const *>();
   std::cout << "Compute ghost "<< std::endl;
       GEOS_LOG_RANK( "in compute ghost size "<< this->size());
   std::vector< newAggregate > ghostAggregates;
   std::unordered_map<localIndex, localIndex> mapNewAggregate;
   localIndex offset = this->size();
   localIndex oldSize = this->size();
   elementRegion->forElementSubRegions([&]( auto const * const subRegion )-> void
   {
     GEOS_LOG_RANK_0("global to local lsize " << subRegion->m_globalToLocalMap.size());
     GEOS_LOG_RANK_0("local to global " << subRegion->m_localToGlobalMap.size());
     std::cout << "sub region size " << subRegion->size() << std::endl;
     std::cout << "sub region ghost element " << subRegion->GetNumberOfGhosts() << std::endl;
     localIndex oldFineToCoarseSize = m_fineToCoarse.size();
     m_fineToCoarse.resize(  oldFineToCoarseSize  + subRegion->GetNumberOfGhosts());
      auto & aggregateIndex =
        subRegion->template getWrapper< array1d< localIndex > > ("aggregateIndex")->reference();
      auto & aggregateCenter =
        subRegion->template getWrapper< array1d< R1Tensor > > ("aggregateCenter")->reference();
      auto & aggregateVolume =
        subRegion->template getWrapper< array1d< real64 > > ("aggregateVolume")->reference();
       GEOS_LOG_RANK(subRegion->m_globalToLocalMap.size());
       GEOS_LOG_RANK(subRegion->m_localToGlobalMap.size());
       GEOS_LOG_RANK("aggregateIndex ssize " << aggregateIndex.size());
       GEOS_LOG_RANK("oldsize " <<oldSize);
       GEOS_LOG_RANK("old fine to coarse size " <<oldFineToCoarseSize);
       GEOS_LOG_RANK("m_nbFineCellsPerCoarseCell " << m_nbFineCellsPerCoarseCell.size());
     for(localIndex i = oldFineToCoarseSize ; i < subRegion->size();i++)
     {
       newAggregate curAggregate;
       curAggregate.indexFrom = aggregateIndex[i];
     }
     for(localIndex i = oldFineToCoarseSize ; i < m_fineToCoarse.size();i++)
     {
       newAggregate curAggregate;
       curAggregate.indexFrom = aggregateIndex[i];
       auto it = std::find ( ghostAggregates.begin(), ghostAggregates.end(), curAggregate);
       if( it == ghostAggregates.end() )
       {
         curAggregate.center = aggregateCenter[i];
         curAggregate.volume = aggregateVolume[i];
         m_fineToCoarse[i] = offset;
         curAggregate.fineGhostCells.push_back(i);
         curAggregate.newIndex = offset++;
         curAggregate.rankFrom = subRegion->GhostRank()[i];
         ghostAggregates.push_back(curAggregate);
       }
       else
       {
         m_fineToCoarse[i] = it->newIndex;
         it->fineGhostCells.push_back(i);
       }
     }
   });
  int mpiSize;
  int mpiRank;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );

  localIndex nbAggregates = this->size();
  array1d< localIndex > OldNbAggregatesPerRank(mpiSize);
  MPI_Allgather( &nbAggregates, 1, MPI_LONG_LONG,OldNbAggregatesPerRank.data(), 1, MPI_LONG_LONG, MPI_COMM_GEOSX);
   this->resize( this->size() + ghostAggregates.size());
nbAggregates = this->size();
  array1d< localIndex > NewNbAggregatesPerRank(mpiSize);
  MPI_Allgather( &nbAggregates, 1, MPI_LONG_LONG,NewNbAggregatesPerRank.data(), 1, MPI_LONG_LONG, MPI_COMM_GEOSX);
  globalIndex offsetForGlobalIndex = 0;
  for( localIndex i = 0; i < mpiSize; i++)
  {
    offsetForGlobalIndex += OldNbAggregatesPerRank[i];
  }
  for( localIndex i = 0; i < mpiRank; i++)
  {
    offsetForGlobalIndex += NewNbAggregatesPerRank[i] - OldNbAggregatesPerRank[i];
  }
   m_nbFineCellsPerCoarseCell.resize( m_nbFineCellsPerCoarseCell.size() + ghostAggregates.size());
   localIndex counter= 0;
   for(auto curAggregate : ghostAggregates)
   {
     GEOS_LOG_RANK("gloval index "<< offsetForGlobalIndex + counter);
     m_localToGlobalMap[curAggregate.newIndex] = offsetForGlobalIndex + counter;
     m_globalToLocalMap[offsetForGlobalIndex + counter++] = curAggregate.newIndex;
     m_elementVolume[curAggregate.newIndex ] = curAggregate.volume;
     m_elementCenter[curAggregate.newIndex ] = curAggregate.center;
     m_ghostRank[curAggregate.newIndex] = integer_conversion<int>(curAggregate.rankFrom);
     m_nbFineCellsPerCoarseCell[curAggregate.newIndex] = curAggregate.fineGhostCells.size();
     for(localIndex i = 0; i < curAggregate.fineGhostCells.size() ;i++)
     {
       m_fineByAggregates.push_back(curAggregate.fineGhostCells[i]);
     }
   }
   for(auto curAggregate : ghostAggregates)
   {
     m_nbFineCellsPerCoarseCell[curAggregate.newIndex] = m_nbFineCellsPerCoarseCell[curAggregate.newIndex] + m_nbFineCellsPerCoarseCell[curAggregate.newIndex-1];
   }
   GEOS_LOG_RANK("nb new aggrefates" << ghostAggregates.size());
   GEOS_LOG_RANK("new fine to coarse size" << m_fineToCoarse.size());
   GEOS_LOG_RANK("volume size" << m_elementVolume.size());
   GEOS_LOG_RANK("volume first" << m_elementVolume[0]);
   GEOS_LOG_RANK("volume last" << m_elementVolume[this->size()-1]);
   GEOS_LOG_RANK("new local to glocal size " << m_localToGlobalMap.size());
   GEOS_LOG_RANK("new global to local size " << m_globalToLocalMap.size());
}
}
