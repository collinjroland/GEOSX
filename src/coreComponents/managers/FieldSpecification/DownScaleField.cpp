/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * DownScaleright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include "DownScaleField.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/AggregateElementSubRegion.hpp"
#include "FieldSpecificationBase.hpp"

namespace geosx
{

using namespace dataRepository;

DownScaleField::DownScaleField( std::string const & name,
                        ManagedGroup * const parent ):
  TaskBase( name, parent )
{

  RegisterViewWrapper(viewKeyStruct::fieldNameString, &m_fieldName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Name of the field to copy.");
}

DownScaleField::~DownScaleField()
{
}


void DownScaleField::Execute( real64 const time_n,
                         real64 const dt,
                         integer const cycleNumber,
                         integer const eventCounter,
                         real64 const eventProgress,
                         ManagedGroup * domain )
{
  DomainPartition * domainCast = domain->group_cast<DomainPartition*>(domain);
  MeshBody * meshBody = domainCast->getMeshBody(0);
  MeshLevel * meshLevel = meshBody->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  auto pressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( m_fieldName ); //TODO harcoded
  AggregateElementSubRegion * aggregateRegion = elemManager->GetRegion(0)->GetSubRegion("coarse")->group_cast< AggregateElementSubRegion *>();
  CellElementSubRegion * fineRegion = elemManager->GetRegion(0)->GetSubRegion(0)->group_cast< CellElementSubRegion *>();
      auto & aggregateIndexSave =
        fineRegion->template getWrapper< array1d< globalIndex > > (CellElementSubRegion::viewKeyStruct::aggregateIndexString)->reference();
  for( localIndex fineCellIndex  = 0; fineCellIndex < fineRegion->size(); fineCellIndex++)
  {
    localIndex aggregateIndex = aggregateRegion->m_globalToLocalMap.at(aggregateIndexSave[fineCellIndex]);
    pressure[0][0][fineCellIndex] = pressure[0][1][aggregateIndex];
  }
  /*
  for( localIndex aggregateIndex = 0; aggregateIndex < aggregateRegion->size(); aggregateIndex++ )
  {
    aggregateRegion->forFineCellsInAggregate( aggregateIndex,
                                              [&] ( localIndex fineCellIndex )
    {
      pressure[0][0][fineCellIndex] = pressure[0][1][aggregateIndex];
      //pressure[0][0][fineCellIndex] = aggregateIndex;
    });

  }
  */
}

REGISTER_CATALOG_ENTRY( TaskBase, DownScaleField, std::string const &, ManagedGroup * const )

} /* namespace */
