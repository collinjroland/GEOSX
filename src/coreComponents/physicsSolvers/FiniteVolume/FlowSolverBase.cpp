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

/**
 * @file FlowSolverBase.cpp
 */

#include "FlowSolverBase.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/AggregateElementSubRegion.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

FlowSolverBase::FlowSolverBase( std::string const & name,
                                ManagedGroup * const parent )
  : SolverBase( name, parent ),
    m_gravityFlag(1),
    m_fluidName(),
    m_solidName(),
    m_fluidIndex(),
    m_solidIndex(),
    m_poroElasticFlag(0),
    m_numDofPerCell(0),
    m_elemGhostRank(),
    m_volume(),
    m_gravDepth(),
    m_porosityRef()
{
  RegisterViewWrapper( viewKeyStruct::gravityFlagString, &m_gravityFlag, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Flag that enables/disables gravity");

  this->RegisterViewWrapper( viewKeyStruct::discretizationString, &m_discretizationName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of discretization object to use for this solver.");

  this->RegisterViewWrapper( viewKeyStruct::fluidNameString,  &m_fluidName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of fluid constitutive object to use for this solver.");

  this->RegisterViewWrapper( viewKeyStruct::solidNameString,  &m_solidName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of solid constitutive object to use for this solver");

  this->RegisterViewWrapper( viewKeyStruct::fluidIndexString, &m_fluidIndex, false );
  this->RegisterViewWrapper( viewKeyStruct::solidIndexString, &m_solidIndex, false );


}

void FlowSolverBase::RegisterDataOnMesh( ManagedGroup * const MeshBodies )
{
  SolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & subgroup : MeshBodies->GetSubGroups() )
  {
    MeshBody * const meshBody = subgroup.second->group_cast<MeshBody *>();
    MeshLevel * const mesh = meshBody->getMeshLevel(0);

    applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
    {
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::referencePorosityString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->RegisterViewWrapper< array1d<R1Tensor> >( viewKeyStruct::permeabilityString )->setPlotLevel(PlotLevel::LEVEL_0);
      auto toto = subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::gravityDepthString )->setApplyDefaultValue( 0.0 );
     std::cout << "size of the region "<< toto->size() << std::endl;
    });

    mesh->getElemManager()->forElementSubRegions<AggregateElementSubRegion>( [&] (  auto * const aggregateRegion) 
    {
      std::cout << "Aggretion sub region name " << aggregateRegion->getName() << std::endl;
      aggregateRegion->template RegisterViewWrapper< array1d<real64> >( viewKeyStruct::referencePorosityString )->setPlotLevel(PlotLevel::LEVEL_0);
      aggregateRegion->template RegisterViewWrapper< array1d<R1Tensor> >( viewKeyStruct::permeabilityString )->setPlotLevel(PlotLevel::LEVEL_0);
     auto toto = aggregateRegion->template RegisterViewWrapper< array1d<real64> >( viewKeyStruct::gravityDepthString )->setApplyDefaultValue( 0.0 );
     std::cout << "size of the region "<< toto->size() << " " << aggregateRegion->size() << std::endl;
    });

    FaceManager * const faceManager = mesh->getFaceManager();
    faceManager->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::gravityDepthString )->setApplyDefaultValue( 0.0 );
  }
}

void FlowSolverBase::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * fluid  = cm->GetConstitituveRelation<ConstitutiveBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Fluid model " + m_fluidName + " not found" );
  m_fluidIndex = fluid->getIndexInParent();

  ConstitutiveBase const * solid  = cm->GetConstitituveRelation<ConstitutiveBase>( m_solidName );
  GEOS_ERROR_IF( solid == nullptr, "Solid model " + m_solidName + " not found" );
  m_solidIndex = solid->getIndexInParent();
}

void FlowSolverBase::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePostInitialConditions_PreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  ResetViews( domain );

  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData(domain);
}

void FlowSolverBase::PrecomputeData( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();

  R1Tensor const & gravityVector = getGravityVector();

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    arrayView1d<R1Tensor const> const & elemCenter =
      subRegion->getReference<array1d<R1Tensor>>( CellBlock::viewKeyStruct::elementCenterString );

    arrayView1d<real64> const & gravityDepth =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex a )
    {
      gravityDepth[a] = Dot( elemCenter[a], gravityVector );
    } );
  } );

  {
    arrayView1d<R1Tensor const> const & faceCenter =
      faceManager->getReference<array1d<R1Tensor>>(FaceManager::viewKeyStruct::faceCenterString);

    arrayView1d<real64> const & gravityDepth =
      faceManager->getReference<array1d<real64>>(viewKeyStruct::gravityDepthString);

    forall_in_range<elemPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex a )
    {
      gravityDepth[a] = Dot( faceCenter[a], gravityVector );
    } );
  }
}

FlowSolverBase::~FlowSolverBase() = default;

void FlowSolverBase::ResetViews( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_elemGhostRank =
    elemManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  m_volume =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( CellElementSubRegion::viewKeyStruct::elementVolumeString );
  m_gravDepth =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::gravityDepthString );
  m_porosityRef =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::referencePorosityString );
}


} // namespace geosx
