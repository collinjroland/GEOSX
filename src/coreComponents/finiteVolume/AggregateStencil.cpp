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

#include "AggregateStencil.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/AggregateElementSubRegion.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/TwoPointFluxApproximation.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/FiniteVolume/SinglePhaseFlow.hpp"

namespace geosx
{

using namespace dataRepository;

AggregateStencil::AggregateStencil( std::string const & name,
                        ManagedGroup * const parent ):
  TaskBase( name, parent )
{

  RegisterViewWrapper(viewKeyStruct::fromString, &m_from, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Name of the field to copy.");

  RegisterViewWrapper(viewKeyStruct::toString, &m_to, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Name of the field that will be created and filled.");

  RegisterViewWrapper(viewKeyStruct::targetRegionsString, &m_targetRegions, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Regions on which the copy will be done.");
}

AggregateStencil::~AggregateStencil()
{
}


void AggregateStencil::Execute( real64 const time_n,
                         real64 const dt,
                         integer const cycleNumber,
                         integer const eventCounter,
                         real64 const eventProgress,
                         ManagedGroup * domain )
{
  NumericalMethodsManager * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fineFluxApprox = fvManager->getFluxApproximation( "singlePhaseTPFA" ); // TODO hardcoded.
  TwoPointFluxApproximation * coarseFluxApprox = fvManager->CreateChild("TwoPointFluxApproximation", "coarseSinglePhaseTPFA")->group_cast<TwoPointFluxApproximation*>(); // TODO hardcoded.

  // Compute the coarse stencil
  coarseFluxApprox->computeCoarsetencil(domain->group_cast<DomainPartition*>(), fineFluxApprox->getStencil(),coarseFluxApprox->getStencil(),"pressure1","pressure2","pressure3"); // TODO hardcoded.

  // Upscale porosity and initial Pressure
  MeshLevel * const mesh = domain->group_cast< DomainPartition* >()->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  auto porosity =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( "referencePorosity" ); // TODO hardcoded
  auto elementSubRegionVolumes = elemManager->GetRegion(0)->GetSubRegion(0)->getElementVolume(); // TODO hardcoded
  auto aggregateVolumes = elemManager->GetRegion(0)->GetSubRegion("coarse")->getElementVolume(); // TODO hardcoded
  elemManager->forElementSubRegions<AggregateElementSubRegion>( [&] (  auto * const aggregateRegion) 
  {
    for( localIndex aggregateIndex = 0; aggregateIndex < aggregateRegion->size() ; aggregateIndex++ )
    {
        aggregateRegion->forFineCellsInAggregate( aggregateIndex,
                                                  [&] ( localIndex fineCellIndex )
        {
          porosity[0][1][aggregateIndex] += porosity[0][0][fineCellIndex] * elementSubRegionVolumes[fineCellIndex];
        });
        porosity[0][1][aggregateIndex] /= aggregateVolumes[aggregateIndex];
    }
  });

  auto * solver =
    domain->getParent()->GetGroup<PhysicsSolverManager>( "Solvers" )->GetGroup("SinglePhaseFlow"); // TODO hardcoded

  solver->group_cast< SinglePhaseFlow *>()->switchAggregateMode(true);
  solver->group_cast< SinglePhaseFlow *>()->updateSolid();
  solver->group_cast< SinglePhaseFlow *>()->updateFluid();
  solver->group_cast< SinglePhaseFlow *>()->InitializeAfterAggreg(domain);

  /*
  constitutive::ConstitutiveManager * constitutiveManager = domain->GetGroup<constitutive::ConstitutiveManager>(keys::ConstitutiveManager);
 auto toto =  constitutiveManager->GetGroup("water");
 auto titi =  constitutiveManager->GetGroup("water");
 */
   

  
  
  //const FluxApproximationBase::CellStencil & cellStencil = fluxApprox->getStencil();
  /*
  fluxApprox->forCellStencils( [&] ( FluxApproximationBase::CellStencil const & stencilCollection )
  {
     std::cout << "UPSCALING" << std::endl;
  });
  */
}

REGISTER_CATALOG_ENTRY( TaskBase, AggregateStencil, std::string const &, ManagedGroup * const )

} /* namespace */
