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
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/NumericalMethodsManager.hpp"

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
  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( "singlePhaseTPFA" ); // TODO hardcoded.
  fluxApprox->forCellStencils( [&] ( FluxApproximationBase::CellStencil const & stencilCollection )
  {
     std::cout << "UPSCALING" << std::endl;
  });
}

REGISTER_CATALOG_ENTRY( TaskBase, AggregateStencil, std::string const &, ManagedGroup * const )

} /* namespace */
