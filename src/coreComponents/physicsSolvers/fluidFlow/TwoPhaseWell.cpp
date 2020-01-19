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
 * @file TwoPhaseWell.cpp
 */

#include "TwoPhaseWell.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/DomainPartition.hpp"
#include "wells/PerforationData.hpp"
#include "wells/WellElementSubRegion.hpp"
#include "wells/WellControls.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

TwoPhaseWell::TwoPhaseWell( const string & name,
                            Group * const parent )
  :
  WellSolverBase( name, parent )
{
  m_numDofPerWellElement = 1;  
}

void TwoPhaseWell::RegisterDataOnMesh(Group * const meshBodies )
{
  WellSolverBase::RegisterDataOnMesh(meshBodies);

  MeshLevel * const meshLevel = meshBodies->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );

    PerforationData * const perforationData = subRegion->GetPerforationData();
    perforationData->registerWrapper<array1d<real64>>( viewKeyStruct::perforationRateString );
    perforationData->registerWrapper<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
  });
  
}
  
void TwoPhaseWell::InitializePreSubGroups( Group * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {

    PerforationData * const perforationData = subRegion->GetPerforationData();
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString ).resizeDimension<1>(2);

  });
}
  
void TwoPhaseWell::UpdateState( WellElementSubRegion * const GEOSX_UNUSED_ARG( subRegion ) )
{
}

void TwoPhaseWell::InitializeWells( DomainPartition * const domain )
{
  localIndex constexpr numPhases = 2;//TwoPhaseBase::NUM_PHASES;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resPressure  = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resSat       = m_resWettingPhaseSat;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseDens = m_resPhaseDens;

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion * const subRegion )
  {
    WellControls const * const wellControls = GetWellControls( subRegion );
    PerforationData const * const perforationData = subRegion->GetPerforationData();

    // get the info stored on well elements
    arrayView1d<real64 const> const & wellElemGravDepth =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );

    // get well primary variables on well elements
    arrayView1d<real64> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );
    
    // define a reservoir pressure used for initialization
    real64 resPres = ( wellControls->GetType() == WellControls::Type::PRODUCER )
                   ? 1e20 : 0;

    // 1) Loop over all perforations to compute an average density
    real64 avgMixtureDensity = 0;
    for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
    {
      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      stackArray1d<real64, numPhases> phaseSat(numPhases);
      phaseSat[m_ipw]  = resSat[er][esr][ei];
      phaseSat[m_ipnw] = 1.0 - resSat[er][esr][ei];
      
      // increment the average mixture density
      for (localIndex ip = 0; ip < numPhases; ++ip)
      {
        real64 const resDensity = resPhaseDens[er][esr][m_resFluidIndex][ei][0][ip];
        avgMixtureDensity += phaseSat[ip] * resDensity;
      }

      // save min pressure for producer
      if ( wellControls->GetType() == WellControls::Type::PRODUCER &&
           resPres > resPressure[er][esr][ei] )
      { 
        resPres = resPressure[er][esr][ei];
      }
      // save max pressure for injector
      else if ( wellControls->GetType() == WellControls::Type::INJECTOR &&
                resPres < resPressure[er][esr][ei] )
      { 
        resPres = resPressure[er][esr][ei];
      }
    }

    // communicate the pressures to the ranks without perforations
    // this will be used to initialize the pressure, starting by the owner rank
    if ( wellControls->GetType() == WellControls::Type::PRODUCER )
    { 
      resPres = MpiWrapper::Min( resPres );
    }
    else if ( wellControls->GetType() == WellControls::Type::INJECTOR )
    { 
      resPres = MpiWrapper::Max( resPres );
    }
    
    avgMixtureDensity = MpiWrapper::Sum( avgMixtureDensity );

    globalIndex const numPerforationsGlobal = perforationData->GetNumPerforationsGlobal();
    avgMixtureDensity /= numPerforationsGlobal;

    real64 pressureControl = 0.0;
    real64 gravDepthControl = 0.0;
    if (subRegion->IsLocallyOwned())
    {

      // get the reference data for this well
      localIndex const iwelemControl = wellControls->GetReferenceWellElementIndex();
      gravDepthControl = wellElemGravDepth[iwelemControl];

      // 2) Initialize the reference pressure
      real64 const & targetBHP = wellControls->GetTargetBHP();
      if (wellControls->GetControl() == WellControls::Control::BHP)
      {
        // if pressure constraint, set the ref pressure at the constraint
        pressureControl = targetBHP;
      }
      else // rate control
      {
        // if rate constraint, set the ref pressure slightly 
        // above/below the target pressure depending on well type
        pressureControl = 
          (wellControls->GetType() == WellControls::Type::PRODUCER)
          ? 0.5 * resPres 
          : 2.0 * resPres;
      }

      wellElemPressure[iwelemControl] = pressureControl;
    }

    // TODO optimize
    MpiWrapper::Broadcast( pressureControl, subRegion->GetTopRank() );
    MpiWrapper::Broadcast( gravDepthControl, subRegion->GetTopRank() );

    GEOSX_ERROR_IF( pressureControl <= 0, "Invalid well initialization: negative pressure was found" );

    // 3) Estimate the pressures in the well elements using this avgDensity
    integer const gravityFlag = m_gravityFlag;

    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] = pressureControl
        + ( gravityFlag 
          ? avgMixtureDensity * ( wellElemGravDepth[iwelem] - gravDepthControl ) 
          : 0 );
    });

    // 4) Recompute the pressure-dependent properties
    UpdateState( subRegion );

  });
}


void TwoPhaseWell::SetupDofs( DomainPartition const * const domain, 
                              DofManager & dofManager ) const
{
  MeshLevel const * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  array1d<string> regions;
  elemManager->forElementRegions<WellElementRegion>( [&]( WellElementRegion const * const region )
  {
    regions.push_back( region->getName() );
  } );

  dofManager.addField( WellElementDofName(),
                       DofManager::Location::Elem,
                       DofManager::Connectivity::Node,
                       NumDofPerWellElement(),
                       regions );

}

void TwoPhaseWell::AssembleFluxTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
                                      real64 const GEOSX_UNUSED_ARG( dt ),
                                      DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                      DofManager const * const GEOSX_UNUSED_ARG( dofManager ),
                                      ParallelMatrix * const GEOSX_UNUSED_ARG( matrix ),
                                      ParallelVector * const GEOSX_UNUSED_ARG( rhs ) )
{
  // nothing to do here 
}


void TwoPhaseWell::AssemblePerforationTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
                                             real64 const GEOSX_UNUSED_ARG( dt ),
                                             DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                             DofManager const * const GEOSX_UNUSED_ARG( dofManager ),
                                             ParallelMatrix * const GEOSX_UNUSED_ARG( matrix ),
                                             ParallelVector * const GEOSX_UNUSED_ARG( rhs ) )
{
}


void TwoPhaseWell::FormPressureRelations( DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                          DofManager const * const GEOSX_UNUSED_ARG( dofManager ),
                                          ParallelMatrix * const GEOSX_UNUSED_ARG( matrix ),
                                          ParallelVector * const GEOSX_UNUSED_ARG( rhs ) )
{
}

void TwoPhaseWell::AssembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
                                               real64 const GEOSX_UNUSED_ARG( dt ),
                                               DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                               DofManager const * const GEOSX_UNUSED_ARG( dofManager ),
                                               ParallelMatrix * const GEOSX_UNUSED_ARG( matrix ),
                                               ParallelVector * const GEOSX_UNUSED_ARG( rhs ) )
{
  // nothing to do here
}


void TwoPhaseWell::CheckWellControlSwitch( DomainPartition * const GEOSX_UNUSED_ARG( domain ) )
{
}


real64
TwoPhaseWell::CalculateResidualNorm( DomainPartition const * const domain,
                                     DofManager const & dofManager,
                                     ParallelVector const & rhs )
{
  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager.getKey( WellElementDofName() );

  real64 residualNorm = 0;
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {

    // get the degree of freedom numbers
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        localIndex const lid    = rhs.getLocalRowID( wellElemDofNumber[iwelem] );
        real64 const normalizer = 1; // TODO: compute normalizer
        real64 const val = localResidual[lid] / normalizer;
        residualNorm += val * val;
      }
    }
  });

  // compute global residual norm
  real64 globalResidualNorm;
  MpiWrapper::allReduce(&residualNorm, &globalResidualNorm, 1, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

bool
TwoPhaseWell::CheckSystemSolution( DomainPartition const * const domain,
                                   DofManager const & dofManager,
                                   ParallelVector const & solution,
                                   real64 const scalingFactor )
{
  // get the update
  real64 const * localSolution = solution.extractLocalVector();

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  int isInvalidLocal = 0;

  string const wellDofKey = dofManager.getKey( WellElementDofName() );

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {

    // get the degree of freedom numbers on well elements
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get a reference to the primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {

      if (wellElemGhostRank[iwelem] < 0)
      {
        // extract pressure solution and apply to dP
        localIndex const lid =  solution.getLocalRowID( wellElemDofNumber[iwelem] );
        real64 const newPres = wellElemPressure[iwelem] + dWellElemPressure[iwelem]
                             + scalingFactor * localSolution[lid];

        if (newPres < 0.0)
        {
          isInvalidLocal = 1;
        }
      }
    }
  });

  int isInvalidGlobal;
  MpiWrapper::allReduce(&isInvalidLocal, &isInvalidGlobal, 1, MPI_SUM, MPI_COMM_GEOSX);
 
  bool isValid = (isInvalidGlobal == 0);
  return isValid;
}

void
TwoPhaseWell::ApplySystemSolution( DofManager const & dofManager,
                                   ParallelVector const & solution, 
                                   real64 const scalingFactor,
                                   DomainPartition * const domain )
{
  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    dofManager.addVectorToField( solution,
                                 WellElementDofName(),
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaPressureString,
                                 0, 1 );
  });

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  CommunicationTools::SynchronizeFields(fieldNames,
                                        domain->getMeshBody(0)->getMeshLevel(0),
                                        domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  // update properties
  UpdateStateAll( domain );
}

void TwoPhaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {

    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      dWellElemPressure[iwelem] = 0;
    });
  });

  // call constitutive models
  UpdateStateAll( domain );
}


void TwoPhaseWell::ResetViews(DomainPartition * const domain )
{
  WellSolverBase::ResetViews(domain);

  /*
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_resPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( TwoPhaseFlow::viewKeyStruct::pressureString );

  m_deltaResPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( TwoPhaseFlow::viewKeyStruct::deltaPressureString );

  m_resSat =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( TwoPhaseFlow::viewKeyStruct::wettingPhaseSatString );

  m_deltaResSat =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( TwoPhaseFlow::viewKeyStruct::deltaWettingPhaseSatString );
  
  m_resPhaseDens =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                          constitutiveManager );
  m_dResPhaseDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDens_dPresString,
                                                                                          constitutiveManager );
  */
}

void TwoPhaseWell::ComputeAllPerforationRates( WellElementSubRegion const * const GEOSX_UNUSED_ARG( subRegion ) )
{
}


void TwoPhaseWell::FormControlEquation( DomainPartition const * const GEOSX_UNUSED_ARG( domain ), 
 			                DofManager const * const GEOSX_UNUSED_ARG( dofManager ),
					ParallelMatrix * const GEOSX_UNUSED_ARG( matrix ),
                                        ParallelVector * const GEOSX_UNUSED_ARG( rhs ) )
{
}


void TwoPhaseWell::ImplicitStepComplete( real64 const & GEOSX_UNUSED_ARG( time ),
                                         real64 const & GEOSX_UNUSED_ARG( dt ),
                                         DomainPartition * const domain )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {
    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & wellElemPressure  =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
    }

    // TODO: improve well data output
    /*
    int mpiSize = CommunicationTools::Comm_size(MPI_COMM_GEOSX) ;
    if (mpiSize == 1)
    {
      RecordWellData( subRegion );
    }
    */
  });
}

void TwoPhaseWell::RecordWellData( WellElementSubRegion const * const GEOSX_UNUSED_ARG( subRegion ) )
{
}

  
REGISTER_CATALOG_ENTRY(SolverBase, TwoPhaseWell, string const &, Group * const)
}// namespace geosx
