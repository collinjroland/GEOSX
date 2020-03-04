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
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

TwoPhaseWell::TwoPhaseWell( const string & name,
                            Group * const parent )
  :
  WellSolverBase( name, parent )
{
  m_ipw  = -1;
  m_ipnw = -1;
  m_injectedPhase = -1;
  m_numDofPerWellElement = 1;

  this->registerWrapper( viewKeyStruct::resRelPermNameString,  &m_resRelPermName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the relative permeability constitutive model to use");

  this->registerWrapper( viewKeyStruct::resRelPermIndexString, &m_resRelPermIndex, false );
}

void TwoPhaseWell::RegisterDataOnMesh(Group * const meshBodies )
{
  WellSolverBase::RegisterDataOnMesh(meshBodies);

  MeshLevel * const meshLevel = meshBodies->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    // register data on the well element subregion
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
    subRegion->registerWrapper<array1d<real64>>( viewKeyStruct::mixtureDensityString );    

    // register data on the perforations
    PerforationData * const perforationData = subRegion->GetPerforationData();
    perforationData->registerWrapper<array2d<real64>>( viewKeyStruct::perforationRateString );
    perforationData->registerWrapper<array3d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
    perforationData->registerWrapper<array2d<real64>>( viewKeyStruct::dPerforationRate_dResSatString );    
  });
  
}
  
void TwoPhaseWell::InitializePreSubGroups( Group * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  RelativePermeabilityBase const * relPerm = cm->GetConstitutiveRelation<RelativePermeabilityBase>( m_resRelPermName );
  GEOSX_ERROR_IF( relPerm == nullptr, "Relative permeability model " + m_resRelPermName + " not found" );
  m_resRelPermIndex = relPerm->getIndexInParent();
  
  MultiFluidBase const * fluid = cm->GetConstitutiveRelation<MultiFluidBase>( m_fluidName );

  // TODO: see if there is a way to grab this from the flow solver  
  
  // determine the indices of the wetting and non-wetting phases
  // we assume that the oil phase is always present
  if ( (fluid->phaseName( 0 ) == "oil" && fluid->phaseName( 1 ) == "gas") ||
       (fluid->phaseName( 1 ) == "oil" && fluid->phaseName( 0 ) == "water") )
  {
    m_ipw  = 0;
    m_ipnw = 1;
  }
  else if ( (fluid->phaseName( 1 ) == "oil" && fluid->phaseName( 0 ) == "gas") ||
            (fluid->phaseName( 0 ) == "oil" && fluid->phaseName( 1 ) == "water") )
  {
    m_ipw  = 1;
    m_ipnw = 0;
  }
  GEOSX_ERROR_IF( m_ipw == -1 || m_ipnw == -1,
                  "TwoPhaseBase: the accepted phase names are water, oil, and gas");

  // we assume that we always inject water (resp. gas) for oil-water (resp. oil-gas) systems
  m_injectedPhase = ( fluid->phaseName( m_ipw ) == "oil" ) ? m_ipnw : m_ipw ;
  
  // fill the array mapping the phase index to the row offset in the residual 
  m_phaseToRow.resize(TwoPhaseBase::NUM_PHASES); 
  m_phaseToRow[m_ipw]  = TwoPhaseBase::RowOffset::WETTING;
  m_phaseToRow[m_ipnw] = TwoPhaseBase::RowOffset::NONWETTING;

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ResizeFields( meshLevel );
  }
}

void TwoPhaseWell::ResizeFields( MeshLevel * const meshLevel )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    // we have one rate per phase at each perforation 
    PerforationData * const perforationData = subRegion->GetPerforationData();
    perforationData->getReference<array2d<real64>>( viewKeyStruct::perforationRateString ).resizeDimension<1>(numPhases);    
    perforationData->getReference<array3d<real64>>( viewKeyStruct::dPerforationRate_dPresString ).resizeDimension<1,2>(numPhases,2);
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dResSatString ).resizeDimension<1>(numPhases);    
  });
}

  
void TwoPhaseWell::UpdateState( WellElementSubRegion * const subRegion )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;  

  //ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resSat  = m_resWettingPhaseSat;
  //ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dResSat = m_deltaResWettingPhaseSat;   
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseDens = m_resPhaseDens;

  // connectivity info 
  arrayView1d<localIndex const> const & nextWellElemIndex =
    subRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );
  
  // mixture density
  arrayView1d<real64> const & wellElemMixtureDensity =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureDensityString );

  PerforationData const * const perforationData = subRegion->GetPerforationData();

  // get the element region, subregion, index
  arrayView1d<localIndex const> const & resElementRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );
  arrayView1d<localIndex const> const & resElementSubRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );
  arrayView1d<localIndex const> const & resElementIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

  // index of the well element that contains this perforation
  arrayView1d<localIndex const> const & perfWellElemIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

  wellElemMixtureDensity = 0.0;
  
  // loop over the perforations to compute the average densities 
  for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
  {     
    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];

    localIndex const iwelem = perfWellElemIndex[iperf];

    // TODO: what follows is temporary. this has to be completed

    // TODO: use perforation rates here instead of saturation 
    stackArray1d<real64, numPhases> weights( numPhases );
    weights[m_ipw]  = 0.5;//resSat[er][esr][ei] + dResSat[er][esr][ei];
    weights[m_ipnw] = 0.5;//1 - weights[m_ipw];
    real64 sumWeights = 0;
    real64 mixtureDensity  = 0;
    // increment the average mixture density
    for (localIndex ip = 0; ip < numPhases; ++ip)
    {
      sumWeights     += weights[ip];      
      mixtureDensity += weights[ip]
                      * resPhaseDens[er][esr][m_resFluidIndex][ei][0][ip];
    }
    mixtureDensity /= sumWeights; 
    
    // TODO: handle the case with more than one perforation per well element
    wellElemMixtureDensity[iwelem] = mixtureDensity;
    
  }

  for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem) 
  {
    if (wellElemMixtureDensity[iwelem] > 0)
    {
      localIndex currentIndex = iwelem;
      localIndex nextIndex = nextWellElemIndex[currentIndex];
      while (nextIndex >= 0)
      {
        if (wellElemMixtureDensity[nextIndex] > 0)
        {
          break;
        }
        wellElemMixtureDensity[nextIndex] = wellElemMixtureDensity[currentIndex];
        currentIndex = nextIndex;
        nextIndex = nextWellElemIndex[nextIndex];
      }
    }
  }      

  forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
  {
    GEOSX_ERROR_IF( wellElemMixtureDensity[iwelem] <= 0, "Invalid well element mixture density" );
  });      
  
}

void TwoPhaseWell::InitializeWells( DomainPartition * const domain )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;
  
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
    arrayView1d<real64 const> const & wellElemGravCoef =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );

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
    real64 gravCoefControl = 0.0;
    if (subRegion->IsLocallyOwned())
    {

      // get the reference data for this well
      localIndex const iwelemControl = wellControls->GetReferenceWellElementIndex();
      gravCoefControl = wellElemGravCoef[iwelemControl];

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
    MpiWrapper::Broadcast( gravCoefControl, subRegion->GetTopRank() );

    GEOSX_ERROR_IF( pressureControl <= 0, "Invalid well initialization: negative pressure was found" );

    // 3) Estimate the pressures in the well elements using this avgDensity
    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] = pressureControl
        + avgMixtureDensity * ( wellElemGravCoef[iwelem] - gravCoefControl );
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
                       NumDofPerWellElement(),
                       regions );
  
  dofManager.addCoupling( WellElementDofName(),
                          WellElementDofName(),
                          DofManager::Connectivity::Node );
}

void TwoPhaseWell::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                      real64 const GEOSX_UNUSED_PARAM( dt ),
                                      DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                      DofManager const * const GEOSX_UNUSED_PARAM( dofManager ),
                                      ParallelMatrix * const GEOSX_UNUSED_PARAM( matrix ),
                                      ParallelVector * const GEOSX_UNUSED_PARAM( rhs ) )
{
  // nothing to do here
  // in this simple well model, the flow inside the well is not represented  
}


void TwoPhaseWell::AssemblePerforationTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                             real64 const dt,
                                             DomainPartition const * const domain,
                                             DofManager const * const dofManager,
                                             ParallelMatrix * const matrix,
                                             ParallelVector * const rhs )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;
  
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager->getKey( WellElementDofName() );
  string const resDofKey  = dofManager->getKey( ResElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > resDofNumberAccessor =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( resDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex const> >::ViewTypeConst resDofNumber =
    resDofNumberAccessor.toViewConst();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {
    PerforationData const * const perforationData = subRegion->GetPerforationData();

    // compute the local rates for this well
    ComputeAllPerforationRates( subRegion );

    // get the degrees of freedom
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    // get well variables on perforations
    arrayView2d<real64 const> const & perfRate =
      perforationData->getReference<array2d<real64>>( viewKeyStruct::perforationRateString );
    arrayView3d<real64 const> const & dPerfRate_dPres =
      perforationData->getReference<array3d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
    arrayView2d<real64 const> const & dPerfRate_dResSat =
      perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dResSatString );

    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );
    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );
    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // local working variables and arrays
    stackArray1d<globalIndex, numPhases> eqnRowIndices( numPhases );
    stackArray1d<globalIndex, numPhases+1> dofColIndices( numPhases+1 );

    stackArray1d<real64, numPhases> localPerf( numPhases );
    stackArray2d<real64, numPhases*(numPhases+1)> localPerfJacobian(numPhases, numPhases+1);

    localIndex const dp = TwoPhaseBase::ColOffset::DPRES;
    localIndex const dS = TwoPhaseBase::ColOffset::DSAT;
   
    // loop over the perforations and add the rates to the residual and jacobian
    for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
    {
      eqnRowIndices = -1;

      localPerf = 0;
      localPerfJacobian = 0;

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem      = perfWellElemIndex[iperf];
      globalIndex const elemOffset = wellElemDofNumber[iwelem];

      // TODO: add tags here 
      
      // column index on reservoir side
      dofColIndices[0] = resDofNumber[er][esr][ei] + dp;
      dofColIndices[1] = resDofNumber[er][esr][ei] + dS;        
      // column index on well side
      dofColIndices[2] = elemOffset;
     
      // eqn row index for the well
      eqnRowIndices[2] = elemOffset; 
      
      for (localIndex ip = 0; ip < numPhases; ++ip)
      { 
        localIndex const rowId = m_phaseToRow[ip];
        
        // row index on reservoir side
        eqnRowIndices[rowId] = resDofNumber[er][esr][ei]
                             + m_phaseToRow[ip];

        // populate local flux vector and derivatives
        localPerf[rowId]  =  dt * perfRate[iperf][ip];
	  
        localPerfJacobian[rowId][0] = dt * dPerfRate_dPres[iperf][ip][SubRegionTag::RES];
        localPerfJacobian[rowId][1] = dt * dPerfRate_dResSat[iperf][ip];
        localPerfJacobian[rowId][2] = dt * dPerfRate_dPres[iperf][ip][SubRegionTag::WELL];

      }

      rhs->add( eqnRowIndices.data(),
                localPerf.data(),
                numPhases );

      matrix->add( eqnRowIndices.data(),
                   dofColIndices.data(),
                   localPerfJacobian.data(),
                   numPhases, numPhases+1 );
    }

  });
}


void TwoPhaseWell::FormPressureRelations( DomainPartition const * const domain,
                                          DofManager const * const dofManager,
                                          ParallelMatrix * const matrix,
                                          ParallelVector * const rhs )
{

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion const * const subRegion )
  {

    WellControls const * const wellControls = GetWellControls( subRegion );

    // get the degrees of freedom numbers, depth, next well elem index
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );
    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d<real64 const> const & wellElemGravCoef =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      subRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );
    arrayView1d<real64 const> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get secondary data on well elements    
    arrayView1d<real64 const> const & wellElemMixtureDensity =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureDensityString );

    // loop over the well elements to compute the pressure relations between well elements
    forall_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const iwelem )
    {

      if (wellElemGhostRank[iwelem] >= 0)
      {
        return;
      }

      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      if ( iwelemNext >= 0 )  // if iwelemNext < 0, form control equation, not momentum
      {

        // local working variables and arrays
        stackArray1d<globalIndex, 2> dofColIndices( 2 );
        stackArray1d<real64, 2> localPresRelJacobian( 2 );

        dofColIndices        = -1;
        localPresRelJacobian = 0;

        // compute avg density (lagged in iteration)
        real64 const avgMixtureDensity = 0.5 * ( wellElemMixtureDensity[iwelem]
                                               + wellElemMixtureDensity[iwelemNext] );

        // compute depth diff times acceleration
        real64 const gravD = ( wellElemGravCoef[iwelemNext] - wellElemGravCoef[iwelem] );

        // compute the current pressure in the two well elements
        real64 const pressureCurrent = wellElemPressure[iwelem]     + dWellElemPressure[iwelem];
        real64 const pressureNext    = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];

        // compute a coefficient to normalize the momentum equation
        real64 const & targetBHP  = wellControls->GetTargetBHP();
        real64 const normalizer = targetBHP > std::numeric_limits<real64>::epsilon()
                                ? 1.0 / targetBHP
                                : 1.0;

        // compute momentum flux and derivatives
        real64 const localPresRel = - ( pressureNext - pressureCurrent - avgMixtureDensity * gravD ) * normalizer;
        localPresRelJacobian[ElemTag::NEXT]    = - normalizer;
        localPresRelJacobian[ElemTag::CURRENT] =   normalizer;

        // jacobian indices
        globalIndex const offsetNext    = wellElemDofNumber[iwelemNext];
        globalIndex const offsetCurrent = wellElemDofNumber[iwelem];
        globalIndex const eqnRowIndex   = offsetCurrent;
        dofColIndices[ElemTag::NEXT]    = offsetNext;
        dofColIndices[ElemTag::CURRENT] = offsetCurrent;

        rhs->add( &eqnRowIndex,
                  &localPresRel,
                  1 );

        matrix->add( &eqnRowIndex,
                     dofColIndices.data(),
                     localPresRelJacobian.data(),
                     1, 2 );
      }
    });
  });
  
}

void TwoPhaseWell::AssembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                               real64 const GEOSX_UNUSED_PARAM( dt ),
                                               DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                               DofManager const * const GEOSX_UNUSED_PARAM( dofManager ),
                                               ParallelMatrix * const GEOSX_UNUSED_PARAM( matrix ),
                                               ParallelVector * const GEOSX_UNUSED_PARAM( rhs ) )
{
  // nothing to do here
  // in this simple well model, the flow inside the well is not represented
}


void TwoPhaseWell::CheckWellControlSwitch( DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
  // this is not implemented yet
  // TODO: implement this  
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
  // update properties: lagged mixture density
  UpdateStateAll( domain );
  
  dofManager.addVectorToField( solution,
                               WellElementDofName(),
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  CommunicationTools::SynchronizeFields(fieldNames,
                                        domain->getMeshBody(0)->getMeshLevel(0),
                                        domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

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

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_resPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( TwoPhaseBase::viewKeyStruct::pressureString );
  m_deltaResPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( TwoPhaseBase::viewKeyStruct::deltaPressureString );

  m_resWettingPhaseSat =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( TwoPhaseBase::viewKeyStruct::wettingPhaseSatString );
  m_deltaResWettingPhaseSat =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( TwoPhaseBase::viewKeyStruct::deltaWettingPhaseSatString );

  m_resPhaseMob =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( TwoPhaseBase::viewKeyStruct::phaseMobilityString );
  m_dResPhaseMob_dPres =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( TwoPhaseBase::viewKeyStruct::dPhaseMobility_dPressureString );
  m_dResPhaseMob_dSat =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( TwoPhaseBase::viewKeyStruct::dPhaseMobility_dSaturationString );
  
  m_resPhaseDens =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                          constitutiveManager );
  m_dResPhaseDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                          constitutiveManager );
}

void TwoPhaseWell::ComputeAllPerforationRates( WellElementSubRegion const * const subRegion )
{
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & resPressure          = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dResPressure         = m_deltaResPressure;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & resPhaseMob          = m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dResPhaseMob_dPres   = m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dResPhaseMob_dSat    = m_dResPhaseMob_dSat;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseDens        = m_resPhaseDens;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dResPhaseDens_dPres = m_dResPhaseDens_dPres;  

  PerforationData const * const perforationData = subRegion->GetPerforationData();

  // get depth
  arrayView1d<real64 const> const & wellElemGravCoef =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );
    
  // get well primary variables on well elements
  arrayView1d<real64 const> const & wellElemPressure =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dWellElemPressure =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
 
  // get secondary well data on well elements (lagged in iteration)
  arrayView1d<real64 const> const & wellElemMixtureDensity =
    subRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureDensityString );

  // get well variables on perforations
  arrayView1d<real64 const> const & perfGravCoef =
    perforationData->getReference<array1d<real64>>( viewKeyStruct::gravityCoefString );

  arrayView1d<localIndex const> const & perfWellElemIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

  arrayView1d<real64 const> const & perfTransmissibility =
    perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::transmissibilityString );

  arrayView2d<real64> const & perfRate =
    perforationData->getReference<array2d<real64>>( viewKeyStruct::perforationRateString );
  arrayView3d<real64> const & dPerfRate_dPres =
    perforationData->getReference<array3d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
  arrayView2d<real64> const & dPerfRate_dResSat =
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dResSatString );

  // get the element region, subregion, index
  arrayView1d<localIndex const> const & resElementRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );
  arrayView1d<localIndex const> const & resElementSubRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );
  arrayView1d<localIndex const> const & resElementIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

  // local working variables and arrays

  stackArray1d<real64, 2> pressure( 2 );
  stackArray1d<real64, 2> dPressure_dp( 2 );
  stackArray1d<real64, 2> dFlux_dp( 2 );
  stackArray1d<real64, 2> dPotDiff_dp( 2 );  
  stackArray1d<real64, 2> multiplier( 2 );
  
  // loop over the perforations to compute the perforation rates
  for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
  {     

    // clear working arrays
    pressure = 0.0;
    dPressure_dp = 0.0;
        
    // 1) copy the variables from the reservoir and well element
     
    // a) get reservoir variables

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];
    // get the index of the well elem
    localIndex const iwelem = perfWellElemIndex[iperf]; 

    pressure[SubRegionTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dp[SubRegionTag::RES] = 1.0;
    multiplier[SubRegionTag::RES] = 1.0;
    
    // TODO: add a buoyancy term for the reservoir side here 
    
    // b) get well variables

    pressure[SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    dPressure_dp[SubRegionTag::WELL] = 1.0;
    multiplier[SubRegionTag::WELL] = -1.0;

    real64 const gravD = ( perfGravCoef[iperf] - wellElemGravCoef[iwelem] );
    pressure[SubRegionTag::WELL]  += wellElemMixtureDensity[iwelem] * gravD;

    // get transmissibility at the interface
    real64 const trans = perfTransmissibility[iperf]; 
        
    // 2) compute potential difference

    real64 potDiff = 0.0;
    dPotDiff_dp = 0.0;
    for (localIndex i = 0; i < 2; ++i)
    {
      potDiff        += multiplier[i] * trans * pressure[i]; // pressure = pres + dPres
      dPotDiff_dp[i] += multiplier[i] * trans * dPressure_dp[i];
    }
    
    
    // 3) upwinding

    // TODO: double-check what follows. 

    if ( potDiff >= 0 ) // ** reservoir cell is upstream **
    {

      // loop over phases, compute and upwind phase flux
      for (localIndex ip = 0; ip < numPhases; ++ip)
      {
        real64 const phaseCoef     = resPhaseDens[er][esr][m_resFluidIndex][ei][0][ip]
                                   * resPhaseMob[er][esr][ei][ip];
        real64 const dPhaseCoef_dp = dResPhaseDens_dPres[er][esr][m_resFluidIndex][ei][0][ip]
                                   * resPhaseMob[er][esr][ei][ip]
                                   + resPhaseDens[er][esr][m_resFluidIndex][ei][0][ip]
                                   * dResPhaseMob_dPres[er][esr][ei][ip];
        real64 const dPhaseCoef_dS = resPhaseDens[er][esr][m_resFluidIndex][ei][0][ip]
                                   * dResPhaseMob_dSat[er][esr][ei][ip]; 

        // compute the phase flux and derivatives using upstream cell mobility
        perfRate[iperf][ip]                            = phaseCoef * potDiff;
        dPerfRate_dPres[iperf][ip][SubRegionTag::RES]  = dPhaseCoef_dp * potDiff
                                                       + phaseCoef * dPotDiff_dp[SubRegionTag::RES];
        dPerfRate_dPres[iperf][ip][SubRegionTag::WELL] = phaseCoef * dPotDiff_dp[SubRegionTag::WELL];
        dPerfRate_dResSat[iperf][ip]                   = dPhaseCoef_dS * potDiff;
        
      }
    }
    else // ** well is upstream **
    {

      real64 totalMob = 0;
      real64 dTotalMob_dPres = 0;
      real64 dTotalMob_dSat  = 0;
      
      // loop over phases, compute total mobility
      for (localIndex ip = 0; ip < numPhases; ++ip)
      {
        totalMob        += resPhaseMob[er][esr][ei][ip];
        dTotalMob_dPres += dResPhaseMob_dPres[er][esr][ei][ip];
        dTotalMob_dSat  += dResPhaseMob_dSat[er][esr][ei][ip];
      }

      localIndex const ipinj   = m_injectedPhase;
      localIndex const ipother = (m_injectedPhase == m_ipw) ? m_ipnw : m_ipw; 
      
      real64 const phaseCoef     = resPhaseDens[er][esr][m_resFluidIndex][ei][0][ipinj]        * totalMob;
      real64 const dPhaseCoef_dp = dResPhaseDens_dPres[er][esr][m_resFluidIndex][ei][0][ipinj] * totalMob
                                 + resPhaseDens[er][esr][m_resFluidIndex][ei][0][ipinj]        * dTotalMob_dPres;
      real64 const dPhaseCoef_dS = resPhaseDens[er][esr][m_resFluidIndex][ei][0][ipinj]        * dTotalMob_dSat;

      // compute the injected phase flux and derivatives using upstream cell mobility
      perfRate[iperf][ipinj] = phaseCoef * potDiff;
      dPerfRate_dPres[iperf][ipinj][SubRegionTag::RES] = dPhaseCoef_dp * potDiff
        + phaseCoef * dPotDiff_dp[SubRegionTag::RES];
      dPerfRate_dPres[iperf][ipinj][SubRegionTag::WELL] = phaseCoef * dPotDiff_dp[SubRegionTag::WELL];
	dPerfRate_dResSat[iperf][ipinj] = dPhaseCoef_dS * potDiff;

      // set the perforation rate of the other phase to zero
      perfRate[iperf][ipother] = 0;
      dPerfRate_dPres[iperf][ipother][SubRegionTag::RES]  = 0;
      dPerfRate_dPres[iperf][ipother][SubRegionTag::WELL] = 0;
      dPerfRate_dResSat[iperf][ipother] = 0;
      
    }      
  }
}


void TwoPhaseWell::FormControlEquation( DomainPartition const * const domain, 
                                        DofManager const * const dofManager,
                                        ParallelMatrix * const matrix,
                                        ParallelVector * const rhs )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion const * const subRegion )
  {

    if (! subRegion->IsLocallyOwned())
    {
      return;
    }

    WellControls const * const wellControls = GetWellControls( subRegion );

    // get the degrees of freedom 
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );
    
    // get the index of the well element where the control is enforced
    localIndex const iwelemControl = wellControls->GetReferenceWellElementIndex();

    // get well control
    WellControls::Control const control = wellControls->GetControl();

    // BHP control
    if (control == WellControls::Control::BHP)
    {
      // get primary variables on well elements
      arrayView1d<real64 const> const & wellElemPressure =
        subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );
      arrayView1d<real64 const> const & dWellElemPressure =
        subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // get the pressure and compute normalizer
      real64 const currentBHP = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];  
      real64 const & targetBHP = wellControls->GetTargetBHP();
      real64 const normalizer = targetBHP > std::numeric_limits<real64>::epsilon()
                              ? 1.0 / targetBHP
                              : 1.0;

      // control equation is a normalized difference 
      // between current pressure and target pressure
      real64 const controlEqn = ( currentBHP - targetBHP ) * normalizer;
      real64 const dControlEqn_dPres = normalizer;

      globalIndex const eqnRowIndex = wellElemDofNumber[iwelemControl];
      globalIndex const dofColIndex = wellElemDofNumber[iwelemControl];

      rhs->add( &eqnRowIndex,
                &controlEqn,
                1 );

      matrix->add( &eqnRowIndex,
                   &dofColIndex,
                   &dControlEqn_dPres,
                   1, 1 );

    }
    else if (control == WellControls::Control::LIQUIDRATE) // liquid rate control
    {
      GEOSX_ERROR( "This constraint is not supported yet" );
    }
    else
    {
      GEOSX_ERROR_IF( (control != WellControls::Control::BHP) 
                   && (control != WellControls::Control::LIQUIDRATE),
                    "Phase rate contraints for CompositionalMultiphaseWell will be implemented later" );
    }
  });
}


void TwoPhaseWell::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                         real64 const & GEOSX_UNUSED_PARAM( dt ),
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

void TwoPhaseWell::RecordWellData( WellElementSubRegion const * const GEOSX_UNUSED_PARAM( subRegion ) )
{
}

  
REGISTER_CATALOG_ENTRY(SolverBase, TwoPhaseWell, string const &, Group * const)
}// namespace geosx
