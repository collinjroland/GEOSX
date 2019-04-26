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
 * @file TwoPointFluxApproximation.cpp
 *
 */
#include "TwoPointFluxApproximation.hpp"

#include "meshUtilities/ComputationalGeometry.hpp"
#include "mesh/AggregateElementSubRegion.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include <unordered_set>
#include <ctime>

namespace geosx
{

using namespace dataRepository;

TwoPointFluxApproximation::TwoPointFluxApproximation(std::string const &name,
                                                     ManagedGroup *const parent)
  : FluxApproximationBase(name, parent)
{

}

namespace
{

void makeFullTensor(R1Tensor const & values, R2SymTensor & result)
{
  result = 0.0;
  R1Tensor axis;
  R2SymTensor temp;

  // assemble full tensor from eigen-decomposition
  for (unsigned icoord = 0; icoord < 3; ++icoord)
  {
    // assume principal axis aligned with global coordinate system
    axis = 0.0;
    axis(icoord) = 1.0;

    // XXX: is there a more elegant way to do this?
    temp.dyadic_aa(axis);
    temp *= values(icoord);
    result += temp;
  }
}

}

void TwoPointFluxApproximation::computeCellStencil( DomainPartition const & domain,
                                                    CellStencil & stencil )
{
  MeshBody const * const meshBody = domain.getMeshBody(0);
  MeshLevel const * const mesh = meshBody->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  arrayView2d<localIndex const> const & elemRegionList     = faceManager->elementRegionList();
  arrayView2d<localIndex const> const & elemSubRegionList  = faceManager->elementSubRegionList();
  arrayView2d<localIndex const> const & elemList           = faceManager->elementList();
  arrayView1d<R1Tensor const>   const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const elemCenter =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(
                                        CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const coefficient =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >( m_coeffName );

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  // make a list of region indices to be included
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  constexpr localIndex numElems = 2;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight, faceWeightInv;

  array1d<CellDescriptor> stencilCells(numElems);
  array1d<real64> stencilWeights(numElems);

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve(faceManager->size(), 2);
  real64 const areaTolerance = pow( meshBody->getGlobalLengthScale() * this->m_areaRelTol, 2 );

  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    if (faceGhostRank[kf] >= 0 || elemRegionList[kf][0] == -1 || elemRegionList[kf][1] == -1)
      continue;

    if (!regionFilter.empty() && !(regionFilter.contains(elemRegionList[kf][0]) && regionFilter.contains(elemRegionList[kf][1])))
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );

    if( faceArea < areaTolerance )
      continue;

    faceWeightInv = 0.0;

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      if (elemRegionList[kf][ke] != -1)
      {
        localIndex const er  = elemRegionList[kf][ke];
        localIndex const esr = elemSubRegionList[kf][ke];
        localIndex const ei  = elemList[kf][ke];

        cellToFaceVec = faceCenter;
        cellToFaceVec -= elemCenter[er][esr][ei];

        if (ke == 1)
          cellToFaceVec *= -1.0;

        real64 const c2fDistance = cellToFaceVec.Normalize();

        // assemble full coefficient tensor from principal axis/components
        for(int i = 0; i < 3; i++)
        {
          if(coefficient[er][esr][ei][i] < 1e-30)
          {
            GEOS_LOG_RANK("detect : " << coefficient[er][esr][ei][i]);
            coefficient[er][esr][ei][i] = 1e-15;
          }
        }
        makeFullTensor(coefficient[er][esr][ei], coefTensor);

        faceConormal.AijBj(coefTensor, faceNormal);
        real64 const ht = Dot(cellToFaceVec, faceConormal) * faceArea / c2fDistance;

        faceWeightInv += 1.0 / ht; // XXX: safeguard against div by zero?
      }
    }

    faceWeight = 1.0 / faceWeightInv; // XXX: safeguard against div by zero?

    // ensure consistent normal orientation
    if (Dot(cellToFaceVec, faceNormal) < 0)
      faceWeight *= -1;

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      stencilCells[ke] = { elemRegionList[kf][ke], elemSubRegionList[kf][ke], elemList[kf][ke] };
      stencilWeights[ke] = faceWeight * (ke == 0 ? 1 : -1);
    }
    stencil.add(stencilCells.data(), stencilCells, stencilWeights, kf);
  }
  stencil.compress();
}

void TwoPointFluxApproximation::computeCoarsetencil( DomainPartition * domain,
                                                     CellStencil const & fineStencil,
                                                     CellStencil & coarseStencil,
                                                     std::string const & elementaryPressure1Name,
                                                     std::string const & elementaryPressure2Name,
                                                     std::string const & elementaryPressure3Name )
{
  time_t  start_time = clock();      
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ElementRegion * const elemRegion = elemManager->GetRegion(0); // TODO : still one region / elemsubregion
  AggregateElementSubRegion * const aggregateElement = elemRegion->GetSubRegion("coarse")->group_cast< AggregateElementSubRegion * >();
  auto aggregateGlobalIndex = elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( CellElementSubRegion::viewKeyStruct::aggregateIndexString );
  //const array1d< localIndex > & fineToCoarse = aggregateElement->GetFineToCoarseMap();
  // TODO : will it work for fractures ? (no)

  array1d<CellDescriptor> stencilCells(2);
  array1d<real64> stencilWeights(2);

  std::set< std::pair< localIndex, localIndex > > interfaces;
  fineStencil.forAll( [&] ( StencilCollection<CellDescriptor, real64>::Accessor stencil ) //TODO maybe find a clever way to iterate between coarse interfaces ?
  {
    localIndex const stencilSize = stencil.size();
    if( stencil.size() == 2)
    {
      CellDescriptor const & cell1 = stencil.connectedIndex(0);
      CellDescriptor const & cell2 = stencil.connectedIndex(1);
      localIndex aggregateNumber1 = aggregateElement->m_globalToLocalMap.at(aggregateGlobalIndex[cell1.region][cell1.subRegion][cell1.index]);
      localIndex aggregateNumber2 = aggregateElement->m_globalToLocalMap.at(aggregateGlobalIndex[cell2.region][cell2.subRegion][cell2.index]);
      if( aggregateNumber1 != aggregateNumber2 && interfaces.find(std::make_pair(aggregateNumber1,aggregateNumber2)) == interfaces.end()  && interfaces.find(std::make_pair(aggregateNumber2,aggregateNumber1)) == interfaces.end()) // We find two adjacent aggregates TODO maybe find a clever way to iterate between coarse interfaces ?
      {
        stencilCells[0] = { cell1.region, 1, aggregateNumber1 }; // We can consider here that there is only 1 subregion
        stencilCells[1] = { cell2.region, 1, aggregateNumber2 };

        // Now we compute the transmissibilities
        R1Tensor barycenter1 = aggregateElement->getElementCenter()[aggregateNumber1];
        R1Tensor barycenter2 = aggregateElement->getElementCenter()[aggregateNumber2];
        /*
        std::cout << "======================================"<< std::endl;
        std::cout << "aggregateNumber1 : " << aggregateNumber1 << std::endl;
        std::cout << "aggregateNumber2 : " << aggregateNumber2 << std::endl;
        std::cout << "barycenter1 : " << barycenter1 << std::endl;
        std::cout << "barycenter2 : " << barycenter2 << std::endl;
        */
        barycenter1 -= barycenter2; // normal between the two aggregates
        barycenter1.Normalize();
        //std::cout << "vector between the two aggregates : " << barycenter1 << std::endl;

        int systemSize = integer_conversion< int >(aggregateElement->GetNbCellsPerAggregate( aggregateNumber1 )
                                                 + aggregateElement->GetNbCellsPerAggregate( aggregateNumber2 ));
        Teuchos::LAPACK< int, real64 > lapack;
        Teuchos::SerialDenseMatrix< int, real64 > A(systemSize, 4);
        Teuchos::SerialDenseVector< int, real64 > pTarget(systemSize);

        /// Get the elementary pressures
        ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure1 =
          elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure1Name );
        ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure2 =
          elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure2Name );
        ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure3 =
          elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure3Name );
        int count =0;
        /*
        std::cout << "@@@@@@@@@@@@@@@" << std::endl;
        std::cout << pressure3[0].size() << std::endl;
        for(int i = 0; i < pressure3[0][0].size(); i++)
        {
          std::cout << pressure3[0][0][i] <<std::endl;
        }
        std::cout << "@@@@@@@@@@@@@@@" << std::endl;
        */
        aggregateElement->forGlobalFineCellsInAggregate( aggregateNumber1,
                                                   [&] ( globalIndex fineCellIndexGlobal )
        {
          localIndex fineCellIndex = elemRegion->GetSubRegion(cell1.subRegion)->m_globalToLocalMap.at(fineCellIndexGlobal);
          A(count,0) = pressure1[cell1.region][cell1.subRegion][fineCellIndex];
          A(count,1) = pressure2[cell1.region][cell1.subRegion][fineCellIndex];
          A(count,2) = pressure3[cell1.region][cell1.subRegion][fineCellIndex];
          A(count,3) = 1.;
          /*
 std::cout<< "fine cell index : " << fineCellIndex << " "  <<pressure1[cell1.region][cell1.subRegion][fineCellIndex] << " " 
          << pressure2[cell1.region][cell1.subRegion][fineCellIndex] << " "
          << pressure3[cell1.region][cell1.subRegion][fineCellIndex] << std::endl;
          */
          GEOS_ERROR_IF(fineCellIndex >= elemRegion->GetSubRegion(cell1.subRegion)->size(),"error");
          R1Tensor barycenterFineCell = elemRegion->GetSubRegion(cell1.subRegion)->getElementCenter()[fineCellIndex];
          //std::cout << "fine cell center : " <<  barycenterFineCell << std::endl;
          pTarget(count++) = barycenterFineCell[0]*barycenter1[0]
                             + barycenterFineCell[1]*barycenter1[1]
                             + barycenterFineCell[2]*barycenter1[2];
        });
        aggregateElement->forGlobalFineCellsInAggregate( aggregateNumber2,
                                                   [&] ( globalIndex fineCellIndexGlobal )
        {
          localIndex fineCellIndex = elemRegion->GetSubRegion(cell2.subRegion)->m_globalToLocalMap.at(fineCellIndexGlobal);
          A(count,0) = pressure1[cell2.region][cell2.subRegion][fineCellIndex];
          A(count,1) = pressure2[cell2.region][cell2.subRegion][fineCellIndex];
          A(count,2) = pressure3[cell2.region][cell2.subRegion][fineCellIndex];
          A(count,3) = 1.;
          /*
 std::cout<< "fine cell index : " << fineCellIndex << " "  <<pressure1[cell2.region][cell2.subRegion][fineCellIndex] << " " 
          << pressure2[cell2.region][cell2.subRegion][fineCellIndex] << " "
          << pressure3[cell2.region][cell2.subRegion][fineCellIndex] << std::endl;
          */
          R1Tensor barycenterFineCell = elemRegion->GetSubRegion(cell2.subRegion)->getElementCenter()[fineCellIndex];
          //std::cout << "fine cell center : " <<  barycenterFineCell << std::endl;
          pTarget(count++) = barycenterFineCell[0]*barycenter1[0]
                             + barycenterFineCell[1]*barycenter1[1]
                             + barycenterFineCell[2]*barycenter1[2];
        });

        int info;
        real64  rwork1;
        real64 svd[4];
        int rank;
        /*
        std::cout << "==A==" << std::endl;
        A.print(std::cout);
        std::cout << "==b==" << std::endl;
        pTarget.print(std::cout);
        */

        // Solve the least square system
        lapack.GELSS(systemSize,4,1,A.values(),A.stride(),pTarget.values(),pTarget.stride(),svd,-1,&rank,&rwork1,-1,&info);
        int lwork = static_cast< int > ( rwork1 );
        real64 * rwork = new real64[lwork];
        lapack.GELSS(systemSize,4,1,A.values(),A.stride(),pTarget.values(),pTarget.stride(),svd,-1,&rank,rwork,lwork,&info);
        /*
        std::cout << "==Solution==" << std::endl;
        pTarget.print(std::cout);
        */

        // Computation of coarse-grid flow parameters
        real64 coarseAveragePressure1 = 0.;
        real64 coarseAveragePressure2 = 0.;
        array1d<real64> coarseFlowRate(2);

        aggregateElement->forGlobalFineCellsInAggregate( aggregateNumber1,
                                                   [&] ( globalIndex fineCellIndexGlobal )
        {
          localIndex fineCellIndex = elemRegion->GetSubRegion(cell1.subRegion)->m_globalToLocalMap.at(fineCellIndexGlobal);
          coarseAveragePressure1 += ( pTarget[0] * pressure1[cell1.region][cell1.subRegion][fineCellIndex]
                                   + pTarget[1] * pressure2[cell1.region][cell1.subRegion][fineCellIndex]
                                   + pTarget[2] * pressure3[cell1.region][cell1.subRegion][fineCellIndex]
                                   + pTarget[3] )* elemRegion->GetSubRegion(cell1.subRegion)->getElementVolume()[fineCellIndex];
          /*
          std::cout <<"=========================================" << std::endl;
          std::cout <<  "p target : "<<pTarget[0] << " " << pTarget[1] << " "<< pTarget[2] << " " << pTarget[3]
 << std::endl;
 std::cout<<  "fine pressure : " <<pressure1[cell1.region][cell1.subRegion][fineCellIndex] << " " 
          << pressure2[cell1.region][cell1.subRegion][fineCellIndex] << " "
          << pressure3[cell1.region][cell1.subRegion][fineCellIndex] << std::endl;
          std::cout << "fine volume : "<< elemRegion->GetSubRegion(cell1.subRegion)->getElementVolume()[fineCellIndex] << std::endl;
          std::cout << "to be summed: " << pTarget[0] * pressure1[cell1.region][cell1.subRegion][fineCellIndex] + pTarget[1] * pressure2[cell1.region][cell1.subRegion][fineCellIndex] + pTarget[2] * pressure3[cell1.region][cell1.subRegion][fineCellIndex] + pTarget[3] << std::endl;
          std::cout << "cur Coarse pressure " << coarseAveragePressure1 << std::endl;
          */


        });
        aggregateElement->forGlobalFineCellsInAggregate( aggregateNumber2,
                                                   [&] ( globalIndex fineCellIndexGlobal )
        {
          localIndex fineCellIndex = elemRegion->GetSubRegion(cell2.subRegion)->m_globalToLocalMap.at(fineCellIndexGlobal);
          coarseAveragePressure2 += (pTarget[0] * pressure1[cell2.region][cell2.subRegion][fineCellIndex]
                                   + pTarget[1] * pressure2[cell2.region][cell2.subRegion][fineCellIndex]
                                   + pTarget[2] * pressure3[cell2.region][cell2.subRegion][fineCellIndex]
                                   + pTarget[3] )* elemRegion->GetSubRegion(cell2.subRegion)->getElementVolume()[fineCellIndex];
        });

        /*
        std::cout << "coarse average pressure1 " << coarseAveragePressure1 << std::endl;
        std::cout << "coarse average pressure2 " << coarseAveragePressure2 << std::endl;
        */
        coarseAveragePressure1 /= aggregateElement->getElementVolume()[aggregateNumber1];
        coarseAveragePressure2 /= aggregateElement->getElementVolume()[aggregateNumber2];
        /*
        std::cout << "coarse average pressure1 " << coarseAveragePressure1 << std::endl;
        std::cout << "coarse average pressure2 " << coarseAveragePressure2 << std::endl;
        */

        fineStencil.forAll( [&] ( StencilCollection<CellDescriptor, real64>::Accessor stencilBis )
        {
          localIndex const stencilSizeBis = stencilBis.size();
          if( stencilBis.size() == 2)
          {
            CellDescriptor const & cell1Bis = stencilBis.connectedIndex(0);
            CellDescriptor const & cell2Bis = stencilBis.connectedIndex(1);
            localIndex aggregateNumber1Bis = aggregateElement->m_globalToLocalMap.at(aggregateGlobalIndex[cell1Bis.region][cell1Bis.subRegion][cell1Bis.index]);
            localIndex aggregateNumber2Bis = aggregateElement->m_globalToLocalMap.at(aggregateGlobalIndex[cell2Bis.region][cell2Bis.subRegion][cell2Bis.index]);
            if( aggregateNumber1Bis == aggregateNumber1 &&
                aggregateNumber2Bis == aggregateNumber2 ) // TODO clever way to find the interfaces of adjacent fine cells within aggregates?
            {
              real64 finePressure1 =   pTarget[0] * pressure1[cell1.region][cell1.subRegion][cell1Bis.index]
                                     + pTarget[1] * pressure2[cell1.region][cell1.subRegion][cell1Bis.index]
                                     + pTarget[2] * pressure3[cell1.region][cell1.subRegion][cell1Bis.index]
                                     + pTarget[3];

              real64 finePressure2 =   pTarget[0] * pressure1[cell2.region][cell2.subRegion][cell2Bis.index]
                                     + pTarget[1] * pressure2[cell2.region][cell2.subRegion][cell2Bis.index]
                                     + pTarget[2] * pressure3[cell2.region][cell2.subRegion][cell2Bis.index]
                                     + pTarget[3];
              stencilBis.forAll([&] (CellDescriptor const & cell, real64 w, localIndex const i) -> void
              {
                  coarseFlowRate[i] += w * ( finePressure1 - finePressure2 ); //TODO sign ?
              });
            }
          }
        });
        for( localIndex i = 0; i < 2; i++ )
        {
          stencilWeights[i] = std::fabs(coarseFlowRate[i] / ( coarseAveragePressure1 - coarseAveragePressure2 )) * std::pow(-1,i) ; // TODO sign ?
        //  stencilWeights[i] = 1e-13* std::pow(-1,i) ; // TODO sign ?
          //std::cout << "trans : " <<  stencilWeights[i] << std::endl;
        }
        coarseStencil.add(stencilCells.data(), stencilCells, stencilWeights, 0.);
        interfaces.insert(std::make_pair(aggregateNumber1, aggregateNumber2));
        interfaces.insert(std::make_pair(aggregateNumber2, aggregateNumber1));
      }
    }
});

GEOS_ERROR_IF(true,coarseStencil.numConnections());
  float time1 = static_cast<float>( (clock() - start_time)) / CLOCKS_PER_SEC; 
  GEOS_LOG_RANK("upscaling in "<< time1);
}


void TwoPointFluxApproximation::computeBestCoarsetencil( DomainPartition * domain,
                                                     CellStencil const & fineStencil,
                                                     CellStencil & coarseStencil,
                                                     std::string const & elementaryPressure1Name,
                                                     std::string const & elementaryPressure2Name,
                                                     std::string const & elementaryPressure3Name )
{
  MeshBody * const meshBody = domain->getMeshBody(0);
  MeshLevel * const mesh = meshBody->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ElementRegion * const elemRegion = elemManager->GetRegion(0); // TODO harcoded;
  AggregateElementSubRegion *  aggregateElement = elemManager->GetRegion(0)->GetSubRegion("coarse")->group_cast< AggregateElementSubRegion * >(); // harcoded
  auto aggregateGlobalIndexes = elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( CellElementSubRegion::viewKeyStruct::aggregateIndexString );
  //auto halfTrans = elemManager->ConstructViewAccessor<array2d< real64 > , arrayView2d< real64 > >( ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString );
  //auto halfTrans = elemManager->ConstructViewAccessor<array1d< array1d < real64 > >, arrayView1d< array1d <real64 > > >( ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString );
  //auto halfTrans = elemManager->ConstructViewAccessor<array1d< real64* > , arrayView1d< real64*  > >( ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString );
  //auto nbNeighbors = elemManager->ConstructViewAccessor<array1d< localIndex > , arrayView1d< localIndex  > >("nbNeighbors" );
  auto halfTrans = elemManager->ConstructViewAccessor<array1d< map< globalIndex, real64 >>, arrayView1d< map< globalIndex, real64> > >( ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString );
  //auto neighbors = elemManager->ConstructViewAccessor<array2d< real64 > , arrayView2d< real64 > >( ElementSubRegionBase::viewKeyStruct::neighborsString );

  arrayView2d<localIndex const> const & elemRegionList     = faceManager->elementRegionList();
  arrayView2d<localIndex const> const & elemSubRegionList  = faceManager->elementSubRegionList();
  arrayView2d<localIndex const> const & elemList           = faceManager->elementList();
  arrayView1d<R1Tensor const>   const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const elemCenter =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(
                                        CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const coefficient =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >( "permeability" );

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  // make a list of region indices to be included
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  constexpr localIndex numElems = 2;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight, faceWeightInv;

  array1d<CellDescriptor> stencilCells(numElems);
  array1d<real64> stencilWeights(numElems);

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  real64 const areaTolerance = pow( meshBody->getGlobalLengthScale() * this->m_areaRelTol, 2 );

  std::vector< AggregateCouple > adjacentAggregates;
  adjacentAggregates.reserve(faceManager->size());

  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    // Face on boundary ?
    if ( elemRegionList[kf][0] == -1 || elemRegionList[kf][1] == -1)
      continue;

    if (!regionFilter.empty() && !(regionFilter.contains(elemRegionList[kf][0]) && regionFilter.contains(elemRegionList[kf][1])))
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );

    if( faceArea < areaTolerance )
      continue;
      
    GEOS_ERROR_IF(elemRegionList[kf][0] == -1 || elemRegionList[kf][1] == -1, "Should be never reach");
    localIndex const er0  = elemRegionList[kf][0];
    localIndex const esr0 = elemSubRegionList[kf][0];
    localIndex const ei0  = elemList[kf][0];

    localIndex const er1  = elemRegionList[kf][1];
    localIndex const esr1 = elemSubRegionList[kf][1];
    localIndex const ei1  = elemList[kf][1];

    // These ghost faces at the domain boundary will NEVER be a boundary of two aggragates within this domain (ghost or owned)
    if(ei0 == -1 || ei1 == -1)
      continue;

    Aggregate aggregate0;
    Aggregate aggregate1;

    aggregate0.aggregateGlobalIndex = aggregateGlobalIndexes[er0][esr0][ei0];
    aggregate1.aggregateGlobalIndex = aggregateGlobalIndexes[er1][esr1][ei1];

    // Theses faces are not an aggregate boundary
    if ( aggregate0.aggregateGlobalIndex == aggregate1.aggregateGlobalIndex )
      continue;

    AggregateCouple coupleToFind;
    coupleToFind.aggregate0 = aggregate0;
    coupleToFind.aggregate1 = aggregate1;

    AggregateCouple * couple = nullptr;
    auto aggCoupleIterator = std::find( adjacentAggregates.begin(),
                                        adjacentAggregates.end(),
                                        coupleToFind);
    
    if( aggCoupleIterator == adjacentAggregates.end() )
    {
      // Couple was not added in the vector
      adjacentAggregates.push_back( coupleToFind );
      couple = &adjacentAggregates.back();

      couple->aggregate0.ghostRank = elemManager->GetRegion(er0)
                                              ->GetSubRegion(esr0)
                                              ->GhostRank()[ei0];
      couple->aggregate1.ghostRank = elemManager->GetRegion(er1)
                                              ->GetSubRegion(esr1)
                                              ->GhostRank()[ei1];

      couple->aggregate0.aggregateLocalIndex = aggregateElement->m_globalToLocalMap.at(aggregateGlobalIndexes[er0][esr0][ei0]);
      couple->aggregate1.aggregateLocalIndex = aggregateElement->m_globalToLocalMap.at(aggregateGlobalIndexes[er1][esr1][ei1]);

      couple->aggregate0.er = er0;
      couple->aggregate1.er = er1;

      couple->aggregate0.esr = esr0;
      couple->aggregate1.esr = esr1;
    }
    else
    {
      // Couple already exist in the vector
      couple = &(*aggCoupleIterator);
    }


    // Complete the couple
    couple->faceIndicies.push_back(kf);
    couple->aggregate0.cellBoundInterface.push_back(ei0);
    couple->aggregate1.cellBoundInterface.push_back(ei1);
  }

  // Compute the half transmissibilities
  int mpiSize;
  int mpiRank;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );
  
  /// Get the elementary pressures
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure1 =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure1Name );
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure2 =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure2Name );
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure3 =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure3Name );

  for(auto & aggregateCouple : adjacentAggregates)
  {
    // We first check if the couple of aggregates is not a relation ghost <---> ghost relation
    if( aggregateCouple.GhostGhostRelation() )
      continue;

    Aggregate aggregate0;
    Aggregate aggregate1;

    // aggregate0 will ALWAY be OWNED
    if( aggregateCouple.aggregate0.ghostRank >= 0 )
    {
      aggregate0 = aggregateCouple.aggregate1;
      aggregate1 = aggregateCouple.aggregate0;
    }
    else if( aggregateCouple.aggregate1.ghostRank >= 0 )
    {
      aggregate0 = aggregateCouple.aggregate0;
      aggregate1 = aggregateCouple.aggregate1;
    }
    else
    {
      aggregate0 = aggregateCouple.aggregate0;
      aggregate1 = aggregateCouple.aggregate1;
    }

    computeCoarseHT(domain, 
                   elementaryPressure1Name,
                   elementaryPressure2Name,
                   elementaryPressure3Name,
                   aggregateCouple,
                   aggregate0,
                   aggregate1);
    if( aggregate1.ghostRank < 0 )
    {
      computeCoarseHT(domain, 
                     elementaryPressure1Name,
                     elementaryPressure2Name,
                     elementaryPressure3Name,
                     aggregateCouple,
                     aggregate1,
                     aggregate0);
    }
  }
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator>>( domain->viewKeys.neighbors );

  // We synchronize the half transmissibilities
  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );
  // Now we can build the coarse stencil !!!!!
  for(auto & aggregateCouple : adjacentAggregates)
  {
    // We then check if the couple of aggregates is not a relation ghost <---> owned with rank(ghost) < rank(owned).
    if( aggregateCouple.aggregate0.ghostRank >= 0  && mpiRank > aggregateCouple.aggregate0.ghostRank )
      continue;

    if( aggregateCouple.aggregate1.ghostRank >= 0  && mpiRank > aggregateCouple.aggregate1.ghostRank )
      continue;

    real64 ht0 = halfTrans[0][1][aggregateCouple.aggregate1.aggregateLocalIndex].at(aggregateCouple.aggregate0.aggregateGlobalIndex);
    real64 ht1 = halfTrans[0][1][aggregateCouple.aggregate0.aggregateLocalIndex].at(aggregateCouple.aggregate1.aggregateGlobalIndex);
    real64 transmissibility = (ht0*ht1)/(ht0 + ht1);

    stencilCells[0] = { 0, 1, aggregateCouple.aggregate0.aggregateLocalIndex};
    stencilCells[1] = { 0, 1, aggregateCouple.aggregate1.aggregateLocalIndex};
    stencilWeights[0] = transmissibility;
    stencilWeights[1] = -transmissibility;
    coarseStencil.add(stencilCells.data(), stencilCells, stencilWeights, 0.);
  }
}

void TwoPointFluxApproximation::computeCoarseHT( DomainPartition * domain,
                                                 std::string const & elementaryPressure1Name,
                                                 std::string const & elementaryPressure2Name,
                                                 std::string const & elementaryPressure3Name,
                                                 const AggregateCouple& aggregateCouple,
                                                 const Aggregate& aggregate0,
                                                 const Aggregate& aggregate1)
{
  MeshBody * const meshBody = domain->getMeshBody(0);
  MeshLevel * const mesh = meshBody->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ElementRegion * const elemRegion = elemManager->GetRegion(0); // TODO harcoded;
  AggregateElementSubRegion *  aggregateElement = elemManager->GetRegion(0)->GetSubRegion("coarse")->group_cast< AggregateElementSubRegion * >(); // harcoded
  //auto halfTrans = elemManager->ConstructViewAccessor<array2d< real64 > , arrayView2d< real64 > >( ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString );
  auto halfTrans = elemManager->ConstructViewAccessor<array1d < map< globalIndex, real64 > >, arrayView1d< map< globalIndex, real64> > >( ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString );
//  auto halfTrans = elemManager->ConstructViewAccessor<array1d< array1d < real64 > >, arrayView1d< array1d <real64 > > >( ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString );
  //auto neighbors = elemManager->ConstructViewAccessor<array2d< real64 > , arrayView2d< real64 > >( ElementSubRegionBase::viewKeyStruct::neighborsString );

  arrayView2d<localIndex const> const & elemList           = faceManager->elementList();
  arrayView1d<R1Tensor const>   const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const elemCenter =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(
                                        CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const coefficient =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >( "permeability" );

  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight, faceWeightInv;
  real64 const areaTolerance = pow( meshBody->getGlobalLengthScale() * this->m_areaRelTol, 2 );

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure1 =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure1Name );
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure2 =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure2Name );
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure3 =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure3Name );

  R1Tensor aggregateCenter0 = aggregateElement->getElementCenter()[aggregate0.aggregateLocalIndex];
  R1Tensor aggregateCenter1 = aggregateElement->getElementCenter()[aggregate1.aggregateLocalIndex];

  real64 aggregateVolume0 = aggregateElement->getElementVolume()[aggregate0.aggregateLocalIndex];
  real64 aggregateVolume1 = aggregateElement->getElementVolume()[aggregate1.aggregateLocalIndex];

  R1Tensor agg0toAgg1Direction = aggregateCenter0;
  agg0toAgg1Direction -= aggregateCenter1;
  agg0toAgg1Direction.Normalize();

  // We compute half transmssibility aggregate 0 --- > aggregate 1

  // Least square system on the fine cells of the aggregate0
  int systemSize = integer_conversion< int >( aggregateElement->GetNbCellsPerAggregate( aggregate0.aggregateLocalIndex ) 
      + aggregate1.cellBoundInterface.size());
  GEOS_LOG_RANK("nb cells by agg "<< aggregateElement->GetNbCellsPerAggregate( aggregate0.aggregateLocalIndex ) );
  GEOS_LOG_RANK("nb cells bound "<<aggregate1.cellBoundInterface.size()  );
  Teuchos::LAPACK< int, real64 > lapack;
  Teuchos::SerialDenseMatrix< int, real64 > A(systemSize, 4);
  Teuchos::SerialDenseVector< int, real64 > pTarget(systemSize);

  int count = 0;
  aggregateElement->forGlobalFineCellsInAggregate( aggregate0.aggregateLocalIndex,
      [&] ( globalIndex fineCellIndexGlobal )
      {
      localIndex fineCellIndex = elemRegion->GetSubRegion(0)->m_globalToLocalMap.at(fineCellIndexGlobal); //TODO hardcoded
      A(count,0) = pressure1[0][0][fineCellIndex];
      A(count,1) = pressure2[0][0][fineCellIndex];
      A(count,2) = pressure3[0][0][fineCellIndex];
      A(count,3) = 1.;

      R1Tensor barycenterFineCell = elemRegion->GetSubRegion(aggregate0.esr)->getElementCenter()[fineCellIndex];
      pTarget(count++) = barycenterFineCell[0]*agg0toAgg1Direction[0]
      + barycenterFineCell[1]*agg0toAgg1Direction[1]
      + barycenterFineCell[2]*agg0toAgg1Direction[2];
      });
  for( localIndex cellInOtherSide : aggregate1.cellBoundInterface )
  {
    A(count,0) = pressure1[0][0][cellInOtherSide];
    A(count,1) = pressure2[0][0][cellInOtherSide];
    A(count,2) = pressure3[0][0][cellInOtherSide];
    A(count,3) = 1.;

    R1Tensor barycenterFineCell = elemRegion->GetSubRegion(0)->getElementCenter()[cellInOtherSide];
    pTarget(count++) = barycenterFineCell[0]*agg0toAgg1Direction[0]
      + barycenterFineCell[1]*agg0toAgg1Direction[1]
      + barycenterFineCell[2]*agg0toAgg1Direction[2];
  }
  A.print(std::cout);
  pTarget.print(std::cout);

  // Solve the least square system
  int info;
  real64  rwork1;
  real64 svd[4];
  int rank;
  lapack.GELSS(systemSize,4,1,A.values(),A.stride(),pTarget.values(),pTarget.stride(),svd,-1,&rank,&rwork1,-1,&info);
  int lwork = static_cast< int > ( rwork1 );
  real64 * rwork = new real64[lwork];
  lapack.GELSS(systemSize,4,1,A.values(),A.stride(),pTarget.values(),pTarget.stride(),svd,-1,&rank,rwork,lwork,&info);

  // Compute the coarseAveragePressure in aggregate 0
  real64 coarseAveragePressure = 0.;
  aggregateElement->forGlobalFineCellsInAggregate( aggregate0.aggregateLocalIndex,
      [&] ( globalIndex fineCellIndexGlobal )
      {
      localIndex fineCellIndex = elemRegion->GetSubRegion(0)->m_globalToLocalMap.at(fineCellIndexGlobal); //TODO hardcoded
      coarseAveragePressure += ( pTarget[0] * pressure1[0][0][fineCellIndex]
          + pTarget[1] * pressure2[0][0][fineCellIndex]
          + pTarget[2] * pressure3[0][0][fineCellIndex]
          + pTarget[3] )* elemRegion->GetSubRegion(0)->getElementVolume()[fineCellIndex];
      });
  coarseAveragePressure /= aggregateVolume0;

  // Compute the Pressure at the interface
  real64 pressureAtTheInterface = 0.;
  real64 volumeOfCellsAtTheInterface = 0;

  // Compute the coarse flow rate and the pressure at the interfaces
  real64 coarseFlowRate = 0.;
  for( localIndex kf : aggregateCouple.faceIndicies )
  {
    faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );
    faceWeightInv = 0.0;
    for( localIndex ke = 0; ke < 2 ; ke++)
    {
      localIndex const ei  = elemList[kf][ke];

      real64 pressureEi = pTarget[0] * pressure1[0][0][ei]
        + pTarget[1] * pressure2[0][0][ei]
        + pTarget[2] * pressure3[0][0][ei]
        + pTarget[3];

      real64 volumeEi = elemRegion->GetSubRegion(0)->getElementVolume()[ei];

      pressureAtTheInterface += pressureEi * volumeEi;
      volumeOfCellsAtTheInterface += volumeEi;
      cellToFaceVec = faceCenter;
      cellToFaceVec -= elemCenter[0][0][ei];

      if (ke == 1)
        cellToFaceVec *= -1.0;

      real64 const c2fDistance = cellToFaceVec.Normalize();

      makeFullTensor(coefficient[0][0][ei], coefTensor);

      faceConormal.AijBj(coefTensor, faceNormal);
      real64 const ht = Dot(cellToFaceVec, faceConormal) * faceArea / c2fDistance;
      faceWeightInv += 1.0 / ht; // XXX: safeguard against div by zero?
    }
    faceWeight = 1.0 / faceWeightInv; // XXX: safeguard against div by zero?

    // ensure consistent normal orientation
    if (Dot(cellToFaceVec, faceNormal) < 0)
      faceWeight *= -1;
    for( localIndex ke = 0; ke < 2 ; ke++)
    {
      localIndex const ei  = elemList[kf][ke];

      real64 pressureEi = pTarget[0] * pressure1[0][0][ei]
        + pTarget[1] * pressure2[0][0][ei]
        + pTarget[2] * pressure3[0][0][ei]
        + pTarget[3];
      coarseFlowRate += faceWeight * std::pow(-1,ke) *  pressureEi;
    }
  }
  pressureAtTheInterface /= volumeOfCellsAtTheInterface;
  real64 coarseHalfTransmissibility = std::fabs( coarseFlowRate / ( coarseAveragePressure - pressureAtTheInterface) );
  halfTrans[0][1][aggregate0.aggregateLocalIndex][aggregate1.aggregateGlobalIndex]=coarseHalfTransmissibility;
  //aggregateElement->getReference<array2d<real64>>(ElementSubRegionBase::viewKeyStruct::halfTransmissibilitiesString).rei
  //halfTransMaps[0][1][aggregate0.aggregateGlobalIndex][aggregate1.aggregateGlobalIndex]= coarseHalfTransmissibility;
}


void TwoPointFluxApproximation::computeFractureStencil( DomainPartition const & domain,
                                                        CellStencil & fractureStencil,
                                                        CellStencil & cellStencil )
{

  MeshLevel const * const mesh = domain.getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  EdgeManager const * const edgeManager = mesh->getEdgeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const elemCenter =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const coefficient =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(m_coeffName);

  arrayView1d<real64 const>   const & faceArea   = faceManager->faceArea();
  arrayView1d<R1Tensor const> const & faceCenter = faceManager->faceCenter();
  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();

  arrayView1d<R1Tensor const> const & X = nodeManager->referencePosition();

  // make a list of region indices to be included
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  elemManager->forElementSubRegionsComplete<FaceElementSubRegion>( [&] ( localIndex fractureRegionIndex,
                                                                         localIndex fractureSubRegionIndex,
                                                                         ElementRegion const * const fractureRegion,
                                                                         FaceElementSubRegion const * const fractureSubRegion)
  {
    if (!regionFilter.empty() && !regionFilter.contains(fractureRegionIndex))
      return;

    FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion->faceList();

    array1d<localIndex> const & fractureConnectorIndices =
      fractureRegion->getReference< array1d<localIndex > >( ElementRegion::viewKeyStruct::fractureConnectorIndicesString );

    array1d<array1d<localIndex> > const & fractureConnectors =
      fractureRegion->getReference< array1d<array1d<localIndex> > >( ElementRegion::viewKeyStruct::fractureElementConnectorString );

    FixedToManyElementRelation const & fractureCellConnectors =
      fractureRegion->getReference< FixedToManyElementRelation >( ElementRegion::viewKeyStruct::fractureToCellConnectorString );

    arrayView1d< real64 const > const & aperture = fractureSubRegion->getElementAperture();

    // connections between FaceElements
    for( localIndex fci=0 ; fci<fractureConnectors.size() ; ++fci )
    {
      localIndex const numElems = fractureConnectors[fci].size();
      localIndex const edgeIndex = fractureConnectorIndices[fci];

      array1d<CellDescriptor> stencilCells(numElems);
      array1d<real64> stencilWeights(numElems);

      R1Tensor edgeCenter, edgeLength;
      edgeManager->calculateCenter( edgeIndex, X, edgeCenter );
      edgeManager->calculateLength( edgeIndex, X, edgeLength );

      for( localIndex kfe=0 ; kfe<numElems ; ++kfe )
      {
        localIndex const fractureElementIndex = fractureConnectors[fci][kfe];
        R1Tensor cellCenterToEdgeCenter = edgeCenter;
        cellCenterToEdgeCenter -= faceCenter[ faceMap[fractureElementIndex][0] ];
        stencilCells[kfe] = { fractureRegionIndex, 0, fractureElementIndex };
        // TODO stenciWeights will mean something else once you take out the aperture.
        // We won't be doing the harmonic mean here...etc.
        stencilWeights[kfe] = pow( -1 , kfe ) * pow( aperture[fractureElementIndex], 3) / 12.0 * edgeLength.L2_Norm() / cellCenterToEdgeCenter.L2_Norm();
      }
      fractureStencil.add( stencilCells.data(), stencilCells, stencilWeights, edgeIndex );
    }

    // add connections for FaceElements to/from CellElements.
    {
      array2d< CellDescriptor > cellStencilZeros( fractureCellConnectors.size(0), 2 );

      arrayView2d<localIndex const> const & elemRegionList = fractureCellConnectors.m_toElementRegion;
      arrayView2d<localIndex const> const & elemSubRegionList = fractureCellConnectors.m_toElementSubRegion;
      arrayView2d<localIndex const> const & elemList = fractureCellConnectors.m_toElementIndex;
      for( localIndex kfe=0 ; kfe<fractureCellConnectors.size(0) ; ++kfe )
      {
        localIndex const numElems = fractureCellConnectors.size(1);

        array1d<CellDescriptor> stencilCells(numElems);
        array1d<real64> stencilWeights(numElems);

        R2SymTensor coefTensor;
        R1Tensor cellToFaceVec;
        R1Tensor faceConormal;

        for (localIndex ke = 0; ke < numElems; ++ke)
        {
          localIndex const faceIndex = faceMap[kfe][ke];
          localIndex const er  = elemRegionList[kfe][ke];
          localIndex const esr = elemSubRegionList[kfe][ke];
          localIndex const ei  = elemList[kfe][ke];

          cellToFaceVec = faceCenter[faceIndex];
          cellToFaceVec -= elemCenter[er][esr][ei];

          real64 const c2fDistance = cellToFaceVec.Normalize();

          // assemble full coefficient tensor from principal axis/components
          makeFullTensor(coefficient[er][esr][ei], coefTensor);

          faceConormal.AijBj(coefTensor, faceNormal[faceIndex]);
          real64 const ht = Dot( cellToFaceVec, faceConormal ) * faceArea[faceIndex] / c2fDistance;

          stencilCells[0] = { er, esr, ei};
          stencilWeights[0] = pow(-1,ke) * ht ;

          stencilCells[1] = { fractureRegionIndex, 0, kfe};
          stencilWeights[1] = -pow(-1,ke) * ht ;

          fractureStencil.add( stencilCells.data(), stencilCells, stencilWeights, faceIndex );
        }


        // remove cell connectors from original stencil
        cellStencilZeros[kfe][0] = { elemRegionList[kfe][0], elemSubRegionList[kfe][0], elemList[kfe][0] };
        cellStencilZeros[kfe][1] = { elemRegionList[kfe][1], elemSubRegionList[kfe][1], elemList[kfe][1] };

        cellStencil.zero( faceMap[kfe][0], cellStencilZeros.data() );
      }
    }

  } );

  fractureStencil.compress();

}


void TwoPointFluxApproximation::computeBoundaryStencil( DomainPartition const & domain,
                                                        set<localIndex> const & faceSet,
                                                        BoundaryStencil & stencil )
{
  MeshBody const * const meshBody = domain.getMeshBody(0);
  MeshLevel const * const mesh = meshBody->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  array2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList           = faceManager->elementList();
  r1_array const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const elemCenter =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(
                                        CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const coefficient =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(m_coeffName);

  integer_array const & faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::
                                                                                 viewKeyStruct::
                                                                                 ghostRankString);

  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  // make a list of region indices to be included
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  constexpr localIndex numElems = 2;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight;

  array1d<PointDescriptor> stencilPoints(numElems);
  array1d<real64> stencilWeights(numElems);

  real64 const areaTolerance = pow( meshBody->getGlobalLengthScale() * this->m_areaRelTol, 2 );

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve(faceSet.size(), 2);
  for (localIndex kf : faceSet)
  {
    if (faceGhostRank[kf] >= 0)
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      if (elemRegionList[kf][ke] < 0)
        continue;

      if (!regionFilter.empty() && !regionFilter.contains(elemRegionList[kf][ke]))
        continue;

      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];
      localIndex const ei  = elemList[kf][ke];

      cellToFaceVec = faceCenter;
      cellToFaceVec -= elemCenter[er][esr][ei];

      real64 const c2fDistance = cellToFaceVec.Normalize();

      // assemble full coefficient tensor from principal axis/components
      makeFullTensor(coefficient[er][esr][ei], coefTensor);

      faceConormal.AijBj(coefTensor, faceNormal);
      faceWeight = Dot(cellToFaceVec, faceConormal) * faceArea / c2fDistance;

      // ensure consistent normal orientation
      if (Dot(cellToFaceVec, faceNormal) < 0)
        faceWeight *= -1;

      stencilPoints[0].tag = PointDescriptor::Tag::CELL;
      stencilPoints[0].cellIndex = { er, esr, ei };
      stencilWeights[0] = faceWeight;

      stencilPoints[1].tag = PointDescriptor::Tag::FACE;
      stencilPoints[1].faceIndex = kf;
      stencilWeights[1] = -faceWeight;

      stencil.add(stencilPoints.data(), stencilPoints, stencilWeights, kf );
    }
  }
  stencil.compress();
}


REGISTER_CATALOG_ENTRY(FluxApproximationBase, TwoPointFluxApproximation, std::string const &, ManagedGroup * const)

}
