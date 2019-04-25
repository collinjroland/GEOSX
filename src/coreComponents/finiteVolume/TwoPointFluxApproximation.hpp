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
 * @file TwoPointFluxApproximation.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_

#include "finiteVolume/FluxApproximationBase.hpp"

namespace geosx
{

class TwoPointFluxApproximation : public FluxApproximationBase
{
public:

  static std::string CatalogName() { return "TwoPointFluxApproximation"; }

  TwoPointFluxApproximation() = delete;

  TwoPointFluxApproximation(std::string const & name, dataRepository::ManagedGroup * const parent);

  void computeCoarsetencil( DomainPartition * domain,
                            CellStencil const & fineStencil,
                            CellStencil & coarseStencil,
                            std::string const & elementaryPressure1Name,
                            std::string const & elementaryPressure2Name,
                            std::string const & elementaryPressure3Name );

  void computeBestCoarsetencil( DomainPartition * domain,
                            CellStencil const & fineStencil,
                            CellStencil & coarseStencil,
                            std::string const & elementaryPressure1Name,
                            std::string const & elementaryPressure2Name,
                            std::string const & elementaryPressure3Name );

protected:

  virtual void computeCellStencil( DomainPartition const & domain,
                                   CellStencil & stencil ) override;


  virtual void computeFractureStencil( DomainPartition const & domain,
                                       CellStencil & fractureStencil,
                                       CellStencil & cellStencil ) override;

  virtual void computeBoundaryStencil( DomainPartition const & domain,
                                       set<localIndex> const & faceSet,
                                       BoundaryStencil & stencil ) override;
private:
  struct Aggregate
  {
    localIndex aggregateLocalIndex;
    globalIndex  aggregateGlobalIndex;
    localIndex ghostRank;
    localIndex er;
    localIndex esr;
    array1d< localIndex > cellBoundInterface;
    bool  operator==( const Aggregate & rhs) const
    {
      return aggregateGlobalIndex == rhs.aggregateGlobalIndex;
    }
  };

  struct AggregateCouple
  {
    Aggregate aggregate0;
    Aggregate aggregate1;
    array1d< localIndex > faceIndicies;

    bool operator==( const AggregateCouple & rhs) const
    {
      return (aggregate0 == rhs.aggregate0 && aggregate1 == rhs.aggregate1) || (aggregate0 == rhs.aggregate1 && aggregate1 == rhs.aggregate0);
    }

    /*
    bool GhostOwnedRelation() const
    {
      return aggregate0.ghostRank >= 0 || aggregate1.ghostRank >= 0;
    }
    */

    bool GhostGhostRelation() const
    {
      return aggregate0.ghostRank >= 0 && aggregate1.ghostRank >= 0;
    }
  };
  void computeCoarseHT( DomainPartition * domain,
                        std::string const & elementaryPressure1Name,
                        std::string const & elementaryPressure2Name,
                        std::string const & elementaryPressure3Name,
                        const AggregateCouple& aggregateCouple,
                        const Aggregate& aggregate0,
                        const Aggregate& aggregate1);
                         
};

}


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
