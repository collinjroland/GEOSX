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
 * @file TwoPhaseWell.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEWELL_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEWELL_HPP_

#include "WellSolverBase.hpp"
#include "physicsSolvers/fluidFlow/TwoPhaseBase.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}
  
class WellElementSubRegion;

/**
 * @class TwoPhaseWell
 *
 * A simple two-phase well solver
 */
class TwoPhaseWell : public WellSolverBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  TwoPhaseWell( const string& name,
                Group * const parent );

  /// deleted default constructor
  TwoPhaseWell() = delete;

  /// deleted copy constructor
  TwoPhaseWell( TwoPhaseWell const & ) = delete;

  /// default move constructor
  TwoPhaseWell( TwoPhaseWell && ) = default;

  /// deleted assignment operator
  TwoPhaseWell & operator=( TwoPhaseWell const & ) = delete;

  /// deleted move operator
  TwoPhaseWell & operator=( TwoPhaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~TwoPhaseWell() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "TwoPhaseWell"; }

  virtual void RegisterDataOnMesh(Group * const meshBodies) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;
  
  virtual bool
  CheckSystemSolution( DomainPartition const * const domain,
                       DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;
  
  virtual void 
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void 
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition * const domain ) override;

  virtual void
  SetupDofs( DomainPartition const * const domain,
             DofManager & dofManager ) const override;

  /**@}*/

  virtual string WellElementDofName() const override { return viewKeyStruct::dofFieldString; }

  virtual string ResElementDofName() const override { return TwoPhaseBase::viewKeyStruct::elemDofFieldString; }

  virtual localIndex NumFluidComponents() const override { return 1; }

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models) on the well
   * @param subRegion the well subRegion containing the well elements and their associated fields
   */
  virtual void UpdateState( WellElementSubRegion * subRegion ) override;
 
  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void AssembleFluxTerms( real64 const time_n,
                          real64 const dt,
                          DomainPartition const * const domain,
                          DofManager const * const dofManager,
                          ParallelMatrix * const matrix,
                          ParallelVector * const rhs ) override;


  /**
   * @Brief assembles the perforation rate terms 
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void AssemblePerforationTerms( real64 const time_n,
                                         real64 const dt,
                                         DomainPartition const * const domain,
                                         DofManager const * const dofManager,
                                         ParallelMatrix * const matrix,
                                         ParallelVector * const rhs ) override;
  
  /**
   * @brief assembles the volume balance terms for all well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void AssembleVolumeBalanceTerms( real64 const time_n,
                                           real64 const dt,
                                           DomainPartition const * const domain,
                                           DofManager const * const dofManager,
                                           ParallelMatrix * const matrix,
                                           ParallelVector * const rhs ) override;

  /**
   * @brief assembles the pressure relations at all connections between well elements except at the well head
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void FormPressureRelations( DomainPartition const * const domain,
                                      DofManager const * const dofManager,
                                      ParallelMatrix * const matrix,
                                      ParallelVector * const rhs ) override;

  /**
   * @brief assembles the control equation for the well head (first connection)
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void FormControlEquation( DomainPartition const * const domain,
                                    DofManager const * const dofManager,
                                    ParallelMatrix * const matrix,
                                    ParallelVector * const rhs ) override;

  
  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    static constexpr auto dofFieldString = "twoPhaseWellVars";

    static constexpr auto resRelPermNameString  = "wellRelPermName";
    static constexpr auto resRelPermIndexString = "elementRelPermIndex";

    // primary solution field
    static constexpr auto pressureString      = TwoPhaseBase::viewKeyStruct::pressureString;
    static constexpr auto deltaPressureString = TwoPhaseBase::viewKeyStruct::deltaPressureString;

    // density (lagged in iteration)
    static constexpr auto mixtureDensityString = "wellElementMixtureDensity";

    // perforation rates
    static constexpr auto perforationRateString          = "perforationRate";
    static constexpr auto dPerforationRate_dPresString   = "dPerforationRate_dPressure";
    static constexpr auto dPerforationRate_dResSatString = "dPerforationRate_dResWettingPhaseSaturation";    
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey resRelPermName  = { resRelPermNameString };
    ViewKey resRelPermIndex = { resRelPermIndexString };
    
    // primary solution field
    ViewKey pressure      = { pressureString };
    ViewKey deltaPressure = { deltaPressureString };
    
    // perforation rates
    ViewKey perforationRate          = { perforationRateString };
    ViewKey dPerforationRate_dPres   = { dPerforationRate_dPresString };
    ViewKey dPerforationRate_dResSat = { dPerforationRate_dResSatString };    
    
  } viewKeysTwoPhaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysTwoPhaseWell;

protected:

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

private:

  /**
   * @brief Setup stored reservoir views into domain data for the current step
   * @param domain the domain containing the well manager to access individual wells
   */
  void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Resize the allocated multidimensional fields
   * @param domain the domain containing the mesh and fields
   *
   * Resize fields along dimensions 1 and 2 (0 is the size of containing object, i.e. element subregion)
   * once the number of phases/components is known (e.g. component fractions)
   */
  void ResizeFields( MeshLevel * const meshLevel );
  
  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void InitializeWells( DomainPartition * const domain ) override;

  /**
   * @brief Check if the controls are viable; if not, switch the controls
   * @param domain the domain containing the well manager to access individual wells
   */
  void CheckWellControlSwitch( DomainPartition * const domain ) override;

  /**
   * @brief Compute all the perforation rates for this well
   * @param well the well with its perforations
   */
  void ComputeAllPerforationRates( WellElementSubRegion const * const subRegion );

  /**
   * @brief Save all the rates and pressures in the well for reporting purposes
   * @param well the well with its perforations
   */
  void RecordWellData( WellElementSubRegion const * const subRegion );

  /// map from the phase indices to the row indices
  array1d<localIndex> m_phaseToRow;
  
  /// index of the wetting phase in the MaterialViewAccessors
  localIndex m_ipw;

  /// index of the non-wetting phase in the MaterialViewAccessors
  localIndex m_ipnw;

  /// index of the phase injected into the reservoir
  localIndex m_injectedPhase;

  /// name of the rel perm constitutive model
  string m_resRelPermName;

  /// index of the rel perm constitutive model in the flow solver
  localIndex m_resRelPermIndex;
  
  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaResPressure;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_resWettingPhaseSat;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaResWettingPhaseSat;

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_dResPhaseMob_dSat;
  
  /// views into reservoir material fields

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_resPhaseDens;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dResPhaseDens_dPres;
};

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEWELL_HPP_
