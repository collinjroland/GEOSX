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
  * @file SlurryFluidBase.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace constitutive
{

class SlurryFluidBase : public ConstitutiveBase
{
public:

  SlurryFluidBase( std::string const & name, Group * const parent );

  virtual ~SlurryFluidBase() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override = 0;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  static constexpr localIndex MAX_NUM_COMPONENTS = 3;
  
  // *** SlurryFluidBase-specific interface

  virtual void PointUpdate( real64 const & pressure, real64 const & proppantConcentration, arraySlice1d<real64 const> const & Componentconcentration, real64 const & shearRate, localIndex const k, localIndex const q ) = 0;

  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure, arrayView1d<real64 const> const & proppantConcentration,  arrayView2d<real64 const> const & componentConcentration, arrayView1d<real64 const> const & shearRate) = 0;

  virtual void PointUpdateFluidProperty( real64 const & pressure, arraySlice1d<real64 const> const & Componentconcentration, real64 const & shearRate, localIndex const k, localIndex const q ) = 0;  
  
  localIndex numFluidComponents() const;

  string const & componentName( localIndex ic ) const;

  array1d<real64> const & nIndex() const { return m_nIndices; }
  array1d<real64> const & KIndex() const { return m_Ks; }  

  array2d<real64> const & density() const { return m_density; }
  array2d<real64>       & density()       { return m_density; }

  array2d<real64> const & dDensity_dPressure() const { return m_dDens_dPres; }
  array2d<real64>       & dDensity_dPressure()       { return m_dDens_dPres; }

  array2d<real64> const & dDensity_dProppantConcentration() const { return m_dDens_dProppantConc; }
  array2d<real64>       & dDensity_dProppantConcentration()       { return m_dDens_dProppantConc; }  

  array3d<real64> const & dDensity_dComponentConcentration() const { return m_dDens_dCompConc; }
  array3d<real64>       & dDensity_dComponentConcentration()       { return m_dDens_dCompConc; }  
  
  array2d<real64> const & fluidDensity() const { return m_fluidDensity; }
  array2d<real64>       & fluidDensity()       { return m_fluidDensity; }  

  array2d<real64> const & dFluidDensity_dPressure() const { return m_dFluidDens_dPres; }
  array2d<real64>       & dFluidDensity_dPressure()       { return m_dFluidDens_dPres; }
  
  array3d<real64> const & dFluidDensity_dComponentConcentration() const { return m_dFluidDens_dCompConc; }
  array3d<real64>       & dFluidDensity_dComponentConcentration()       { return m_dFluidDens_dCompConc; }  


  array2d<real64> const & fluidViscosity() const { return m_fluidViscosity; }
  array2d<real64>       & fluidViscosity()       { return m_fluidViscosity; }  

  array2d<real64> const & dFluidViscosity_dPressure() const { return m_dFluidVisc_dPres; }
  array2d<real64>       & dFluidViscosity_dPressure()       { return m_dFluidVisc_dPres; }
  
  array3d<real64> const & dFluidViscosity_dComponentConcentration() const { return m_dFluidVisc_dCompConc; }
  array3d<real64>       & dFluidViscosity_dComponentConcentration()       { return m_dFluidVisc_dCompConc; }  

  
  array2d<real64> const & viscosity() const { return m_viscosity; }
  array2d<real64>       & viscosity()       { return m_viscosity; }

  array2d<real64> const & dViscosity_dPressure() const { return m_dVisc_dPres; }
  array2d<real64>       & dViscosity_dPressure()       { return m_dVisc_dPres; }

  array2d<real64> const & dViscosity_dProppantConcentration() const { return m_dVisc_dProppantConc; }
  array2d<real64>       & dViscosity_dProppantConcentration()       { return m_dVisc_dProppantConc; }  

  array3d<real64> const & dViscosity_dComponentConcentration() const { return m_dVisc_dCompConc; }
  array3d<real64>       & dViscosity_dComponentConcentration()       { return m_dVisc_dCompConc; }

  bool isNewtonianFluid() const { return m_isNewtonianFluid; }  

  // *** Data repository keys

  struct viewKeyStruct
  {

    static constexpr auto componentNamesString       = "componentNames";    
    static constexpr auto defaultDensityString      = "defaultDensity";
    static constexpr auto defaultCompressibilityString      = "defaultCompressibility";
    static constexpr auto defaultViscosityString      = "defaultViscosity";    
    
    static constexpr auto densityString      = "density";
    static constexpr auto dDens_dPresString  = "dDens_dPres";
    static constexpr auto dDens_dProppantConcString  = "dDens_dProppantConc";    
    static constexpr auto dDens_dCompConcString  = "dDens_dCompConc";

    static constexpr auto componentDensityString      = "componentDensity";
    static constexpr auto dCompDens_dPresString  = "dCompDens_dPres";
    static constexpr auto dCompDens_dCompConcString  = "dCompDens_dCompConc";

    static constexpr auto fluidDensityString      = "FluidDensity";
    static constexpr auto dFluidDens_dPresString  = "dFluidDens_dPres";
    static constexpr auto dFluidDens_dCompConcString  = "dFluidDens_dCompConc";        

    static constexpr auto fluidViscosityString      = "FluidViscosity";
    static constexpr auto dFluidVisc_dPresString  = "dFluidVisc_dPres";
    static constexpr auto dFluidVisc_dCompConcString  = "dFluidVisc_dCompConc";        
    
    static constexpr auto viscosityString    = "viscosity";
    static constexpr auto dVisc_dPresString  = "dVisc_dPres";
    static constexpr auto dVisc_dProppantConcString  = "dVisc_dProppantConc";
    static constexpr auto dVisc_dCompConcString  = "dVisc_dCompConc";        
    static constexpr auto flowBehaviorIndexString   = "flowBehaviorIndex";

    static constexpr auto flowConsistencyIndexString   = "flowConsistencyIndex";                

    using ViewKey = dataRepository::ViewKey;

    ViewKey density     = { densityString };
    ViewKey dDens_dPres = { dDens_dPresString };
    ViewKey dDens_dProppantConc = { dDens_dProppantConcString };
    ViewKey dDens_dCompConc = { dDens_dCompConcString };    

    ViewKey fluidDensity = { fluidDensityString };    
    ViewKey dFluidDens_dPres = { dFluidDens_dPresString };
    ViewKey dFluidDens_dCompConc = { dFluidDens_dCompConcString };

    ViewKey fluidViscosity = { fluidViscosityString };    
    ViewKey dFluidVisc_dPres = { dFluidVisc_dPresString };
    ViewKey dFluidVisc_dCompConc = { dFluidVisc_dCompConcString };
    
    ViewKey viscosity   = { viscosityString };
    ViewKey dVisc_dPres = { dVisc_dPresString };
    ViewKey dVisc_dProppantConc = { dVisc_dProppantConcString };
    ViewKey dVisc_dCompConc = { dVisc_dCompConcString };        

    ViewKey flowBehaviorIndex = { flowBehaviorIndexString };
    ViewKey flowConsistencyIndex = { flowConsistencyIndexString };
    
  } viewKeysSlurryFluidBase;

protected:

  virtual void PostProcessInput() override;

  string_array    m_componentNames;

  array1d<real64> m_defaultDensity;
  array1d<real64> m_defaultCompressibility;
  array1d<real64> m_defaultViscosity;  

  array1d<real64> m_nIndices;
  array1d<real64> m_Ks;  

  array2d<real64> m_density;
  array2d<real64> m_dDens_dPres;
  array2d<real64> m_dDens_dProppantConc;  
  array3d<real64> m_dDens_dCompConc;  

  array3d<real64> m_componentDensity;
  array3d<real64> m_dCompDens_dPres;
  array4d<real64> m_dCompDens_dCompConc;  
  
  array2d<real64> m_fluidDensity;  
  array2d<real64> m_dFluidDens_dPres;
  array3d<real64> m_dFluidDens_dCompConc;

  array2d<real64> m_fluidViscosity;  
  array2d<real64> m_dFluidVisc_dPres;
  array3d<real64> m_dFluidVisc_dCompConc;      
  
  array2d<real64> m_viscosity;
  array2d<real64> m_dVisc_dPres;
  array2d<real64> m_dVisc_dProppantConc;  
  array3d<real64> m_dVisc_dCompConc;

  bool m_isNewtonianFluid;

};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP