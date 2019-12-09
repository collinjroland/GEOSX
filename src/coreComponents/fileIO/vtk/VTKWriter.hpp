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
 * @file VTKWriter.hpp
 */

#ifndef GEOSX_FILEIO_VTK_VTKWRITER_HPP_
#define GEOSX_FILEIO_VTK_VTKWRITER_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/RestartFlags.hpp" 
#include "mesh/InterObjectRelation.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

class DomainPartition;

class VTKWriter
{
  public:
  VTKWriter() = delete;

  /*!
   * @brief Initialize the VTK file writer
   * @details This constructor will construct the root file for
   * the output (.pvd)
   * @param[in] name the name of the pvd file
   */
  VTKWriter( string const & name );

  /*!
   * @brief Set the plot level
   * @param[in] plotLevel the plot level. All fields flagged with
   * a plot level inferior or equal to this prescribed value will
   * be output
   */
  void setPlotLevel( const int plotLevel )
  {
    m_plotLevel = dataRepository::IntToPlotLevel(plotLevel);
  }

  /*!
   * @brief Set the binary mode
   * @param[in] binary the binary mode
   */
  void setBinaryMode( const bool binary )
  {
    m_binary = binary;
  }

  /*!
   * @brief Output a file for one time step
   * @param[in] cycle the cycle number
   */
  void write( double const timeStep, 
              DomainPartition const & domain );

  private:
    /// Root file ( .pvd )
    xmlWrapper::xmlDocument m_rootFile;

    /// Unstructured file gathering all vtu files for a time step ( .pvtu )
    pugi::xml_document m_pvtuFile;

    /// Plot level
    dataRepository::PlotLevel m_plotLevel;

    /// Base name of the output
    string m_baseName;

    /// Tells wether or not the output is binary
    bool m_binary;
};
} /* namespace geosx */
#endif /* GEOSX_FILEIO_VTK_VTKWRITER_HPP_ */
