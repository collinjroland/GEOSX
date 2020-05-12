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

#include "LASWellGenerator.hpp"
#include "fileIO/las/LASFile.hpp"


namespace geosx
{
using namespace dataRepository;

LASWellGenerator::LASWellGenerator( string const & name, Group * const parent ):
  WellGeneratorBase( name, parent )
{
  registerWrapper( keys::fileName, &m_fileName )->
    setInputFlag( InputFlags::REQUIRED )->
    setSizedFromParent( 0 )->
    setDescription( "Path to the las file" );

  registerWrapper( keys::geometryLogIndexInFile, &m_logIndexToTakeForGeometry )->
    setInputFlag( InputFlags::OPTIONAL )->
    setApplyDefaultValue( -1 )->
    setDescription( "Position of the log to take if there are several log sections defined in the LAS file " );
}

LASWellGenerator::~LASWellGenerator()
{}

void LASWellGenerator::PostProcessInput()
{}

void LASWellGenerator::GeneratePolyLine()
{
  LASFile lasFile;
  lasFile.Load( m_fileName );
  if( lasFile.HasLog( "X" ) && lasFile.HasLog( "Y" ) )
  {
    GeneratePolyLineFromXYZ( lasFile );
  }
  else
  {
    GeneratePolyLineFromDepth( lasFile );
  }

}

void LASWellGenerator::GeneratePolyLineFromXYZ( LASFile const & lasFile )
{
  auto X = lasFile.GetLog( "X" );
  auto Y = lasFile.GetLog( "Y" );
  auto Z = lasFile.GetLog( "TVDSS" );
  localIndex nbLogEntries = lasFile.LogSize( "X" );
  GEOSX_ASSERT( nbLogEntries == lasFile.LogSize( "Y" ) );
  GEOSX_ASSERT( nbLogEntries == lasFile.LogSize( "TVDSS" ) );
  m_polyNodeCoords.resize( nbLogEntries );

  auto lasLineXInfo = lasFile.GetLASLines< LASCurveInformationSection >( "X" );
  auto lasLineYInfo = lasFile.GetLASLines< LASCurveInformationSection >( "Y" );
  auto lasLineZInfo = lasFile.GetLASLines< LASCurveInformationSection >( "TVDSS" );
  GEOSX_ERROR_IF( lasLineXInfo.size() > 1, "X curve info is defined in more than one Curve Information Section" );
  GEOSX_ERROR_IF( lasLineYInfo.size() > 1, "Y curve info is defined in more than one Curve Information Section" );
  GEOSX_ERROR_IF( lasLineZInfo.size() > 1, "TVDSS curve info is defined in more than one Curve Information Section" );
  double factorX = GetFactor( *lasLineXInfo[0] );
  double factorY = GetFactor( *lasLineYInfo[0] );
  double factorZ = GetFactor( *lasLineZInfo[0] );
  m_segmentToPolyNodeMap.resize( nbLogEntries - 1, 2 );
  for( localIndex i = 0; i < nbLogEntries; i++ )
  {
    m_polyNodeCoords[i] = { X[i]*factorX, Y[i]*factorY, Z[i]*factorZ };
    if( i < nbLogEntries - 1 )
    {
      m_segmentToPolyNodeMap[i][0] = i;
      m_segmentToPolyNodeMap[i][1] = i +1;
    }
  }
}

void LASWellGenerator::GeneratePolyLineFromDepth( LASFile const & lasFile )
{

  auto Xs = lasFile.GetLASLines< LASWellInformationSection >( "XCOORD" );
  auto Ys = lasFile.GetLASLines< LASWellInformationSection >( "YCOORD" );
  auto STARTs = lasFile.GetLASLines< LASWellInformationSection >( "STRT" );
  auto STOPs = lasFile.GetLASLines< LASWellInformationSection >( "STOP" );
  GEOSX_ASSERT( Ys.size() == Xs.size() );
  GEOSX_ASSERT( Ys.size() == STARTs.size() );
  GEOSX_ASSERT( Ys.size() == STOPs.size() );
  localIndex wellSectionIndex = 0;
  if( Xs.size() > 1 && m_logIndexToTakeForGeometry == -1 )
  {
    GEOSX_ERROR_IF( Xs.size() > 1 && m_logIndexToTakeForGeometry == -1,
                    "Warning : " << this->getName() << " corresponding LAS file has more than 1 log section defined "
                                 << "please specify the index of the log section you want to take into account to write the well "
                                 << "into the GEOSX data structure. You have to use the keyword " << keys::geometryLogIndexInFile
                                 << ". Taking the first one by default." );
  }
  else if( Xs.size() > 1 && m_logIndexToTakeForGeometry != -1 )
  {
    wellSectionIndex = m_logIndexToTakeForGeometry;
  }

  real64 elev = 0.;
  if( lasFile.HasLine< LASWellInformationSection >( "ELEV" ) )
  {
    auto topDepths = lasFile.GetLASLines< LASWellInformationSection >( "ELEV" );
    elev = topDepths[wellSectionIndex]->GetDataAsReal64();
  }
  real64 const topX = Xs[wellSectionIndex]->GetDataAsReal64();
  real64 const topY = Ys[wellSectionIndex]->GetDataAsReal64();
  real64 const stop = STOPs[wellSectionIndex]->GetDataAsReal64();
  real64 const start = STARTs[wellSectionIndex]->GetDataAsReal64();
  R1Tensor top( topX, topY, start - elev );
  R1Tensor bottom( topX, topY, stop - elev );
  if( start < elev ) // Upward positive z axis
  {
    top[2] *= -1;
    bottom[2] *= -1;
  }
  m_polyNodeCoords.resize( 2 );
  m_polyNodeCoords[0] = top;
  m_polyNodeCoords[1] = bottom;
  m_segmentToPolyNodeMap.resize( 1, 2 );
  m_segmentToPolyNodeMap[0][0] = 0;
  m_segmentToPolyNodeMap[0][1] = 1;
}

real64 LASWellGenerator::GetFactor( LASLine const & lasLine )
{
  real64 factor = 0.;
  if( lasLine.GetUnit() == "F" ||
      lasLine.GetUnit() == "ft" ||
      lasLine.GetUnit() == "feets" ||
      lasLine.GetUnit() == "feet" )
  {
    factor = 0.3048;
  }
  else if( lasLine.GetUnit() == "M" ||
           lasLine.GetUnit() == "meters" ||
           lasLine.GetUnit() == "meter"  ||
           lasLine.GetUnit() == "metre" ||
           lasLine.GetUnit() == "metres" )
  {
    factor = 1.;
  }
  else
  {
    GEOSX_ERROR( "Unit : " << lasLine.GetUnit() << " is not valid, valid units are FEETS and METERS." );
  }
  return factor;
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, LASWellGenerator, std::string const &, Group * const )

} // namespace
