/*
 * MPICommData.cpp
 *
 *  Created on: Oct 17, 2018
 *      Author: settgast
 */

#include "MPICommData.hpp"

namespace geosx
{

MPICommData::MPICommData():
  size(0),
  commID(-1),
  sizeCommID(-1),
  fieldNames(),
  m_sendBufferSize(),
  m_receiveBufferSize(),
  m_sendBuffer(),
  m_receiveBuffer(),
  mpiSendBufferRequest(),
  mpiRecvBufferRequest(),
  mpiSendBufferStatus(),
  mpiRecvBufferStatus(),
  mpiSizeSendBufferRequest(),
  mpiSizeRecvBufferRequest(),
  mpiSizeSendBufferStatus(),
  mpiSizeRecvBufferStatus()
{
  commID = reserveCommID();
  sizeCommID = reserveCommID();
}

MPICommData::~MPICommData()
{
  if( commID >= 0 )
  {
    releaseCommID(commID);
  }

  if( sizeCommID >= 0 )
  {
    releaseCommID(sizeCommID);
  }

}

std::set<int> & MPICommData::getFreeCommIDs()
{
  static std::set<int> commIDs;
  static bool isInitialized = false;

  if( !isInitialized )
  {
    for( int a = 0 ; a < 100 ; ++a )
    {
      commIDs.insert( a );
    }
    isInitialized = true;
  }

  return commIDs;
}


int MPICommData::reserveCommID()
{
  std::set<int> & commIDs = getFreeCommIDs();

  int rval = *( commIDs.begin() );
  commIDs.erase( rval );
  return rval;
}

void MPICommData::releaseCommID( int & ID )
{
  std::set<int> & commIDs = getFreeCommIDs();

  if( commIDs.count( ID ) > 0 )
  {
    GEOS_ERROR( "Attempting to release commID that is already free" );
  }
  commIDs.insert( ID );
  ID = -1;
}

void MPICommData::resize( localIndex numMessages )
{
  m_sendBufferSize.resize( numMessages );
  m_receiveBufferSize.resize( numMessages );
  m_sendBuffer.resize( numMessages );
  m_receiveBuffer.resize( numMessages );
  mpiSendBufferRequest.resize( numMessages );
  mpiRecvBufferRequest.resize( numMessages );
  mpiSendBufferStatus.resize( numMessages );
  mpiRecvBufferStatus.resize( numMessages );
  mpiSizeSendBufferRequest.resize( numMessages );
  mpiSizeRecvBufferRequest.resize( numMessages );
  mpiSizeSendBufferStatus.resize( numMessages );
  mpiSizeRecvBufferStatus.resize( numMessages );
  size = static_cast<int>(numMessages);
}


} /* namespace geosx */
