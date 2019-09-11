/*
 * MPICommData.hpp
 *
 *  Created on: Oct 17, 2018
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_MPI_COMMUNICATIONS_MPICOMMDATA_HPP_
#define SRC_CORECOMPONENTS_MPI_COMMUNICATIONS_MPICOMMDATA_HPP_

#include <set>

#include "common/DataTypes.hpp"



namespace geosx
{

class MPICommData
{
public:

  MPICommData();
  ~MPICommData();


  static std::set<int> & getFreeCommIDs();
  static int reserveCommID();
  static void releaseCommID( int & ID );

  void resize( localIndex numMessages );

  int size;
  int commID;
  int sizeCommID;
  std::map<string, string_array > fieldNames;

  array1d<int> m_sendBufferSize;
  array1d<int> m_receiveBufferSize;

  array1d<buffer_type> m_sendBuffer;
  array1d<buffer_type> m_receiveBuffer;

  array1d<MPI_Request> mpiSendBufferRequest;
  array1d<MPI_Request> mpiRecvBufferRequest;
  array1d<MPI_Status>  mpiSendBufferStatus;
  array1d<MPI_Status>  mpiRecvBufferStatus;

  array1d<MPI_Request> mpiSizeSendBufferRequest;
  array1d<MPI_Request> mpiSizeRecvBufferRequest;
  array1d<MPI_Status>  mpiSizeSendBufferStatus;
  array1d<MPI_Status>  mpiSizeRecvBufferStatus;

};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_MPI_COMMUNICATIONS_MPICOMMDATA_HPP_ */
