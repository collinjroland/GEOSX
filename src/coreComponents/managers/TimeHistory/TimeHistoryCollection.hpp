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
 * @file TimeHistoryCollection.hpp
 */

#ifndef GEOSX_TimeHistoryCollection_HPP_
#define GEOSX_TimeHistoryCollection_HPP_

#include "managers/Tasks/TaskBase.hpp"
#include "managers/TimeHistory/HistoryDataSpec.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/ProblemManager.hpp"
#include "dataRepository/BufferOpsDevice.hpp"

#include <functional>

namespace geosx
{
using namespace dataRepository;

/**
 * @class HistoryCollection
 *
 * A task class for serializing time history data into a buffer for later I/O.
 */
class HistoryCollection : public TaskBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group(std::string const & name, Group * const parent)
  HistoryCollection( string const & name, Group * parent ):
    TaskBase( name, parent ),
    m_buffer_call()
  {   }

  /**
   * @brief Get the metadata for what this collector collects.
   * @param problem_group The problem manager cast to a group.
   * @return A HistoryMetadata object describing  the history data being collected by this collector.
   */
  virtual HistoryMetadata GetMetadata( Group * problem_group )
  {
    GEOSX_UNUSED_VAR( problem_group );
    return HistoryMetadata( "null", 0, std::type_index( typeid(nullptr)) );
  }

  /**
   * @brief Get the number of collectors of meta-information (set indices, etc) writing time-independent information during initialization.
   * @return The number of collectors of meta-information for this collector.
   */
  virtual localIndex GetNumMetaCollectors( ) const { return 0; }

  /**
   * @brief Get a pointer to a collector of meta-information for this collector.
   * @param problem_group The ProblemManager cast to a group.
   * @param meta_idx Which of the meta-info collectors to return. (see HistoryCollection::GetNumMetaCollectors( ) ).
   * @param meta_rank_offset The offset for this rank for the meta-info collector, used to number index metadata consistently across the
   * simulation.
   * @return A unique pointer to the HistoryCollection object used for meta-info collection. Intented to fall out of scope and desctruct
   * immediately
   *         after being used to perform output during simulation initialization.
   */
  virtual std::unique_ptr< HistoryCollection > GetMetaCollector( Group * problem_group, localIndex meta_idx, globalIndex meta_rank_offset )
  {
    GEOSX_UNUSED_VAR( problem_group );
    GEOSX_UNUSED_VAR( meta_idx );
    GEOSX_UNUSED_VAR( meta_rank_offset );
    return std::unique_ptr< HistoryCollection >( nullptr );
  }

  /**
   * @brief Collect history information into the provided buffer. Typically called from HistoryCollection::Execute .
   * @param domain The ProblemDomain cast to a group.
   * @param time_n The current simulation time.
   * @param dt The current simulation time delta.
   * @param buffer A properly-sized buffer to serialize history data into.
   */
  virtual void Collect( Group * domain, real64 const time_n, real64 const dt, buffer_unit_type * & buffer )
  {
    GEOSX_UNUSED_VAR( domain );
    GEOSX_UNUSED_VAR( time_n );
    GEOSX_UNUSED_VAR( dt );
    GEOSX_UNUSED_VAR( buffer );
  }

  /**
   * @brief Collects history data.
   * @copydoc EventBase::Execute()
   */
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        Group * domain ) override
  {
    GEOSX_UNUSED_VAR( cycleNumber );
    GEOSX_UNUSED_VAR( eventCounter );
    GEOSX_UNUSED_VAR( eventProgress );
    // std::function defines the == and =! comparable against nullptr_t to check the
    //  function pointer is actually assigned (an error would be thrown on the call attempt even so)
    GEOSX_ERROR_IF( m_buffer_call == nullptr,
                    "History collection buffer retrieval function is unassigned, did you declare a related TimeHistoryOutput event?" );
    // using GEOSX_ERROR_IF_EQ causes type issues since the values are used in iostreams
    buffer_unit_type * buffer = m_buffer_call();
    Collect( domain, time_n, dt, buffer );

    int rank = MpiWrapper::Comm_rank();
    if( rank == 0 && m_time_buffer_call )
    {
      buffer_unit_type * time_buffer = m_time_buffer_call();
      memcpy( time_buffer, &time_n, sizeof(decltype(time_n)) );
    }
  }

  /**
   * @brief Register a callback that gives the current head of the time history data buffer.
   * @param buffer_call A functional that when invoked returns a pointer to the head of a buffer at least large enough to
   *                    serialize one timestep of history data into.
   * @note This is typically meant to callback to BufferedHistoryIO::GetBufferHead( )
   */
  void RegisterBufferCall( std::function< buffer_unit_type *() > buffer_call )
  {
    m_buffer_call = buffer_call;
  }

  /**
   * @brief Get a metadata object relating the the Time variable itself.
   * @return A HistroyMetadata object describing the Time variable.
   */
  HistoryMetadata GetTimeMetadata( ) const
  {
    return HistoryMetadata( "Time", 1, std::type_index( typeid(real64)));
  }

  /**
   * @brief Register a callback that gives the current head of the time data buffer.
   * @param time_buffer_call A functional that when invoked returns a pointer to the head of a buffer at least large enough to
   *                    serialize one instance of the Time variable into.
   * @note This is typically meant to callback to BufferedHistoryIO::GetBufferHead( )
   */
  void RegisterTimeBufferCall( std::function< buffer_unit_type *() > time_buffer_call )
  {
    m_time_buffer_call = time_buffer_call;
  }

private:
  /// Callbacks to get the current time buffer head to write time data into
  std::function< buffer_unit_type *() > m_time_buffer_call;
  /// Callbacks to get the current buffer head to write history data into
  std::function< buffer_unit_type *() > m_buffer_call;
};

}

#endif