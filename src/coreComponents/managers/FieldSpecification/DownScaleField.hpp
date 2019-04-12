
/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * DownScaleright (c) 2019, Lawrence Livermore National Security, LLC.
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

#ifndef DOWNSCALEFIELD_HPP_
#define DOWNSCALEFIELD_HPP_



#include <string>
#include <limits>

#include <managers/Tasks/TaskBase.hpp>

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}

class DownScaleField : public TaskBase
{
public:

  explicit DownScaleField( std::string const & name,
                      ManagedGroup * const parent );

  DownScaleField( DownScaleField && ) = default;

  virtual ~DownScaleField() override;

  DownScaleField() = delete;
  DownScaleField( DownScaleField const & ) = delete;
  DownScaleField& operator=( DownScaleField const & ) = delete;
  DownScaleField& operator=( DownScaleField&& ) = delete;

  static string CatalogName() { return "DownScaleField"; }
  void Execute( real64 const time_n,
                real64 const dt,
                integer const cycleNumber,
                integer const eventCounter,
                real64 const eventProgress,
                ManagedGroup * domain ) override;

  struct viewKeyStruct {
    static constexpr auto fieldNameString = "fieldName";
  };
private:
  /// Name of the field
  string m_fieldName;

  /// Regions on which the copy will be done
  string_array m_targetRegions;
};
} /* namespace geosx */


#endif /* COPYFIELD_HPP_ */
