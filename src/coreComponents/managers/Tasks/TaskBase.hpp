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

#ifndef TASKBASE_HPP_
#define TASKBASE_HPP_

#include <string>
#include <limits>

#include "dataRepository/ExecutableGroup.hpp"
#include "common/DataTypes.hpp"
namespace geosx
{

class TaskBase : public ExecutableGroup
{
public:

  explicit TaskBase( std::string const & name,
                     Group * const parent );

  TaskBase( TaskBase && ) = default;

  virtual ~TaskBase() override;

  TaskBase() = delete;
  TaskBase( TaskBase const & ) = delete;
  TaskBase& operator=( TaskBase const & ) = delete;
  TaskBase& operator=( TaskBase&& ) = delete;

  static string CatalogName() { return "TaskBase"; }

  using CatalogInterface = dataRepository::CatalogInterface< TaskBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  void PostProcessInput() override;


};

} /* namespace */

#endif /* TASKBASE_HPP_ */
