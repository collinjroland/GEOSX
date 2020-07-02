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

#include "TaskBase.hpp"

namespace geosx
{

using namespace dataRepository;

TaskBase::TaskBase( std::string const & name,
                    Group * const parent ):
  ExecutableGroup( name, parent )
{
}

TaskBase::~TaskBase()
{ }

TaskBase::CatalogInterface::CatalogType& TaskBase::GetCatalog()
{
  static TaskBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


void TaskBase::PostProcessInput()
{ }


} /* namespace */