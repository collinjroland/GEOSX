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

#ifndef SRC_COMPONENTS_CORE_MESHUTILITIES_MESHTRANSLATOR_HPP_
#define SRC_COMPONENTS_CORE_MESHUTILITIES_MESHTRANSLATOR_HPP_

namespace geosx
{

namespace meshTranslator
{

namespace VTKToInternal
{
  static const integer tetra[4] = { 0, 1, 2, 3 };
  static const integer hex[8]   = { 0, 1, 3, 2, 4, 5, 7, 6 }
}

}
}

#endif
