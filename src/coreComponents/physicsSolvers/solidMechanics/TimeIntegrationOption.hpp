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

#pragma once

/**
 * @file TimeIntegrationOption.hpp
 */

#include "common/Enum.hpp"

namespace geosx
{

/**
 * @enum TimeIntegrationOption
 *
 * The options for time integration
 */
enum class TimeIntegrationOption : int
{
  QuasiStatic,    //!< QuasiStatic
  ImplicitDynamic,//!< ImplicitDynamic
  ExplicitDynamic //!< ExplicitDynamic
};

///@cond DO_NOT_DOCUMENT
namespace TimeIntegrationOptionStrings
{
char constexpr QuasiStatic[] = "QuasiStatic";
char constexpr ImplicitDynamic[] = "ImplicitDynamic";
char constexpr ExplicitDynamic[] = "ExplicitDynamic";
}
///@endcond

/**
 * @brief Type used to handle TimeIntegrationOption input values.
 */
struct TimeIntegrationOptionInput : public Enum< TimeIntegrationOption,
                                                 TimeIntegrationOptionStrings::QuasiStatic,
                                                 TimeIntegrationOptionStrings::ImplicitDynamic,
                                                 TimeIntegrationOptionStrings::ExplicitDynamic >
{ using Enum::Enum; };

} /// namespace geosx
