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
 * @file Enum.hpp
 *
 * Collection of utilities to facilitate I/O of enumeration types.
 * Provides a macro definition that allows associating string names
 * with enumeration constants and a set of functions that make use
 * of these strings, like stream insertion/extraction operators.
 */

#ifndef GEOSX_COMMON_ENUM_HPP_
#define GEOSX_COMMON_ENUM_HPP_

#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "common/Logger.hpp"

#include <iostream>
#include <sstream>
#include <array>
#include <type_traits>
#include <algorithm>

namespace geosx
{

namespace internal
{

/**
 * @brief Simple compile-time variadic function that counts its arguments.
 * @tparam ARGS variadic pack of argument types
 * @return the number of arguments passed
 */
template< typename ... ARGS >
constexpr int countArgs( ARGS ... )
{
  return sizeof...( ARGS );
}

/**
 * @brief Provides some useful properties (like string representation) for enumerations.
 * @tparam T the enumeration type
 * @note This template is meant to be specialized using ENUM_STRINGS macro.
 */
template< typename T >
struct EnumStringsBase
{
  /// Flag used to distinguish enumeration types that have strings.
  static constexpr bool defined = false;
};

}

/**
 * @brief Associate a list of string names with enumeration values.
 * @param ENUM the enumeration type
 * @param ... list of names (C-string literals)
 *
 * Conditions (not enforced but won't work correctly if violated):
 *  - the macro must be called in geosx namespace (not any nested namespace therein)
 *  - the number and order of string arguments passed must match the number of enum values
 *  - enumeration constants must not have custom values assigned
 *
 * After the macro has been called, template instantiation EnumStrings<ENUM>
 * may be used to get access to strings at runtime. While not strictly necessary,
 * it is recommended that macro call immediately follows the enum definition.
 */
#define ENUM_STRINGS( ENUM, ... )                                     \
namespace internal                                                    \
{                                                                     \
template<>                                                            \
struct EnumStringsBase< ENUM >                                        \
{                                                                     \
  static_assert( std::is_enum< ENUM >::value, "Not an enumeration" ); \
  static constexpr bool defined = true;                               \
                                                                      \
  static auto const & get()                                           \
  {                                                                   \
    constexpr int N = internal::countArgs( __VA_ARGS__ );             \
    static constexpr std::array< char const *, N > ss{ __VA_ARGS__ }; \
    return ss;                                                        \
  }                                                                   \
};                                                                    \
}

/**
 * @brief Provides enum <-> string conversion facilities.
 * @tparam ENUM the enumeration type
 */
template< typename ENUM >
struct EnumStrings : private internal::EnumStringsBase< ENUM >
{
private:

  /// Alias for private base type containing the strings
  using Base = internal::EnumStringsBase< ENUM >;

public:

  /// Alias for the enumeration type
  using enum_type = ENUM;

  /// Alias for enum's underlying fundamental type
  using base_type = std::underlying_type_t< ENUM >;

  using Base::get;

  /**
   * @brief Get a list of valid options as a delimited string.
   * @param delim delimiter (defaults to single space)
   * @return the string containing all valid strings for this type
   */
  static string concat( string const & delim = " " )
  {
    auto const & strings = get();
    return ::geosx::stringutilities::strjoin( begin( strings ), end( strings ), delim );
  }

  /**
   * @brief Convert enum to string.
   * @param e the enum value to convert
   * @return the corresponding string
   *
   * An error is raised if enum's numerical value is greater of equal than the number of strings.
   */
  static string toString( enum_type const & e )
  {
    auto const & strings = get();
    base_type const index = static_cast< base_type >( e );
    GEOSX_ERROR_IF_GE_MSG( index, LvArray::integerConversion< base_type >( strings.size() ),
                           "Unrecognized value of " << TypeName< ENUM >::brief() << ": " << index << ".\n" <<
                           "Valid range is 0.." << strings.size() - 1 );
    return strings[ index ];
  }

  /**
   * @brief Convert string to enum
   * @param s the string to convert
   * @return the corresponding enum value
   */
  static enum_type fromString( string const & s )
  {
    auto const & names = get();
    auto const it = std::find( begin( names ), end( names ), s );
    GEOSX_ERROR_IF( it == names.end(),
                    "'" << s << "' is not a recognized value of " << TypeName< enum_type >::brief() << ".\n" <<
                        "Valid options are: \n  " << concat( "\n  " ) );
    enum_type const e = static_cast< enum_type >( LvArray::integerConversion< base_type >( std::distance( names.begin(), it ) ) );
    return e;
  }
};

/**
 * @brief Specialization of TypeRegex for enumeration types with strings attached (pun intended).
 * @tparam ENUM the type of enumeration
 */
template< typename ENUM >
struct TypeRegex< ENUM, std::enable_if_t< internal::EnumStringsBase< ENUM >::defined > >
{
  /**
   * @brief @return Regex for validating enumeration inputs for @p ENUM type.
   */
  static std::string get()
  {
    return EnumStrings< ENUM >::concat( "|" );
  }
};

/**
 * @brief Stream insertion operator for enumeration types with strings.
 * @tparam ENUM the enumeration type
 * @param os the output stream
 * @param e the enumeration value
 * @return @p os
 */
template< typename ENUM >
std::enable_if_t< internal::EnumStringsBase< ENUM >::defined, std::ostream & >
operator<<( std::ostream & os, ENUM const & e )
{
  os << EnumStrings< ENUM >::toString( e );
  return os;
}

/**
 * @brief Stream extraction operator for enumeration types with strings.
 * @tparam ENUM the enumeration type
 * @param is the input stream
 * @param e the enumeration value
 * @return @p is
 */
template< typename ENUM >
std::enable_if_t< internal::EnumStringsBase< ENUM >::defined, std::istream & >
operator>>( std::istream & is, ENUM & e )
{
  string s;
  is >> s;
  e = EnumStrings< ENUM >::fromString( s );
  return is;
}

} // namespace geosx

#endif //GEOSX_COMMON_ENUM_HPP_
