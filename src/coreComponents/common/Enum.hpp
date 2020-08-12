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
 */
#ifndef GEOSX_COMMON_ENUM_HPP_
#define GEOSX_COMMON_ENUM_HPP_

#include "Logger.hpp"

#include <iostream>
#include <array>
#include <type_traits>

namespace geosx
{

/**
 * @brief Generic printable/parsable enumeration class.
 * @tparam T underlying c++ type (strong/weak enum)
 * @tparam SS list of compile-time strings for recognized enum values
 *
 * @warning The enum type used must not have non-default (user-assigned) values for its constants.
 */
template< typename T, char const * ... SS >
struct Enum
{
  static_assert( std::is_enum< T >::value, "T must be an enumeration type" );
  static_assert( sizeof...( SS ) > 0, "Must provide at least one recognized string label" );

  using underlying_type = std::underlying_type_t< T >;

  GEOSX_HOST_DEVICE
  Enum() = default;

  GEOSX_HOST_DEVICE
  Enum( T const v )
    : m_val( v )
  { }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  operator T() const
  {
    return m_val;
  }

  static string validOptions( string const & delim )
  {
    std::ostringstream os;
    os << *s_labels.begin();
    for( auto it = std::next( std::begin( s_labels ) ); it != std::end( s_labels ); ++it )
    {
      os << delim << *it;
    }
    return os.str();
  }

  static string typeRegex()
  {
    return validOptions( "|" );
  }

private:

  static constexpr std::array< char const *, sizeof...( SS ) > s_labels{ SS... };

  static string typeName()
  {
    string const full_name = LvArray::system::demangle( typeid(T).name() );
    string::size_type const pos = full_name.find_last_of( "::" );
    return ( pos == string::npos ) ? full_name : full_name.substr( pos );
  }

  friend std::ostream & operator<<( std::ostream & os, Enum const & e )
  {
    underlying_type const index = static_cast< underlying_type >( e.m_val );
    GEOSX_ERROR_IF_GE_MSG( index, static_cast< underlying_type >( s_labels.size() ),
                           "Unrecognized value of " << typeName() << ": " << index << ".\n"
                           "Valid range is 0.." << s_labels.size()-1 );
    os << Enum::s_labels[ index ];
    return os;
  }

  friend std::istream & operator>>( std::istream & is, Enum & e )
  {
    std::string s;
    is >> s;
    auto const it = std::find( Enum::s_labels.begin(), Enum::s_labels.end(), s );
    GEOSX_ERROR_IF( it == Enum::s_labels.end(),
                    "'" << s << "' is not a recognized value of " << typeName() << ".\n"
                    "Valid options are: \n  " << validOptions( "\n  " ) );
    e.m_val = static_cast< T >( std::distance( Enum::s_labels.begin(), it ) );
    return is;
  }

  T m_val;
};

template< typename T, char const * ... SS >
constexpr std::array< char const *, sizeof...( SS ) > Enum< T, SS... >::s_labels;

} // namespace geosx

#endif //GEOSX_COMMON_ENUM_HPP_
