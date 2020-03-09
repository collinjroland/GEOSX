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

/*
 * @file SimpleGeometricObjectBase.hpp
 */

#ifndef SIMPLEGEOMETRICOBJECTS_H_
#define SIMPLEGEOMETRICOBJECTS_H_

//#include "common/Common.h"
#include "dataRepository/Group.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ObjectCatalog.hpp"

class Function;

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const geometricObjects("GeometricObjects");
}
}


class SimpleGeometricObjectBase : public dataRepository::Group
{
public:

  explicit SimpleGeometricObjectBase( std::string const & name,
                                      Group * const parent );

  virtual ~SimpleGeometricObjectBase();

  static string CatalogName() { return "SimpleGeometricObjectBase"; }

  virtual void PostProcessInput() override;

  virtual bool IsCoordInObjectSub( const R1Tensor& coord ) const = 0;

  bool IsCoordInObject( const R1Tensor& coord ) const
  {
    return IsCoordInObjectSub( coord ) ^ m_invertNodeset;
  }

  using CatalogInterface = dataRepository::CatalogInterface< SimpleGeometricObjectBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType& GetCatalog();

private:
  bool m_invertNodeset;
  integer m_invertNodesetInput;

  struct viewKeyStruct
  {
    static constexpr auto invertString = "invert";
    dataRepository::ViewKey invert = { "invert" };
  } viewKeys;

};


}
#endif /* SIMPLEGEOMETRICOBJECTS_H_ */
