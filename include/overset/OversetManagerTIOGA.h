/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef OVERSETMANAGERTIOGA_H
#define OVERSETMANAGERTIOGA_H

#include "overset/OversetManager.h"
#include "overset/TiogaSTKIface.h"

#include <stk_mesh/base/FieldBase.hpp>

namespace sierra {
namespace nalu {

class Realm;
class OversetInfo;
struct OversetUserData;

/** Overset Connectivity Algorithm using TIOGA third-party library
 *
 *  This class is a thin Nalu-TIOGA wrapper to provide compatibility with Nalu's
 *  built-in STK based overset connectivity algorithm. The heavy lifting is done
 *  by the TiogaSTKIface class. Please refer to the documentation of that class
 *  for actual implementation details.
 */
class OversetManagerTIOGA : public OversetManager
{
public:
  OversetManagerTIOGA(Realm&, const OversetUserData&);

  virtual ~OversetManagerTIOGA();

  virtual void setup();

  virtual void initialize();

  /// Instance holding all the data from input files
  const OversetUserData& oversetUserData_;

  /// Tioga-STK interface instance that performs the necessary translation
  /// between TIOGA and STK data structures.
  tioga_nalu::TiogaSTKIface tiogaIface_;

  /// Flag tracking initialization phase for part registration
  bool isInit_{true};
};

}  // nalu
}  // sierra

#endif /* OVERSETMANAGERTIOGA_H */
