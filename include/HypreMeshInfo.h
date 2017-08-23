/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef HYPREMESHINFO_H
#define HYPREMESHINFO_H

#include "FieldTypeDef.h"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Entity.hpp"

#include <map>

namespace sierra {
namespace nalu {

/** Track STK to Hypre mapping information
 *
 */
class HypreMeshInfo
{
public:
  ScalarIntFieldType* hypreGlobalId_;

  int iLower_{0};
  int iUpper_{-1};
};

}  // nalu
}  // sierra


#endif /* HYPREMESHINFO_H */
