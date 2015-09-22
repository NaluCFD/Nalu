/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PropertyEvaluator_h
#define PropertyEvaluator_h

#include <stk_mesh/base/Entity.hpp>

#include <vector>

namespace sierra{
namespace nalu{

class PropertyEvaluator
{
public:

  PropertyEvaluator() {}
  virtual ~PropertyEvaluator() {}
  
  virtual double execute(
    double *indVarList,
    stk::mesh::Entity node = stk::mesh::Entity()) = 0;
  
};

} // namespace nalu
} // namespace Sierra

#endif
