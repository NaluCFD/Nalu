/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ConstantPropertyEvaluator_h
#define ConstantPropertyEvaluator_h

#include <property_evaluator/PropertyEvaluator.h>

namespace stk {
namespace mesh {
struct Entity;
}
}

namespace sierra{
namespace nalu{

class ConstantPropertyEvaluator : public PropertyEvaluator
{
 public:

  ConstantPropertyEvaluator(const double & value);
  virtual ~ConstantPropertyEvaluator();

  double execute(double *indVarList,
                 stk::mesh::Entity node);

  double value_;

};


} // namespace nalu
} // namespace Sierra

#endif
