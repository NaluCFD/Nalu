/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef GenericPropAlgorithm_h
#define GenericPropAlgorithm_h

#include <Algorithm.h>

namespace stk {
namespace mesh {
class FieldBase;
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;
class PropertyEvaluator;

class GenericPropAlgorithm : public Algorithm
{
public:

  GenericPropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop,
    PropertyEvaluator *propEvaluator);

  virtual ~GenericPropAlgorithm() {}

  virtual void execute();

  stk::mesh::FieldBase *prop_;
  PropertyEvaluator *propEvaluator_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
