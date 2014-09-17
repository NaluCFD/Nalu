/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TemperaturePropAlgorithm_h
#define TemperaturePropAlgorithm_h

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

class TemperaturePropAlgorithm : public Algorithm
{
public:

  TemperaturePropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop,
    stk::mesh::FieldBase * temperature,
    PropertyEvaluator *propEvaluator);

  virtual ~TemperaturePropAlgorithm() {}

  virtual void execute();

  stk::mesh::FieldBase *prop_;
  stk::mesh::FieldBase *temperature_;
  PropertyEvaluator *propEvaluator_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
