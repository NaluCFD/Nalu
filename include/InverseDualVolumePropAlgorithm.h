/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef InverseDualVolumePropAlgorithm_h
#define InverseDualVolumePropAlgorithm_h

#include <Algorithm.h>
#include <FieldTypeDef.h>

namespace stk {
namespace mesh {
class FieldBase;
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class InverseDualVolumePropAlgorithm : public Algorithm
{
public:

  InverseDualVolumePropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop);
  
  virtual ~InverseDualVolumePropAlgorithm();
  
  virtual void execute();
  
  stk::mesh::FieldBase *prop_;
  ScalarFieldType *dualNodalVolume_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
