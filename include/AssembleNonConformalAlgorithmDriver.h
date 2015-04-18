/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNonConformalAlgorithmDriver_h
#define AssembleNonConformalAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <FieldTypeDef.h>

namespace stk{
namespace mesh{
class FieldBase;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleNonConformalAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleNonConformalAlgorithmDriver(
    const Realm &realm,
    stk::mesh::FieldBase *ncNormalFlux,
    stk::mesh::FieldBase *ncPenalty,
    ScalarFieldType *ncArea,
    const int fluxFieldSize);
  ~AssembleNonConformalAlgorithmDriver();

  void pre_work();
  void post_work();

  // fields on which this AlgDriver operates
  stk::mesh::FieldBase *ncNormalFlux_;
  stk::mesh::FieldBase *ncPenalty_;
  stk::mesh::FieldBase *ncArea_;
  
  const unsigned fluxFieldSize_;
};
  

} // namespace nalu
} // namespace Sierra

#endif
