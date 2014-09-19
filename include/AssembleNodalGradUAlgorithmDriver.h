/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradUAlgorithmDriver_h
#define AssembleNodalGradUAlgorithmDriver_h

#include<AlgorithmDriver.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradUAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleNodalGradUAlgorithmDriver(
    const Realm &realm);
  virtual ~AssembleNodalGradUAlgorithmDriver() {}

  virtual void pre_work();
  virtual void post_work();
};

} // namespace nalu
} // namespace Sierra

#endif
