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
    Realm &realm,
    const std::string dudxName);
  virtual ~AssembleNodalGradUAlgorithmDriver() {}

  virtual void pre_work();
  virtual void post_work();

  const std::string dudxName_;
};

} // namespace nalu
} // namespace Sierra

#endif
