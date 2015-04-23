/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradAlgorithmDriver_h
#define AssembleNodalGradAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleNodalGradAlgorithmDriver(
    Realm &realm,
    const std::string & scalarQName,
    const std::string & dqdxName);
  ~AssembleNodalGradAlgorithmDriver();

  void pre_work();
  void post_work();

  const std::string scalarQName_;
  const std::string dqdxName_;
  
};
  

} // namespace nalu
} // namespace Sierra

#endif
