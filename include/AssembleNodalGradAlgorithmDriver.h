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
    const std::string & dqdxName,
    const std::string & areaWeightName = "na",
    const bool areaWeight = false);
  ~AssembleNodalGradAlgorithmDriver();

  void pre_work();
  void post_work();
  void normalize_by_area();

  const std::string scalarQName_;
  const std::string dqdxName_;
  const std::string areaWeightName_;  
  const bool areaWeight_;  
};
  

} // namespace nalu
} // namespace Sierra

#endif
