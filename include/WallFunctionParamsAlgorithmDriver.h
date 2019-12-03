/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef WallFunctionParamsAlgorithmDriver_h
#define WallFunctionParamsAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class WallFunctionParamsAlgorithmDriver : public AlgorithmDriver
{
public:

  WallFunctionParamsAlgorithmDriver(
    Realm &realm);
  ~WallFunctionParamsAlgorithmDriver();

  void pre_work();
  void post_work();

  ScalarFieldType *assembledWallArea_;
  ScalarFieldType *assembledWallNormalDistance_;
};
  

} // namespace nalu
} // namespace Sierra

#endif
