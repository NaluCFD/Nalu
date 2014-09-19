/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleWallHeatTransferAlgorithmDriver_h
#define AssembleWallHeatTransferAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleWallHeatTransferAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleWallHeatTransferAlgorithmDriver(
    const Realm &realm);
  ~AssembleWallHeatTransferAlgorithmDriver();

  void pre_work();
  void post_work();

  ScalarFieldType *assembledWallArea_;
  ScalarFieldType *referenceTemperature_;
  ScalarFieldType *heatTransferCoefficient_;
  ScalarFieldType *temperature_;

};
  

} // namespace nalu
} // namespace Sierra

#endif
