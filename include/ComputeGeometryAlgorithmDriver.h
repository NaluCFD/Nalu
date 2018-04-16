/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeGeometryAlgorithmDriver_h
#define ComputeGeometryAlgorithmDriver_h

#include<AlgorithmDriver.h>

namespace sierra{
namespace nalu{

class Realm;

class ComputeGeometryAlgorithmDriver : public AlgorithmDriver
{
public:

  ComputeGeometryAlgorithmDriver(
    Realm &realm);
  virtual ~ComputeGeometryAlgorithmDriver() {}

  void pre_work();
  void post_work();
  void check_jacobians();
};
  

} // namespace nalu
} // namespace Sierra

#endif
