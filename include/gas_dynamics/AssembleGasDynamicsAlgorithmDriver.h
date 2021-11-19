/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef AssembleGasDynamicsAlgorithmDriver_h
#define AssembleGasDynamicsAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;

class AssembleGasDynamicsAlgorithmDriver : public AlgorithmDriver
{
public:

  AssembleGasDynamicsAlgorithmDriver(
    Realm &realm);
  ~AssembleGasDynamicsAlgorithmDriver();

  void pre_work();
  void post_work();
};
  

} // namespace nalu
} // namespace Sierra

#endif
