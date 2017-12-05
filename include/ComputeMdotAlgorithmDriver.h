/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotAlgorithmDriver_h
#define ComputeMdotAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>

namespace sierra{
namespace nalu{

class Realm;
class SolutionOptions;

class ComputeMdotAlgorithmDriver : public AlgorithmDriver
{
public:

  ComputeMdotAlgorithmDriver(
    Realm &realm);

  ~ComputeMdotAlgorithmDriver();

  double compute_accumulation();
  void correct_open_mdot(const double finalCorrection);
  void provide_output();

  SolutionOptions &solnOpts_;
  bool hasMass_;
  bool lumpedMass_;

  void pre_work();
  void post_work();
};
  

} // namespace nalu
} // namespace Sierra

#endif
