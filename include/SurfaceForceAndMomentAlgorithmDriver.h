/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SurfaceForceAndMomentAlgorithmDriver_h
#define SurfaceForceAndMomentAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class Realm;

class SurfaceForceAndMomentAlgorithmDriver : public AlgorithmDriver
{
public:

  SurfaceForceAndMomentAlgorithmDriver(
    Realm &realm,
    const int frequency);
  ~SurfaceForceAndMomentAlgorithmDriver();

  const int frequency_;

  std::vector<Algorithm *> algVec_;

  void execute();

  void zero_fields();
  void parallel_assemble_area();
  void parallel_assemble_fields();

  // deal with possible disparate frequencies (last one in wins)
  int get_frequency() {return frequency_;}
};
  

} // namespace nalu
} // namespace Sierra

#endif
