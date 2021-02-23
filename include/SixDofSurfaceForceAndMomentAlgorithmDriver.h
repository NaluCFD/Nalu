/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SixDofSurfaceForceAndMomentAlgorithmDriver_h
#define SixDofSurfaceForceAndMomentAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class Realm;

class SixDofSurfaceForceAndMomentAlgorithmDriver : public AlgorithmDriver
{
public:

  SixDofSurfaceForceAndMomentAlgorithmDriver(
    Realm &realm);
  ~SixDofSurfaceForceAndMomentAlgorithmDriver();

  std::vector<Algorithm *> algVec_;

  void execute();

  void zero_fields();
  void parallel_assemble_area();
  void parallel_assemble_fields();

};
  

} // namespace nalu
} // namespace Sierra

#endif
