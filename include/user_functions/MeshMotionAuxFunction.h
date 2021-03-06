/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MeshMotionAuxFunction_h
#define MeshMotionAuxFunction_h

#include <AuxFunction.h>

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class Realm;

class MeshMotionAuxFunction : public AuxFunction
{
public:

  MeshMotionAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    std::vector<std::string> theStringParams,
    Realm &realm);

  virtual ~MeshMotionAuxFunction();
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;

  void setup(const double time);
  void cross_product(double *c, double *u) const;

private:
  std::vector<double>* omegaMM_;
  std::vector<double>* centroidMM_;
  std::vector<double>* velMM_;
  std::vector<double>* dispMM_;
};

} // namespace nalu
} // namespace Sierra

#endif
