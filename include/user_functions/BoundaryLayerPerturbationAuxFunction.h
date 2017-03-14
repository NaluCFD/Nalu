/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef BoundaryLayerPerturbationAuxFunction_h
#define BoundaryLayerPerturbationAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class BoundaryLayerPerturbationAuxFunction : public AuxFunction
{
public:

  BoundaryLayerPerturbationAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    const std::vector<double> &theParams);

  virtual ~BoundaryLayerPerturbationAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
private:
  double amplitude_;
  double kx_;
  double ky_;
  double thickness_;
  double uInf_;

};

} // namespace nalu
} // namespace Sierra

#endif
