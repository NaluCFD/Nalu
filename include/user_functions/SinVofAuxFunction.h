/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SinVofAuxFunction_h
#define SinVofAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class SinVofAuxFunction : public AuxFunction
{
public:

  SinVofAuxFunction(
    const std::vector<double> &params);

  virtual ~SinVofAuxFunction() {}
  
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
  int coordDirection_;
  int heightDirection_;
  double lambda_;
  double aNot_;
  double pi_;
};

} // namespace nalu
} // namespace Sierra

#endif
