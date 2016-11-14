/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef WindEnergyAuxFunction_h
#define WindEnergyAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class Realm;
class TanhFunction;

class WindEnergyAuxFunction : public AuxFunction
{
public:

  WindEnergyAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    std::vector<double> theParams,
    Realm &realm);

  virtual ~WindEnergyAuxFunction();
  
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

private:
  double centroidX_;
  double centroidY_;
  double omega_;
  double omegaBlend_;
  TanhFunction *tanhFunction_;
};

} // namespace nalu
} // namespace Sierra

#endif
