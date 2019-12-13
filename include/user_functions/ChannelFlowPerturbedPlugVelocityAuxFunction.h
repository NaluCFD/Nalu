/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ChannelFlowPerturbedPlugVelocityAuxFunction_h
#define ChannelFlowPerturbedPlugVelocityAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class ChannelFlowPerturbedPlugVelocityAuxFunction : public AuxFunction
{
public:

  ChannelFlowPerturbedPlugVelocityAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    const std::vector<double> &params);

  virtual ~ChannelFlowPerturbedPlugVelocityAuxFunction() {}
  
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
  double u_m;
  double Lx_;
  double Ly_;
  double Lz_;
  int streamwiseD_;
  int spanwiseD_;
  int normalD_;
  const double pi_;
};

} // namespace nalu
} // namespace Sierra

#endif
