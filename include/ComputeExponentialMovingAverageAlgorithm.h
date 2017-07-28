/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ComputeExponentialMovingAverageAlgorithm_h
#define ComputeExponentialMovingAverageAlgorithm_h

#include <Algorithm.h>
#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class ExponentialMovingAverager
{
public:
  ExponentialMovingAverager(double timeScale);

  double compute_updated_average(double oldAvg, double newVal);
  void compute_and_set_alpha(double dt);

  void init_state(bool init);
private:
  double timeScale_;
  bool isInit_;
  double alpha_;
};


class ComputeExponentialMovingAverageAlgorithm : public Algorithm
{
public:
  ComputeExponentialMovingAverageAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    stk::mesh::FieldBase* field,
    double timeScale);

  virtual void execute();
  
  stk::mesh::FieldBase* field_;
  ExponentialMovingAverager averager_;

  stk::mesh::FieldBase* avgField_;
};

} // namespace nalu
} // namespace Sierra

#endif
