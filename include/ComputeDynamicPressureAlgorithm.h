/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeDynamicPressureAlgorithm_h
#define ComputeDynamicPressureAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeDynamicPressureAlgorithm : public Algorithm
{
public:

  ComputeDynamicPressureAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool useShifted);
  ~ComputeDynamicPressureAlgorithm();

  void execute();

  const bool useShifted_;

  ScalarFieldType *density_;
  GenericFieldType *openMassFlowRate_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *dynamicPressure_;
};

} // namespace nalu
} // namespace Sierra

#endif
