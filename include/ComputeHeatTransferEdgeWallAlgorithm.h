/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeHeatTransferEdgeWallAlgorithm_h
#define ComputeHeatTransferEdgeWallAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeHeatTransferEdgeWallAlgorithm : public Algorithm
{
public:

  ComputeHeatTransferEdgeWallAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeHeatTransferEdgeWallAlgorithm();

  void execute();

  ScalarFieldType *temperature_;
  VectorFieldType *dhdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *thermalCond_;
  ScalarFieldType *specificHeat_;
  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *assembledWallArea_;
  ScalarFieldType *referenceTemperature_;
  ScalarFieldType *heatTransferCoefficient_;
  ScalarFieldType *normalHeatFlux_;
  ScalarFieldType *robinCouplingParameter_;

  double compute_coupling_parameter(const double & kappa,
                                    const double & h,
                                    const double & chi);  
};

} // namespace nalu
} // namespace Sierra

#endif
