/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleGasDynamicsSymmetryAlgorithm_h
#define AssembleGasDynamicsSymmetryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class AssembleGasDynamicsSymmetryAlgorithm : public Algorithm
{
public:

  AssembleGasDynamicsSymmetryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *pressure,
    GenericFieldType *rhsGasDyn);
  virtual ~AssembleGasDynamicsSymmetryAlgorithm() {}

  virtual void execute();

  ScalarFieldType *pressure_;
  GenericFieldType *rhsGasDyn_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
