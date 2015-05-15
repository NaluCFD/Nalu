/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotNonConformalAlgorithm_h
#define ComputeMdotNonConformalAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMdotNonConformalAlgorithm : public Algorithm
{
public:

  ComputeMdotNonConformalAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *pressure,
    ScalarFieldType *ncNormalFlux,
    ScalarFieldType *ncPenalty);

  ~ComputeMdotNonConformalAlgorithm();

  void execute();

  ScalarFieldType *pressure_;
  ScalarFieldType *ncNormalFlux_;
  ScalarFieldType *ncPenalty_;
  
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *ncMassFlowRate_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
