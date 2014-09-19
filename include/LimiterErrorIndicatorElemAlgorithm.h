/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LimiterErrorIndicatorElemAlgorithm_h
#define LimiterErrorIndicatorElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class LimiterErrorIndicatorElemAlgorithm : public Algorithm
{
public:

  LimiterErrorIndicatorElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~LimiterErrorIndicatorElemAlgorithm();

  void execute();
  double van_leer(
    const double &dqm,
    const double &dqp,
    const double &small);

  // extract fields; nodal
  VectorFieldType *velocity_;
  VectorFieldType *coordinates_;
  GenericFieldType *dudx_;
  GenericFieldType *LimiterEI_;

};

} // namespace nalu
} // namespace Sierra

#endif
