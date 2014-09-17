/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotElemAlgorithm_h
#define ComputeMdotElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMdotElemAlgorithm : public Algorithm
{
public:

  ComputeMdotElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool assembleMdotToEdge);
  ~ComputeMdotElemAlgorithm();

  void execute();
  void assemble_edge_mdot();

  const bool assembleMdotToEdge_;

  // extract fields; nodal
  VectorFieldType *velocity_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  GenericFieldType *massFlowRate_;
  ScalarFieldType *edgeMassFlowRate_;

};

} // namespace nalu
} // namespace Sierra

#endif
