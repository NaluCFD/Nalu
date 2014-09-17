/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradUEdgeBoundaryAlgorithm_h
#define AssembleNodalGradUEdgeBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradUEdgeBoundaryAlgorithm : public Algorithm
{
public:

  AssembleNodalGradUEdgeBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *velocity,
    GenericFieldType *dudx);
  virtual ~AssembleNodalGradUEdgeBoundaryAlgorithm() {}

  virtual void execute();

  VectorFieldType *velocity_;
  GenericFieldType *dudx_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
