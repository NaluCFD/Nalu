/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotEdgeContactAlgorithm_h
#define ComputeMdotEdgeContactAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
class FieldBase;
}
}

namespace sierra{
namespace nalu{

class Realm;

class ComputeMdotEdgeContactAlgorithm : public Algorithm
{
public:

  ComputeMdotEdgeContactAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~ComputeMdotEdgeContactAlgorithm() {}
  virtual void execute();

  const bool meshMotion_;

  VectorFieldType *meshVelocity_;
  VectorFieldType *velocity_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  ScalarFieldType *haloMdot_;
  
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
