/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeGeometryInteriorAlgorithm_h
#define ComputeGeometryInteriorAlgorithm_h

#include<Algorithm.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeGeometryInteriorAlgorithm : public Algorithm
{
public:

  ComputeGeometryInteriorAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeGeometryInteriorAlgorithm();
  
  void execute();

  const bool assembleEdgeAreaVec_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
