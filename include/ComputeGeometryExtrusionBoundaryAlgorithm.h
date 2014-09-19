/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeGeometryExtrusionBoundaryAlgorithm_h
#define ComputeGeometryExtrusionBoundaryAlgorithm_h

#include<Algorithm.h>

namespace sierra{
namespace nalu{

class Realm;

class ComputeGeometryExtrusionBoundaryAlgorithm : public Algorithm
{
public:

  ComputeGeometryExtrusionBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~ComputeGeometryExtrusionBoundaryAlgorithm() {}

  virtual void execute();
};

} // namespace nalu
} // namespace Sierra

#endif
