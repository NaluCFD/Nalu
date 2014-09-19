/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeGeometryBoundaryAlgorithm_h
#define ComputeGeometryBoundaryAlgorithm_h

#include<Algorithm.h>

namespace sierra{
namespace nalu{

class Realm;

class ComputeGeometryBoundaryAlgorithm : public Algorithm
{
public:

  ComputeGeometryBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~ComputeGeometryBoundaryAlgorithm() {}

  virtual void execute();
};

} // namespace nalu
} // namespace Sierra

#endif
