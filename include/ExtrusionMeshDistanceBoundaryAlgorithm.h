/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ExtrusionMeshDistanceBoundaryAlgorithm_h
#define ExtrusionMeshDistanceBoundaryAlgorithm_h

#include<Algorithm.h>

namespace sierra{
namespace nalu{

class Realm;

class ExtrusionMeshDistanceBoundaryAlgorithm : public Algorithm
{
public:

  ExtrusionMeshDistanceBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~ExtrusionMeshDistanceBoundaryAlgorithm() {}
  
  virtual void execute();

};

} // namespace nalu
} // namespace Sierra

#endif
