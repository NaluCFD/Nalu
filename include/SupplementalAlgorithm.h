/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SupplementalAlgorithm_h
#define SupplementalAlgorithm_h

#include <vector>

#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { struct Entity; } }

namespace sierra{
namespace nalu{

class Realm;

class SupplementalAlgorithm
{
public:
  
  SupplementalAlgorithm(
    Realm &realm);
  
  virtual ~SupplementalAlgorithm() {}

  virtual void setup() {}

  virtual void elem_execute(
    const int &numScvIntPoints,
    const int &numScsIntPoints,
    double *lhs,
    double *rhs,
    stk::mesh::Entity elem) = 0;
  
  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node) = 0;
  
  Realm &realm_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
