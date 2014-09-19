/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef InversePropAlgorithm_h
#define InversePropAlgorithm_h

#include <Algorithm.h>

namespace stk {
namespace mesh {
class FieldBase;
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class InversePropAlgorithm : public Algorithm
{
public:

  InversePropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * prop,
    stk::mesh::FieldBase * indVar,
    const double primary,
    const double secondary);
  
  virtual ~InversePropAlgorithm();
  
  virtual void execute();
  
  stk::mesh::FieldBase *prop_;
  stk::mesh::FieldBase *indVar_;
  const double primary_;
  const double secondary_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
