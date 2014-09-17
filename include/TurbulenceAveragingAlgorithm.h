/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbulenceAveragingAlgorithm_h
#define TurbulenceAveragingAlgorithm_h

#include <Algorithm.h>
#include <FieldTypeDef.h>

// c++
#include <vector>
#include <utility>

namespace stk {
namespace mesh {
class Part;
class FieldBase;
}
}

namespace sierra{
namespace nalu{

class Realm;

class TurbulenceAveragingAlgorithm : public Algorithm
{
public:

  TurbulenceAveragingAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);

  virtual ~TurbulenceAveragingAlgorithm() {}
  virtual void execute();
  
  std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > favreFieldVecPair_;
  std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > reynoldsFieldVecPair_;
  std::vector<unsigned> favreFieldSize_;
  std::vector<unsigned> reynoldsFieldSize_;

};

} // namespace nalu
} // namespace Sierra

#endif
