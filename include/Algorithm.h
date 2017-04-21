/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Algorithm_h
#define Algorithm_h

#include <vector>

namespace stk {
namespace mesh {
class Part;
typedef std::vector<Part*> PartVector;
}
}
namespace sierra{
namespace nalu{

class Realm;
class MasterElement;
class SupplementalAlgorithm;

class Algorithm
{
public:

  // provide part
  Algorithm(
    Realm &realm,
    stk::mesh::Part *part);

  // provide part vector
  Algorithm(
    Realm &realm,
    stk::mesh::PartVector &partVec);

  virtual ~Algorithm();

  virtual void execute() = 0;

  virtual void pre_work() {}

  Realm &realm_;
  stk::mesh::PartVector partVec_;
  std::vector<SupplementalAlgorithm *> supplementalAlg_;
};

} // namespace nalu
} // namespace Sierra

#endif
