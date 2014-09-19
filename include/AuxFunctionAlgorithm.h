/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef AuxFunctionAlgorithm_h
#define AuxFunctionAlgorithm_h

#include <Algorithm.h>

#include <vector>
#include <stk_mesh/base/Types.hpp>

namespace stk{
namespace mesh{
class Part;
class FieldBase;
class Selector;

typedef std::vector< Part * > PartVector;
}
}

namespace sierra{
namespace nalu{

class AuxFunction;

class AuxFunctionAlgorithm : public Algorithm
{
public:

  AuxFunctionAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * field,
    AuxFunction * auxFunction,
    stk::mesh::EntityRank entityRank);

  virtual ~AuxFunctionAlgorithm();
  virtual void execute();

private:
  stk::mesh::FieldBase * field_;
  AuxFunction *auxFunction_;
  stk::mesh::EntityRank entityRank_;
  
private:
  // make this non-copyable
  AuxFunctionAlgorithm(const AuxFunctionAlgorithm & other);
  AuxFunctionAlgorithm & operator=(const AuxFunctionAlgorithm & other);
};

} // namespace nalu
} // namespace Sierra

#endif
