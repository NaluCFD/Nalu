/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef CopyFieldAlgorithm_h
#define CopyFieldAlgorithm_h

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

class CopyFieldAlgorithm : public Algorithm
{
public:

  CopyFieldAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    stk::mesh::FieldBase * fromField,
    stk::mesh::FieldBase * toField,
    const unsigned beginPos,
    const unsigned endPos,
    const stk::mesh::EntityRank entityRank);
  
  virtual ~CopyFieldAlgorithm() {}
  virtual void execute();

private:
  stk::mesh::FieldBase * fromField_;
  stk::mesh::FieldBase * toField_;

  const unsigned beginPos_;
  const unsigned endPos_;
  const stk::mesh::EntityRank entityRank_;
  
private:
  // make this non-copyable
  CopyFieldAlgorithm(const CopyFieldAlgorithm & other);
  CopyFieldAlgorithm & operator=(const CopyFieldAlgorithm & other);
};

} // namespace nalu
} // namespace Sierra

#endif
