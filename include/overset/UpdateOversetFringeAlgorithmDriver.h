/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef UPDATEOVERSETFRINGEALGORITHMDRIVER_H
#define UPDATEOVERSETFRINGEALGORITHMDRIVER_H

#include "AlgorithmDriver.h"


#include <memory>
#include <vector>

namespace stk {
namespace mesh {
class FieldBase;
}
}

namespace sierra {
namespace nalu {

class Realm;

struct OversetFieldData
{
  OversetFieldData(stk::mesh::FieldBase* field, int sizeRow=1, int sizeCol=1)
    : field_(field),
      sizeRow_(sizeRow),
      sizeCol_(sizeCol)
  {}

  stk::mesh::FieldBase* field_;
  int sizeRow_;
  int sizeCol_;
};

class UpdateOversetFringeAlgorithmDriver : public AlgorithmDriver
{
public:
  UpdateOversetFringeAlgorithmDriver(Realm& realm);

  virtual ~UpdateOversetFringeAlgorithmDriver();

  virtual void pre_work();

  std::vector<std::unique_ptr<OversetFieldData>> fields_;
};

}  // nalu
}  // sierra

#endif /* UPDATEOVERSETFRINGEALGORITHMDRIVER_H */
