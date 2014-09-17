/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TablePropAlgorithm_h
#define TablePropAlgorithm_h

#include <Algorithm.h>

#include <string>
#include <vector>

namespace stk {
namespace mesh {
class FieldBase;
class Part;
class MetaData;
}
}

namespace sierra{
namespace nalu{

class Realm;
class PropertyEvaluator;
class StateTable;

class TablePropAlgorithm : public Algorithm
{
public:

  TablePropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    StateTable *stateTable,
    stk::mesh::FieldBase * prop,
    std::string tablePropName,
    std::vector<std::string> &indVarNameVec,
    std::vector<std::string> &indVarTableNameVec,
    const size_t cIndex,
    const stk::mesh::MetaData &meta_data);

  virtual ~TablePropAlgorithm() {}

  virtual void execute();
  size_t find_index(const double z, size_t iMin, size_t iMax);
  double interp_property(const double *z);

  stk::mesh::FieldBase *prop_;
  std::string tablePropName_;
  std::vector<std::string> indVarTableNameVec_;
  const size_t cIndex_;
  const size_t indVarSize_;
  size_t iMin_;
  size_t iMax_;

  std::vector<std::vector<double > > table_;
  std::vector<stk::mesh::FieldBase *> indVar_;
  std::vector<double *> workIndVar_;
  std::vector<double> workZ_;

};

} // namespace nalu
} // namespace Sierra

#endif
