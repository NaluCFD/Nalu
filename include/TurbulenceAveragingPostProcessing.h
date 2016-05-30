/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbulenceAveragingPostProcessing_h
#define TurbulenceAveragingPostProcessing_h

#include <NaluParsing.h>

#include <string>
#include <vector>
#include <utility>

// stk forwards
namespace stk {
  namespace mesh {
    class BulkData;
    class FieldBase;
    class MetaData;
    class Part;
    typedef std::vector<Part*> PartVector;
    class Selector;
  }
}

namespace sierra{
namespace nalu{

class Realm;
class AveragingInfo;

class TurbulenceAveragingPostProcessing
{
public:
  
  TurbulenceAveragingPostProcessing(
    Realm &realm,
    const YAML::Node &node);
  ~TurbulenceAveragingPostProcessing();
  
  // load all of the options
  void load(
    const YAML::Node & node);

  // setup nodal field registration; parts, fields, etc
  void setup();

  void register_field_from_primitive(
    const std::string primitiveName,
    const std::string averagedName,
    stk::mesh::MetaData &metaData,
    stk::mesh::Part *part);

  void construct_pair(
    const std::string primitiveName,
    const std::string averagedName,
    std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > &fieldVecPair,
    std::vector<unsigned> &fieldSizeVec_,
    stk::mesh::MetaData &metaData);

  void register_field(
    const std::string fieldName,
    const int fieldSize,
    stk::mesh::MetaData &metaData,
    stk::mesh::Part *targetPart);

  void review( 
    const AveragingInfo *avInfo);

  // populate nodal field and output norms (if appropriate)
  void execute();

  // compute tke and stress for each type of operation
  void compute_tke(
    const bool isReynolds,
    const std::string &averageBlockName,
    stk::mesh::Selector s_all_nodes);

  void compute_reynolds_stress(
    const std::string &averageBlockName,
    const double &oldTimeFilter,
    const double &zeroCurrent,
    const double &dt,
    stk::mesh::Selector s_all_nodes);

  void compute_favre_stress(
    const std::string &averageBlockName,
    const double &oldTimeFilter,
    const double &zeroCurrent,
    const double &dt,
    stk::mesh::Selector s_all_nodes);

  // hold the realm
  Realm &realm_;
  
  double currentTimeFilter_; /* provided by restart */
  double timeFilterInterval_; /* user supplied */
  bool forcedReset_; /* allows forhard reset */

  // vector of averaging information
  std::vector<AveragingInfo *> averageInfoVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
