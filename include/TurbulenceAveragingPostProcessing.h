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

  void register_field(
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

  void review( 
    const AveragingInfo *avInfo);

  // populate nodal field and output norms (if appropriate)
  void execute();

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
