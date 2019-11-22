/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef DataProbePostProcessing_h
#define DataProbePostProcessing_h

#include <NaluParsing.h>

#include <string>
#include <vector>
#include <utility>

// stk_mesh/base/fem
#include <stk_mesh/base/Selector.hpp>

// stk forwards
namespace stk {
  namespace mesh {
    class BulkData;
    class FieldBase;
    class MetaData;
    class Part;
    //class Selector; ? why is this?
    struct Entity;
    typedef std::vector< Part * > PartVector;
  }
}

namespace sierra{
namespace nalu{

class Realm;
class Transfer;
class Transfers;

class DataProbeInfo {
public:
   DataProbeInfo() : numProbes_(0) {}
  ~DataProbeInfo() {}
  
  // for each type of probe, e.g., line of site, hold some stuff
  int numProbes_;
  // manage the types of probes
  std::vector<int> isLineOfSite_;
  std::vector<int> isRing_;
  std::vector<std::string> partName_;
  std::vector<int> probeOnThisRank_;
  std::vector<int> numPoints_;
  std::vector<int> numLinePoints_;
  std::vector<int> numTotalPoints_;
  std::vector<int> generateNewIds_;

  // line of nodes
  std::vector<Coordinates> tipCoordinates_;
  std::vector<Coordinates> tailCoordinates_;

  // ring unit normal and origin (vector or vectors to allow for 3D and 2D rotations)
  std::vector<std::array<double, 3> > unitNormal_;
  std::vector<std::array<double, 3> > originCoordinates_;
  std::vector<std::array<double, 3> > seedCoordinates_;
  
  std::vector<std::vector<stk::mesh::Entity> > nodeVector_;
  std::vector<stk::mesh::Part *> part_;
};

class ProbeType {
public:
 ProbeType() : dataProbeInfo_(new DataProbeInfo()) {}
  ~ProbeType() {delete dataProbeInfo_;} 
  DataProbeInfo *dataProbeInfo_;

  void setup();
  void increment();
};

class LineOfSiteProbeType : public ProbeType {
public:
 LineOfSiteProbeType() {}
 ~LineOfSiteProbeType() {}
};

class RingProbeType : public ProbeType {
public:
 RingProbeType() {}
 ~RingProbeType() {}
};

class DataProbeSpecInfo {
public:
  DataProbeSpecInfo();
  ~DataProbeSpecInfo();

  std::string xferName_;
  std::vector<std::string> fromTargetNames_;
  
  // vector of averaging information
  std::vector<DataProbeInfo *> dataProbeInfo_;
 
  // homegeneous collection of fields over each specification
  std::vector<std::pair<std::string, std::string> > fromToName_;
  std::vector<std::pair<std::string, int> > fieldInfo_;

  // vector of probe types
  std::vector<ProbeType *> probeTypeVec_;
};

class DataProbePostProcessing
{
public:
  
  DataProbePostProcessing(
    Realm &realm,
    const YAML::Node &node);
  ~DataProbePostProcessing();
  
  // load all of the options
  void load(
    const YAML::Node & node);

  // setup part creation and nodal field registration (before populate_mesh())
  void setup();

  // setup part creation and nodal field registration (after populate_mesh())
  void initialize();

  void register_field(
    const std::string fieldName,
    const int fieldSize,
    stk::mesh::MetaData &metaData,
    stk::mesh::Part *part);

  void review( 
    const DataProbeInfo *probeInfo);

  // create the transfer and hold the vector in the DataProbePostProcessing class
  void create_transfer();

  // populate nodal field and output norms (if appropriate)
  void execute();

  // output to a file
  void provide_output(const double currentTime);
  
  // general rotation matrix about a unit normal centered at origin (0,0,0)
  void compute_R(const double theta, const std::vector<double> &u, std::vector<double> &R);
  
  // provide a 3x3 * 3x1 multiply
  void mat_vec(const std::vector<double> &coord, const std::vector<double> &R, std::vector<double> &newCoord); 

  // hold the realm
  Realm &realm_;

  // frequency of output
  int outputFreq_;

  // width for output
  int w_;

  // xfer specifications
  std::string searchMethodName_;
  double searchTolerance_;
  double searchExpansionFactor_;

  // vector of specifications
  std::vector<DataProbeSpecInfo *> dataProbeSpecInfo_;
 
  // hold the transfers
  Transfers *transfers_;
};

} // namespace nalu
} // namespace Sierra

#endif
