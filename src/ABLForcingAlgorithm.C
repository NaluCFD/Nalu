
#include "ABLForcingAlgorithm.h"
#include "Realm.h"
#include "xfer/Transfer.h"
#include "xfer/Transfers.h"
#include "utils/LinearInterpolation.h"

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace sierra {
namespace nalu {

ABLForcingAlgorithm::ABLForcingAlgorithm(Realm& realm, const YAML::Node& node)
  : realm_(realm),
    momSrcType_(ABLForcingAlgorithm::OFF),
    tempSrcType_(ABLForcingAlgorithm::OFF),
    alphaMomentum_(1.0),
    alphaTemperature_(1.0),
    velHeights_(0),
    tempHeights_(0),
    velXTimes_(0),
    velYTimes_(0),
    velZTimes_(0),
    tempTimes_(0),
    velX_(0),
    velY_(0),
    velZ_(0),
    temp_(0),
    UmeanCalc_(0),
    rhoMeanCalc_(0),
    USource_(0),
    TmeanCalc_(0),
    TSource_(0),
    searchMethod_("stk_kdtree"),
    searchTolerance_(1.0e-4),
    searchExpansionFactor_(1.5),
    fromTargetNames_(0),
    velPartNames_(),
    tempPartNames_(),
    allPartNames_(),
    allParts_(0),
    inactiveSelector_(),
    transfers_(0),
    velGenPartList_(false),
    tempGenPartList_(false),
    velPartFmt_(""),
    tempPartFmt_(""),
    outputFreq_(10),
    outFileFmt_("abl_%s_sources.dat")
{
  load(node);
}

ABLForcingAlgorithm::~ABLForcingAlgorithm()
{
  if (NULL != transfers_)
    delete transfers_;
}

void
ABLForcingAlgorithm::load(const YAML::Node& node)
{
  get_if_present(node, "search_method", searchMethod_, searchMethod_);
  get_if_present(node, "search_tolerance", searchTolerance_, searchTolerance_);
  get_if_present(node, "search_expansion_factor", searchExpansionFactor_,
                 searchExpansionFactor_);
  get_if_present(node, "output_frequency", outputFreq_, outputFreq_);
  get_if_present(node, "output_format", outFileFmt_, outFileFmt_);

  if (node["momentum"])
    load_momentum_info(node["momentum"]);

  if (node["temperature"])
    load_temperature_info(node["temperature"]);

  if (node["from_target_part"]) {
    const YAML::Node& fromParts = node["from_target_part"];
    if (fromParts.IsSequence()) {
      fromTargetNames_ = fromParts.as<std::vector<std::string>>();
    } else if (fromParts.IsScalar()) {
      fromTargetNames_.resize(1);
      fromTargetNames_[0] = fromParts.as<std::string>();
    }
  } else if (
    (momSrcType_ != ABLForcingAlgorithm::OFF) ||
    (tempSrcType_ != ABLForcingAlgorithm::OFF)) {
    throw std::runtime_error(
      "No from_target_part specified for ABL forcing function");
  }
}

void
ABLForcingAlgorithm::load_momentum_info(const YAML::Node& node)
{
  std::string mom_type = node["type"].as<std::string>();
  if (mom_type == "user_defined") {
    momSrcType_ = ABLForcingAlgorithm::USER_DEFINED;
  } else if (mom_type == "computed") {
    momSrcType_ = ABLForcingAlgorithm::COMPUTED;
  } else {
    throw std::runtime_error(
      "ABLForcingAlgorithm: Invalid type specification for momentum. "
      "Valid types are: [user_defined, computed]");
  }
  get_if_present(node, "relaxation_factor", alphaMomentum_, alphaMomentum_);
  get_required<std::vector<double>>(node, "heights", velHeights_);
  auto nHeights = velHeights_.size();

  if (node["target_parts"]) { // Explicit target parts list provided
    get_required<std::vector<std::string>>(node, "target_parts", velPartNames_);
    ThrowAssertMsg(
      (nHeights == velPartNames_.size()),
      "ABL Forcing: Mismatch between sizes of heights and target parts "
      "provided for momentum source");
    velGenPartList_ = false;
  } else if (node["target_part_format"]) { // Generate part names from printf
                                           // style string
    get_required<std::string>(node, "target_part_format", velPartFmt_);
    velGenPartList_ = true;
    velPartNames_.resize(nHeights);
  } else {
    throw std::runtime_error(
      "ABL Forcing: No target part(s) provided for momentum forcing.");
  }

  // Load momentum source time histories in temporary data structures, test
  // consistency of data with heights and then recast them into 2-D lookup
  // arrays.
  Array2D<double> vxtmp, vytmp, vztmp;
  get_required<Array2D<double>>(node, "velocity_x", vxtmp);
  get_required<Array2D<double>>(node, "velocity_y", vytmp);
  get_required<Array2D<double>>(node, "velocity_z", vztmp);

  create_interp_arrays(nHeights, vxtmp, velXTimes_, velX_);
  create_interp_arrays(nHeights, vytmp, velYTimes_, velY_);
  create_interp_arrays(nHeights, vztmp, velZTimes_, velZ_);

  const int ndim = realm_.spatialDimension_;
  if (momSrcType_ == COMPUTED) {
    UmeanCalc_.resize(nHeights);
    rhoMeanCalc_.resize(nHeights);
  }

  for (size_t i = 0; i < nHeights; i++) {
    if (momSrcType_ == COMPUTED)
      UmeanCalc_[i].resize(ndim);
  }

  USource_.resize(ndim);
  for (int i = 0; i < ndim; i++) {
    USource_[i].resize(nHeights);
  }
}

void
ABLForcingAlgorithm::load_temperature_info(const YAML::Node& node)
{
  std::string temp_type = node["type"].as<std::string>();
  if (temp_type == "user_defined") {
    tempSrcType_ = ABLForcingAlgorithm::USER_DEFINED;
  } else if (temp_type == "computed") {
    tempSrcType_ = ABLForcingAlgorithm::COMPUTED;
  } else {
    throw std::runtime_error(
      "ABLForcingAlgorithm: Invalid type specification for temperature. "
      "Valid types are: [user_defined, computed]");
  }
  get_if_present(node, "relaxation_factor", alphaTemperature_, alphaTemperature_);
  get_required<std::vector<double>>(node, "heights", tempHeights_);
  auto nHeights = tempHeights_.size();

  if (node["target_parts"]) { // Explicit target parts list provided
    get_required<std::vector<std::string>>(
      node, "target_parts", tempPartNames_);
    ThrowAssertMsg(
      (nHeights == tempPartNames_.size()),
      "ABL Forcing: Mismatch between sizes of heights and target parts "
      "provided for momentum source");
    tempGenPartList_ = false;
  } else if (node["target_part_format"]) { // Generate part names from printf
                                           // style string
    get_required<std::string>(node, "target_part_format", tempPartFmt_);
    tempGenPartList_ = true;
    tempPartNames_.resize(nHeights);
  } else {
    throw std::runtime_error(
      "ABL Forcing: No target part(s) provided for momentum forcing.");
  }

  // Load temperature source time histories, check data consistency and create
  // interpolation lookup tables.
  Array2D<double> temp;
  get_required<Array2D<double>>(node, "temperature", temp);

  create_interp_arrays(nHeights, temp, tempTimes_, temp_);

  TSource_.resize(nHeights);
  if (tempSrcType_ == COMPUTED)
    TmeanCalc_.resize(nHeights);
}

void
ABLForcingAlgorithm::create_interp_arrays(
  const std::vector<double>::size_type nHeights,
  const Array2D<double>& inpArr,
  std::vector<double>& outTimes,
  Array2D<double>& outValues)
{
  /* The input vector is shaped [nTimes, nHeights+1]. We transform it to two
   * arrays:
   *    time[nTimes] = inp[nTimes,0], and
   *    value[nHeights, nTimes] -> swap rows/cols from input
   */

  // Check that all timesteps contain values for all the heights
  for (auto vx : inpArr) {
    ThrowAssert((vx.size() == (nHeights + 1)));
  }

  auto nTimes = inpArr.size();
  outTimes.resize(nTimes);
  outValues.resize(nHeights);
  for (std::vector<double>::size_type i = 0; i < nHeights; i++) {
    outValues[i].resize(nTimes);
  }
  for (std::vector<double>::size_type i = 0; i < nTimes; i++) {
    outTimes[i] = inpArr[i][0];
    for (std::vector<double>::size_type j = 0; j < nHeights; j++) {
      outValues[j][i] = inpArr[i][j + 1];
    }
  }
}

void
ABLForcingAlgorithm::setup()
{
  // Momentum sources
  determine_part_names(
    velHeights_, velPartNames_, velGenPartList_, velPartFmt_);
  // Temperature sources
  determine_part_names(
    tempHeights_, tempPartNames_, tempGenPartList_, tempPartFmt_);
  // Register fields
  register_fields();
}

void
ABLForcingAlgorithm::determine_part_names(
  std::vector<double>& heights,
  std::vector<std::string>& nameSet,
  bool flag,
  std::string& nameFmt)
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  if (flag)
    for (size_t i = 0; i < heights.size(); i++)
      nameSet[i] = (boost::format(nameFmt) % heights[i]).str();

  for (auto partName : nameSet) {
    auto it = allPartNames_.find(partName);
    if (it != allPartNames_.end())
      continue;

    stk::mesh::Part* part = meta.get_part(partName);
    if (NULL == part) {
      // TODO: Need "nodeset" creation capability. Integrate with
      // DataProbePostProcessing to minimize code duplication.
      throw std::runtime_error(
        "ABLForcingAlgorithm::setup: Cannot find part " + partName);
    } else {
      // Add this part to the visited part list
      allPartNames_.insert(partName);
    }
  }
}

void
ABLForcingAlgorithm::register_fields()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  int nDim = meta.spatial_dimension();
  int nStates = realm_.number_of_states();

  for (auto key : allPartNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    VectorFieldType& coords = meta.declare_field<VectorFieldType>(
      stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field(coords, *part, nDim);
  }

  for (auto key : velPartNames_) {

    // Part key
    stk::mesh::Part* part = meta.get_part(key);

    // Velocity
    VectorFieldType& vel = meta.declare_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity", nStates);
    stk::mesh::put_field(vel, *part, nDim);

    // Density
    ScalarFieldType& rho = meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "density", nStates);
    stk::mesh::put_field(rho, *part);
  }

  for (auto key : tempPartNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    ScalarFieldType& temp = meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "temperature");
    stk::mesh::put_field(temp, *part);
  }
}

void
ABLForcingAlgorithm::initialize()
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  // We expect all parts to exist, so no creation step here

  // Add all parts to inactive selection
  for (auto key : allPartNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    allParts_.push_back(part);
  }
  inactiveSelector_ = stk::mesh::selectUnion(allParts_);
  create_transfers();

  if (momSrcType_ != OFF) {
    NaluEnv::self().naluOutputP0()
      << "ABL Forcing active for Momentum Equations\n"
      << "\t Number of planes: " << velHeights_.size()
      << "\n\t Number of time steps: " << velXTimes_.size() << std::endl;
  }
  if (tempSrcType_ != OFF) {
    NaluEnv::self().naluOutputP0()
      << "ABL Forcing active for Temperature Equation\n"
      << "\t Number of planes: " << tempHeights_.size()
      << "\n\t Number of time steps: " << tempTimes_.size() << std::endl
      << std::endl;
  }

  // Prepare output files to dump sources when computed during precursor phase
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( momSrcType_ == COMPUTED )) {
    std::string uxname((boost::format(outFileFmt_)%"Ux").str());
    std::string uyname((boost::format(outFileFmt_)%"Uy").str());
    std::string uzname((boost::format(outFileFmt_)%"Uz").str());
    std::fstream uxFile, uyFile, uzFile;
    uxFile.open(uxname.c_str(), std::fstream::out);
    uyFile.open(uyname.c_str(), std::fstream::out);
    uzFile.open(uzname.c_str(), std::fstream::out);

    uxFile << "# Time, " ;
    uyFile << "# Time, " ;
    uzFile << "# Time, " ;
    for (size_t ih = 0; ih < velHeights_.size(); ih++) {
      uxFile << velHeights_[ih] << ", ";
      uyFile << velHeights_[ih] << ", ";
      uzFile << velHeights_[ih] << ", ";
    }
    uxFile << std::endl ;
    uyFile << std::endl ;
    uzFile << std::endl ;
    uxFile.close();
    uyFile.close();
    uzFile.close();
  }
}

void
ABLForcingAlgorithm::create_transfers()
{
  transfers_ = new Transfers(*realm_.root());

  if (momentumForcingOn()) {
    populate_transfer_data("velocity", velPartNames_);
    populate_transfer_data("density", velPartNames_);
  }

  if (temperatureForcingOn())
    populate_transfer_data("temperature", tempPartNames_);

  transfers_->initialize();
}

void
ABLForcingAlgorithm::populate_transfer_data(
  std::string fieldName, const std::vector<std::string>& partNames)
{
  stk::mesh::MetaData& meta = realm_.meta_data();

  Transfer* theXfer = new Transfer(*transfers_);
  theXfer->name_ = "ablforcing_xfer_" + fieldName;
  theXfer->fromRealm_ = &realm_;
  theXfer->toRealm_ = &realm_;
  theXfer->searchMethodName_ = searchMethod_;
  theXfer->searchTolerance_ = searchTolerance_;
  theXfer->searchExpansionFactor_ = searchExpansionFactor_;

  for (auto key : fromTargetNames_) {
    stk::mesh::Part* part = meta.get_part(key);
    theXfer->fromPartVec_.push_back(part);
  }
  for (auto key : partNames) {
    stk::mesh::Part* part = meta.get_part(key);
    theXfer->toPartVec_.push_back(part);
  }
  theXfer->transferVariablesPairName_.push_back(
    std::make_pair(fieldName, fieldName));
  transfers_->transferVector_.push_back(theXfer);
}

void
ABLForcingAlgorithm::execute()
{
  // Map fields from fluidRealm onto averaging planes
  transfers_->execute();

  if (momentumForcingOn())
    compute_momentum_sources();

  if (temperatureForcingOn())
    compute_temperature_sources();
}

void
ABLForcingAlgorithm::compute_momentum_sources()
{
  const double dt = realm_.get_time_step();
  const double currTime = realm_.get_current_time();
 
  if (momSrcType_ == COMPUTED) {
    calc_mean_velocity();
    for (size_t ih = 0; ih < velHeights_.size(); ih++) {
      double xval, yval;
      
      // Interpolate the velocities from the table to the current time
      utils::linear_interp(velXTimes_, velX_[ih], currTime, xval);
      utils::linear_interp(velYTimes_, velY_[ih], currTime, yval);
      
      // Compute the momentum source
      // Momentum source in the x direction
      USource_[0][ih] = rhoMeanCalc_[ih] * (alphaMomentum_ / dt) * 
                          (xval - UmeanCalc_[ih][0]);
      // Momentum source in the y direction
      USource_[1][ih] = rhoMeanCalc_[ih] * (alphaMomentum_ / dt) * 
                          (yval - UmeanCalc_[ih][1]);

      // No momentum source in z-direction
      USource_[2][ih] = 0.0;

    }
    
  } else {
    for (size_t ih = 0; ih < velHeights_.size(); ih++) {
      utils::linear_interp(velXTimes_, velX_[ih], currTime, USource_[0][ih]);
      utils::linear_interp(velYTimes_, velY_[ih], currTime, USource_[1][ih]);
      utils::linear_interp(velZTimes_, velZ_[ih], currTime, USource_[2][ih]);
    }
  }

  const int tcount = realm_.get_time_step_count();
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( momSrcType_ == COMPUTED ) &&
      ( tcount % outputFreq_ == 0)) {
    std::string uxname((boost::format(outFileFmt_)%"Ux").str());
    std::string uyname((boost::format(outFileFmt_)%"Uy").str());
    std::string uzname((boost::format(outFileFmt_)%"Uz").str());
    std::fstream uxFile, uyFile, uzFile;
    uxFile.open(uxname.c_str(), std::fstream::app);
    uyFile.open(uyname.c_str(), std::fstream::app);
    uzFile.open(uzname.c_str(), std::fstream::app);

    uxFile << std::setw(12) << currTime << ", ";
    uyFile << std::setw(12) << currTime << ", ";
    uzFile << std::setw(12) << currTime << ", ";
    for (size_t ih = 0; ih < velHeights_.size(); ih++) {
      uxFile << std::setprecision(6)
             << std::setw(15)
             << USource_[0][ih] << ", ";
      uyFile << std::setprecision(6)
             << std::setw(15)
             << USource_[1][ih] << ", ";
      uzFile << std::setprecision(6)
             << std::setw(15)
             << USource_[2][ih] << ", ";
    }
    uxFile << std::endl;
    uyFile << std::endl;
    uzFile << std::endl;
    uxFile.close();
    uyFile.close();
    uzFile.close();
  }
}

void
ABLForcingAlgorithm::compute_temperature_sources()
{
  const double dt = realm_.get_time_step();
  const double currTime = realm_.get_current_time();

  if (tempSrcType_ == COMPUTED) {
    calc_mean_temperature();
    for (size_t ih = 0; ih < tempHeights_.size(); ih++) {
      double tval;
      utils::linear_interp(tempTimes_, temp_[ih], currTime, tval);
      TSource_[ih] = (alphaTemperature_ / dt) * (tval - TmeanCalc_[ih]);
    }
  } else {
    for (size_t ih = 0; ih < tempHeights_.size(); ih++) {
      utils::linear_interp(tempTimes_, temp_[ih], currTime, TSource_[ih]);
    }
  }

  const int tcount = realm_.get_time_step_count();
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( tempSrcType_ == COMPUTED ) &&
      ( tcount % outputFreq_ == 0)) {
    std::string fname((boost::format(outFileFmt_)%"T").str());
    std::fstream tFile;
    tFile.open(fname.c_str(), std::fstream::app);

    tFile << currTime << ", ";
    for (size_t ih = 0; ih < tempHeights_.size(); ih++) {
      tFile << std::setprecision(6)
            << std::setw(15)
            << TSource_[ih] << ", ";
    }
    tFile << std::endl;
    tFile.close();
  }
}

void
ABLForcingAlgorithm::calc_mean_velocity()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  stk::mesh::BulkData& bulk = realm_.bulk_data();
  const int nDim = meta.spatial_dimension();

  VectorFieldType* velocity =
    meta.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  ScalarFieldType* density =
    meta.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");

  const size_t numPlanes = velHeights_.size();
  // Sum (velocity and density) and number of nodes on this processor 
  // over all planes
  std::vector<double> sumVel(numPlanes * nDim, 0.0);
  std::vector<double> sumRho(numPlanes, 0.0);
  std::vector<unsigned> numNodes(numPlanes, 0);

  // Global sum and nodes for computing global average
  std::vector<double> sumVelGlobal(numPlanes * nDim, 0.0);
  std::vector<double> sumRhoGlobal(numPlanes, 0.0);
  std::vector<unsigned> totalNodes(numPlanes, 0);

  // Loop for all planes
  for (size_t ih = 0; ih < numPlanes; ih++) {
    const int ioff = ih * nDim;
    stk::mesh::Part* part = meta.get_part(velPartNames_[ih]);
    stk::mesh::Selector s_local_part(*part);
    const stk::mesh::BucketVector& node_buckets =
      bulk.get_buckets(stk::topology::NODE_RANK, s_local_part);

    // Calculate sum (velocity and density) for all nodes on this processor
    for (size_t ib = 0; ib < node_buckets.size(); ib++) {
      stk::mesh::Bucket& bukt = *node_buckets[ib];
      double* velField = stk::mesh::field_data(*velocity, bukt);
      double* rhoField = stk::mesh::field_data(*density, bukt);

      for (size_t in = 0; in < bukt.size(); in++) {

        // Compute the offset dimension to access vleocity components
        const int offset = in * nDim;

        // Velocity sum
        for (int i = 0; i < nDim; i++)
          sumVel[ioff + i] += velField[offset + i];

        // Compute the density sum
        sumRho[ih] += rhoField[in];

      }
      numNodes[ih] += bukt.size();
    }
  }

  // Assemble global sum and node count
  // Velocity
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), sumVel.data(), sumVelGlobal.data(),
    numPlanes * nDim);

  // Density
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), sumRho.data(), sumRhoGlobal.data(),
    numPlanes);

  // Revisit this for area or volume weighted averaging.
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), numNodes.data(), totalNodes.data(),
    numPlanes);

  // Compute spatial averages
  for (size_t ih = 0; ih < numPlanes; ih++) {
    const size_t ioff = ih * nDim;
    
    // Mean Velocity divide by number of total nodes
    for (int i = 0; i < nDim; i++) {
      UmeanCalc_[ih][i] = sumVelGlobal[ioff + i] / totalNodes[ih];
    }
    
      // Mean density devide by number of total nodes
      rhoMeanCalc_[ih] = sumRhoGlobal[ih] / totalNodes[ih];

  }
}

void
ABLForcingAlgorithm::calc_mean_temperature()
{
  stk::mesh::MetaData& meta = realm_.meta_data();
  stk::mesh::BulkData& bulk = realm_.bulk_data();
 
  ScalarFieldType* temperature =
    meta.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  const size_t numPlanes = tempHeights_.size();
  std::vector<double> sumTemp(numPlanes, 0.0), sumTempGlobal(numPlanes, 0.0);
  std::vector<unsigned> numNodes(numPlanes, 0), totalNodes(numPlanes, 0);

  for (size_t ih = 0; ih < numPlanes; ih++) {
    stk::mesh::Part* part = meta.get_part(tempPartNames_[ih]);
    stk::mesh::Selector s_local_part(*part);
    const stk::mesh::BucketVector& node_buckets =
      bulk.get_buckets(stk::topology::NODE_RANK, s_local_part);

    for (size_t ib = 0; ib < node_buckets.size(); ib++) {
      stk::mesh::Bucket& bukt = *node_buckets[ib];
      double* tempField = stk::mesh::field_data(*temperature, bukt);

      for (size_t in = 0; in < bukt.size(); in++) {
        sumTemp[ih] += tempField[in];
      }
      numNodes[ih] += bukt.size();
    }
  }

  // Determine global sum and node count
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), sumTemp.data(), sumTempGlobal.data(),
    numPlanes);
  stk::all_reduce_sum(
    NaluEnv::self().parallel_comm(), numNodes.data(), totalNodes.data(),
    numPlanes);

  // Compute spatial average
  for (size_t ih = 0; ih < numPlanes; ih++) {
    TmeanCalc_[ih] = sumTempGlobal[ih] / totalNodes[ih];
  }
}

void
ABLForcingAlgorithm::eval_momentum_source(
  const double zp, std::vector<double>& momSrc)
{
  const int nDim = realm_.spatialDimension_;
  if (velHeights_.size() == 1) {
    // Constant source term throughout the domain
    for (int i = 0; i < nDim; i++) {
      momSrc[i] = USource_[i][0];
    }
  } else {
    // Linearly interpolate source term within the planes, maintain constant
    // source term above and below the heights provided
    for (int i = 0; i < nDim; i++) {
      utils::linear_interp(velHeights_, USource_[i], zp, momSrc[i]);
    }
  }
}

void
ABLForcingAlgorithm::eval_temperature_source(const double zp, double& tempSrc)
{
  if (tempHeights_.size() == 1) {
    tempSrc = TSource_[0];
  } else {
    utils::linear_interp(tempHeights_, TSource_, zp, tempSrc);
  }
}

} // namespace nalu
} // namespace sierra
