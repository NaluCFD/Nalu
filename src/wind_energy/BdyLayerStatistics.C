/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "wind_energy/BdyLayerStatistics.h"
#include "wind_energy/BdyHeightAlgorithm.h"
#include "Realm.h"
#include "TurbulenceAveragingPostProcessing.h"
#include "AveragingInfo.h"
#include "NaluEnv.h"
#include "utils/LinearInterpolation.h"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"

#include <cmath>
#include <fstream>
#include <string>

namespace sierra {
namespace nalu {

BdyLayerStatistics::BdyLayerStatistics(
  Realm& realm,
  const YAML::Node& node
) : realm_(realm),
    nDim_(realm.spatialDimension_)
{
  load(node);
}

BdyLayerStatistics::~BdyLayerStatistics()
{}

void
BdyLayerStatistics::load(const YAML::Node& node)
{
  const auto& partNames = node["target_name"];
  if (partNames.Type() == YAML::NodeType::Scalar) {
    auto pName = partNames.as<std::string>();
    partNames_.push_back(pName);
  } else {
    partNames_ = partNames.as<std::vector<std::string>>();
  }

  // Process algorithm used to determine unique heights
  std::string heightAlg = "rectilinear_mesh";
  get_if_present(node, "height_calc_algorithm", heightAlg, heightAlg);

  if (heightAlg == "rectilinear_mesh") {
    bdyHeightAlg_.reset(new RectilinearMeshHeightAlg(realm_, node));
  } else {
    throw std::runtime_error("BdyLayerStatistics::load(): Incorrect height algorithm.");
  }

  double timeAvgWindow = 3600.0;
  get_if_present(node, "time_filter_interval", timeAvgWindow, timeAvgWindow);
  get_if_present(node, "compute_temperature_statistics", calcTemperatureStats_,
                 calcTemperatureStats_);

  setup_turbulence_averaging(timeAvgWindow);

  get_if_present(node, "output_frequency", outputFrequency_, outputFrequency_);
  get_if_present(node, "time_hist_output_frequency",
                 timeHistOutFrequency_, timeHistOutFrequency_);
}

void
BdyLayerStatistics::setup_turbulence_averaging(
  const double timeAvgWindow)
{
  bool hasTurbAvg = false;
  if (realm_.turbulenceAveragingPostProcessing_ == nullptr) {
    realm_.turbulenceAveragingPostProcessing_ = new TurbulenceAveragingPostProcessing(realm_);
  } else {
    hasTurbAvg = true;
  }

  auto* turbAvg = realm_.turbulenceAveragingPostProcessing_;

  if (hasTurbAvg) {
    const double diff = std::fabs(timeAvgWindow - turbAvg->timeFilterInterval_);
    if (diff > 1.0e-3)
      NaluEnv::self().naluOutputP0()
        << "WARNING:: BdyLayerStatistics: timeFilterInterval inconsistent with that requested for TurbulenceAveragingPostProcessing." << std::endl;
  } else {
    turbAvg->timeFilterInterval_ = timeAvgWindow;
  }

  AveragingInfo* avInfo = new AveragingInfo();

  avInfo->name_ = "abl";
  avInfo->targetNames_ = partNames_;
  avInfo->computeSFSStress_ = true;
  avInfo->computeResolvedStress_ = true;
  avInfo->resolvedFieldNameVec_.push_back("velocity");

  if (calcTemperatureStats_) {
    avInfo->computeTemperatureResolved_ = true;
    avInfo->computeTemperatureSFS_ = true;
    avInfo->resolvedFieldNameVec_.push_back("temperature");
  }

  turbAvg->averageInfoVec_.push_back(avInfo);
}

void
BdyLayerStatistics::setup()
{
  auto& meta = realm_.meta_data();
  const size_t nparts = partNames_.size();
  fluidParts_.resize(nparts);

  for (size_t i=0; i < nparts; i++) {
    auto* part = meta.get_part(realm_.physics_part_name(partNames_[i]));
    if (nullptr == part)
      throw std::runtime_error("BdyLayerStatistics:: Part not found: " + partNames_[i]);
    else
      fluidParts_[i] = part;
  }

  heightIndex_ = &meta.declare_field<ScalarIntFieldType>(
    stk::topology::NODE_RANK, "bdy_layer_height_index_field");
  for (auto* part: fluidParts_)
    stk::mesh::put_field(*heightIndex_, *part);
}

void
BdyLayerStatistics::initialize()
{
  auto& meta = realm_.meta_data();
  stk::mesh::Selector sel = meta.locally_owned_part()
    & stk::mesh::selectUnion(fluidParts_);

  bdyHeightAlg_->calc_height_levels(sel, *heightIndex_, heights_);

  const size_t nHeights = heights_.size();
  sumVol_.resize(nHeights);
  velAvg_.resize(nHeights * nDim_);
  velBarAvg_.resize(nHeights * nDim_);
  uiujAvg_.resize(nHeights * nDim_ * 2);
  sfsAvg_.resize(nHeights * nDim_ * 2);


  if (calcTemperatureStats_) {
    thetaAvg_.resize(nHeights);
    thetaBarAvg_.resize(nHeights * nDim_);
  }

  doInit_ = false;
}

void
BdyLayerStatistics::execute()
{
  if (doInit_) initialize();

  compute_velocity_stats();
  output_velocity_averages();

  if (calcTemperatureStats_) {
    compute_temperature_stats();
    output_temperature_averages();
  }

}

void
BdyLayerStatistics::velocity(
  double height,
  double* velVector)
{
  interpolate_variable(
    realm_.meta_data().spatial_dimension(),
    velAvg_, height, velVector);
}

void
BdyLayerStatistics::time_averaged_velocity(
  double height,
  double* velVector)
{
  interpolate_variable(
    realm_.meta_data().spatial_dimension(),
    velBarAvg_, height, velVector);
}

void
BdyLayerStatistics::interpolate_variable(
  int nComp,
  std::vector<double>& varArray,
  double height,
  double* interpVar)
{
  auto idx = utils::find_index(heights_, height);

  switch (idx.first) {
  case utils::OutOfBounds::LOWLIM: {
    int offset = idx.second * nComp;
    for (int d=0; d < nComp; d++) {
      interpVar[d] = varArray[offset + d];
    }
    break;
  }

  case utils::OutOfBounds::UPLIM: {
    int offset = (idx.second - 1) * nComp;
    for (int d=0; d < nComp; d++) {
      interpVar[d] = varArray[offset + d];
    }
    break;
  }

  case utils::OutOfBounds::VALID: {
    int ih = idx.second;
    int offset = idx.second * nComp;
    double fac = (height - heights_[ih]) / (heights_[ih+1] - heights_[ih]);
    for (int d=0; d < nComp; d++) {
      interpVar[d] = (1.0 - fac) * varArray[offset+d] + fac * varArray[offset + nComp + d];
    }
    break;
  }
  }
}

void
BdyLayerStatistics::compute_velocity_stats()
{
  auto& meta = realm_.meta_data();
  auto& bulk = realm_.bulk_data();
  stk::mesh::Selector sel = meta.locally_owned_part()
    & stk::mesh::selectUnion(fluidParts_)
    & !(realm_.get_inactive_selector());

  const auto bkts = bulk.get_buckets(stk::topology::NODE_RANK, sel);

  ScalarFieldType* density = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
  VectorFieldType* velocity = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "velocity");
  VectorFieldType* velTimeAvg = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "velocity_resa_abl");
  stk::mesh::FieldBase* resStress = meta.get_field(
    stk::topology::NODE_RANK, "resolved_stress");
  stk::mesh::FieldBase* sfsField = meta.get_field(
    stk::topology::NODE_RANK, "sfs_stress");
  stk::mesh::FieldBase* dualVol = meta.get_field(
    stk::topology::NODE_RANK, "dual_nodal_volume");

  const size_t nHeights = heights_.size();
  // Local summation arrays
  std::vector<double> vol(nHeights, 0.0), gVol(nHeights, 0.0);
  std::vector<double> velMean(nHeights * nDim_, 0.0);
  std::vector<double> velBarMean(nHeights * nDim_, 0.0);
  std::vector<double> uiujMean(nHeights * nDim_ * 2, 0.0);
  std::vector<double> sfsMean(nHeights * nDim_ * 2, 0.0);
  std::vector<double> rhoMean(nHeights, 0.0);

  // Global reduction arrays
  std::vector<double> gVelMean(nHeights * nDim_, 0.0);
  std::vector<double> gVelBarMean(nHeights * nDim_, 0.0);
  std::vector<double> gUiujMean(nHeights * nDim_ * 2, 0.0);
  std::vector<double> gSfsMean(nHeights * nDim_ * 2, 0.0);
  std::vector<double> gRhoMean(nHeights);

  // Sum up all the local contributions
  for (auto b: bkts) {
    for (size_t in=0; in < b->size(); in++) {
      auto node = (*b)[in];
      int ih = *stk::mesh::field_data(*heightIndex_, node);
      int offset = ih * nDim_;

      // Volume calcs
      double dVol = *(double*)stk::mesh::field_data(*dualVol, node);
      vol[ih] += dVol;

      // Density calcs
      double rho = *stk::mesh::field_data(*density, node);
      rhoMean[ih] += rho * dVol;

      // Velocity calculations
      {
        double* vel = stk::mesh::field_data(*velocity, node);
        double* velA = stk::mesh::field_data(*velTimeAvg, node);

        for (int d=0; d < nDim_; d++) {
          velMean[offset + d] += vel[d] * rho * dVol;
          velBarMean[offset + d] += velA[d] * dVol;
        }
      }

      // Stress calculations
      offset *= 2;
      {
        double* sfs = static_cast<double*>(stk::mesh::field_data(*sfsField, node));
        double* uiuj = static_cast<double*>(stk::mesh::field_data(*resStress, node));
        int ix = 0;

        for (int i=0; i < nDim_; i++)
          for (int j=0; j < nDim_; j++) {
            sfsMean[offset + ix] += sfs[ix] * dVol;
            uiujMean[offset + ix] += uiuj[ix] * dVol;
            ix++;
          }
      }
    }
  }

  // Global summation
  stk::all_reduce_sum(bulk.parallel(), velMean.data(), gVelMean.data(), nHeights*nDim_);
  stk::all_reduce_sum(bulk.parallel(), velBarMean.data(), gVelBarMean.data(), nHeights*nDim_);
  stk::all_reduce_sum(bulk.parallel(), sfsMean.data(), gSfsMean.data(), nHeights*nDim_*2);
  stk::all_reduce_sum(bulk.parallel(), uiujMean.data(), gUiujMean.data(), nHeights*nDim_*2);
  stk::all_reduce_sum(bulk.parallel(), vol.data(), gVol.data(), nHeights);
  stk::all_reduce_sum(bulk.parallel(), rhoMean.data(), gRhoMean.data(), nHeights);

  // Compute averages
  for (size_t ih=0; ih < nHeights; ih++) {
    int offset = ih * nDim_;

    for (int d=0; d < nDim_; d++) {
      velAvg_[offset + d] = gVelMean[offset + d] / gRhoMean[ih];
      velBarAvg_[offset + d] = gVelBarMean[offset + d] / gVol[ih];
    }

    offset *= 2;
    for (int i=0; i < nDim_ * 2; i++) {
      sfsAvg_[offset + i] = gSfsMean[offset + i] / gVol[ih];
      uiujAvg_[offset + i] = gUiujMean[offset + i] / gVol[ih];
    }

    // Store sum volumes for temperature stats (processed next)
    sumVol_[ih] = gVol[ih];
  }
}

void
BdyLayerStatistics::compute_temperature_stats()
{
  auto& meta = realm_.meta_data();
  auto& bulk = realm_.bulk_data();
  stk::mesh::Selector sel = meta.locally_owned_part()
    & stk::mesh::selectUnion(fluidParts_)
    & !(realm_.get_inactive_selector());

  const auto bkts = bulk.get_buckets(stk::topology::NODE_RANK, sel);
  const size_t nHeights = heights_.size();

  ScalarFieldType* theta = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "temperature");
  ScalarFieldType* thetaA = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "temperature_resa_abl");
  ScalarFieldType* dualVol = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "dual_nodal_volume");

  std::vector<double> tMean(nHeights, 0.0);
  std::vector<double> tAMean(nHeights, 0.0);
  std::vector<double> gTMean(nHeights, 0.0);
  std::vector<double> gTAMean(nHeights, 0.0);

  // Sum up all local contributions
  for (auto b: bkts) {
    for (size_t in=0; in < b->size(); in++) {
      auto node = (*b)[in];
      int ih = *stk::mesh::field_data(*heightIndex_, node);

      // Temperature calculations
      double* temp = stk::mesh::field_data(*theta, node);
      double* tempA = stk::mesh::field_data(*thetaA, node);
      double* dVol = stk::mesh::field_data(*dualVol, node);

      tMean[ih] += temp[0] * dVol[0];
      tAMean[ih] += tempA[0] * dVol[0];

      // TODO: Implement temperature sfs and variance statistics
    }
  }

  // Global summation
  stk::all_reduce_sum(bulk.parallel(), tMean.data(), gTMean.data(), nHeights);
  stk::all_reduce_sum(bulk.parallel(), tAMean.data(), gTAMean.data(), nHeights);

  // Compute averages
  for (size_t ih=0; ih < nHeights; ih++) {
    thetaAvg_[ih] = gTMean[ih] / sumVol_[ih];
    thetaBarAvg_[ih] = gTAMean[ih] / sumVol_[ih];
  }
}

void
BdyLayerStatistics::output_velocity_averages()
{
  const int tStep = realm_.get_time_step_count();
  const int iproc = realm_.bulk_data().parallel_rank();

  // Only output data if at the desired timestep
  if ((iproc != 0) || (tStep % outputFrequency_ != 0)) return;

  std::ofstream velfile;
  std::ofstream uiujfile;
  std::ofstream sfsfile;

  // TODO: Allow customizable filenames?
  velfile.open("abl_velocity_stats.dat", std::ofstream::out);
  uiujfile.open("abl_resolved_stress_stats.dat", std::ofstream::out);
  sfsfile.open("abl_sfs_stress_stats.dat", std::ofstream::out);

  std::string curTime = std::to_string(realm_.get_current_time());
  velfile << "# Time = " << curTime << std::endl;
  uiujfile << "# Time = " << curTime << std::endl;
  sfsfile << "# Time = " << curTime << std::endl;
  velfile << "# Height, <Ux>, <Uy>, <Uz>, Ux, Uy, Uz" << std::endl;
  uiujfile << "# Height, u11, u12, u13, u22, u23, u33" << std::endl;
  sfsfile << "# Height, t11, t12, t13, t22, t23, t33" << std::endl;

  // TODO: Fix format/precision options for output
  const size_t nHeights = heights_.size();
  for (size_t ih=0; ih < nHeights; ih++) {
    int offset = ih * nDim_;

    // Velocity output
    velfile << heights_[ih];
    for (int d=0; d < nDim_; d++)
      velfile << " " << velBarAvg_[offset + d];
    for (int d=0; d < nDim_; d++)
      velfile << " " << velAvg_[offset + d];
    velfile << std::endl;

    // Resolved and SFS stress outputs
    offset *= 2;
    sfsfile << heights_[ih];
    uiujfile << heights_[ih];
    for (int i=0; i < nDim_ * 2; i++) {
      sfsfile << " " << sfsAvg_[offset + i];
      uiujfile << " " << uiujAvg_[offset + i];
    }
    sfsfile << std::endl;
    uiujfile << std::endl;
  }

  velfile.close();
  sfsfile.close();
  uiujfile.close();
}

void
BdyLayerStatistics::output_temperature_averages()
{
  const int tStep = realm_.get_time_step_count();
  const int iproc = realm_.bulk_data().parallel_rank();

  // Only output data if at the desired timestep
  if ((iproc != 0) || (tStep % outputFrequency_ != 0)) return;

  std::ofstream tempfile;
  tempfile.open("abl_temperature_stats.dat", std::ofstream::out);

  std::string curTime = std::to_string(realm_.get_current_time());
  tempfile << "# Time = " << curTime << std::endl;
  tempfile << "# Height, <T>, T" << std::endl;

  const size_t nHeights = heights_.size();
  for (size_t ih=0; ih < nHeights; ih++) {
    // temperature outputs
    tempfile << heights_[ih] << " "
             << thetaBarAvg_[ih] << " "
             << thetaAvg_[ih] << std::endl;
  }

  tempfile.close();
}

}  // nalu
}  // sierra
