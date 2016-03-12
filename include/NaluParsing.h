/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NaluParsing_h
#define NaluParsing_h

#include <BoundaryConditions.h>
#include <Enums.h>
#include <InitialConditions.h>
#include <MaterialPropertys.h>
#include <NaluParsingHelper.h>
#include <NaluEnv.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace sierra {
namespace nalu {

// our data types
struct Velocity {
  double ux_, uy_, uz_;
  Velocity()
    : ux_(0.0), uy_(0.0), uz_(0.0)
  {}
};

struct Coordinates {
  double x_, y_, z_;
  Coordinates()
    : x_(0.0), y_(0.0), z_(0.0)
  {}
};

struct Pressure {
  double pressure_;
  Pressure()
    : pressure_(0.0)
  {}
};

struct TurbKinEnergy {
  double turbKinEnergy_;
  TurbKinEnergy()
    : turbKinEnergy_(0.0)
  {}
};

struct SpecDissRate {
  double specDissRate_;
  SpecDissRate()
    : specDissRate_(0.0)
  {}
};

struct Temperature {
  double temperature_;
  Temperature()
    : temperature_(0.0)
  {}
};

struct MixtureFraction {
  double mixFrac_;
  MixtureFraction()
    : mixFrac_(0.0)
  {}
};

struct MassFraction {
  std::vector<double> massFraction_;
  MassFraction()
  {}
};

struct Emissivity {
  double emissivity_;
  Emissivity()
    : emissivity_(1.0)
  {}
};

struct Irradiation {
  double irradiation_;
  Irradiation()
    : irradiation_(1.0)
  {}
};

struct Transmissivity {
  double transmissivity_;
  Transmissivity()
    : transmissivity_(0.0)
  {}
};

struct EnvironmentalT {
  double environmentalT_;
  EnvironmentalT()
    : environmentalT_(298.0)
  {}
};

struct ReferenceTemperature {
  double referenceTemperature_;
  ReferenceTemperature()
    : referenceTemperature_(298.0)
  {}
};

struct HeatTransferCoefficient {
  double heatTransferCoefficient_;
  HeatTransferCoefficient()
    : heatTransferCoefficient_(0.0)
  {}
};

struct RobinCouplingParameter {
  double robinCouplingParameter_;
  RobinCouplingParameter()
    : robinCouplingParameter_(0.0)
  {}
};

// base class
struct UserData 
{
  std::map<std::string, bool> bcDataSpecifiedMap_;
  std::map<std::string, UserDataType> bcDataTypeMap_;
  std::map<std::string, std::string> userFunctionMap_;
  std::map<std::string, std::vector<double> > functionParams_;

  // FIXME: must elevate temperature due to the temperature_bc_setup method
  Temperature temperature_;
  bool tempSpec_; 
  
UserData() : tempSpec_(false) {}
};

struct NormalHeatFlux {
  double qn_;
  NormalHeatFlux()
    : qn_(0.0)
  {}
};

struct MasterSlave {
  std::string master_;
  std::string slave_;
  MasterSlave() {}
};

// packaged
struct WallUserData : public UserData {
  Velocity u_;
  Velocity dx_;
  TurbKinEnergy tke_;
  MixtureFraction mixFrac_;
  MassFraction massFraction_;
  Emissivity emissivity_;
  Irradiation irradiation_;
  Transmissivity transmissivity_;
  EnvironmentalT environmentalT_;
  NormalHeatFlux q_;
  ReferenceTemperature referenceTemperature_;
  HeatTransferCoefficient heatTransferCoefficient_;
  RobinCouplingParameter robinCouplingParameter_;
  Pressure pressure_;
  
  bool isAdiabatic_;
  bool heatFluxSpec_;
  bool isInterface_;
  bool refTempSpec_;
  bool htcSpec_;
  bool robinParameterSpec_;
  bool irradSpec_;
  bool emissSpec_;

  bool wallFunctionApproach_;

  bool isFsiInterface_;

  WallUserData()
    : UserData(),
      isAdiabatic_(false),
      heatFluxSpec_(false),
      isInterface_(false),
      refTempSpec_(false),
      htcSpec_(false),
      robinParameterSpec_(false),
      irradSpec_(false),
      wallFunctionApproach_(false),
      isFsiInterface_(false) {}
};

struct InflowUserData : public UserData {
  Velocity u_;
  TurbKinEnergy tke_;
  SpecDissRate sdr_;
  MixtureFraction mixFrac_;
  MassFraction massFraction_;
 
  bool uSpec_;
  bool tkeSpec_;
  bool sdrSpec_;
  bool mixFracSpec_;
  bool massFractionSpec_;
  
  InflowUserData()
    : UserData(),
      uSpec_(false), tkeSpec_(false), sdrSpec_(false), mixFracSpec_(false), massFractionSpec_(false)
  {}
};

struct OpenUserData : public UserData {
  Velocity u_;
  Pressure p_;
  TurbKinEnergy tke_;
  SpecDissRate sdr_;
  MixtureFraction mixFrac_;
  MassFraction massFraction_;
 
  bool uSpec_;
  bool pSpec_;
  bool tkeSpec_;
  bool sdrSpec_;
  bool mixFracSpec_;
  bool massFractionSpec_;

  OpenUserData()
    : UserData(),
      uSpec_(false), pSpec_(false), tkeSpec_(false), sdrSpec_(false), mixFracSpec_(false), massFractionSpec_(false)
  {}
};

struct ContactUserData : public UserData {

  double maxSearchRadius_;
  double minSearchRadius_;
  std::vector<std::string> searchBlock_;
  double extrusionDistance_;
  bool useExtrusionAlg_;
  std::string searchMethodName_;
  double expandBoxPercentage_;
  bool clipIsoParametricCoords_;
  bool useHermiteInterpolation_;

  ContactUserData()
    : UserData(),
      maxSearchRadius_(0.0), minSearchRadius_(0.0),
      extrusionDistance_(0.0), useExtrusionAlg_(false), searchMethodName_("na"), expandBoxPercentage_(0.0),
      clipIsoParametricCoords_(false), useHermiteInterpolation_(false)
  {}
};

struct OversetUserData : public UserData {
  // at present, simulation can have one background mesh with multiple, non-interacting overset blocks
  double percentOverlap_;
  bool clipIsoParametricCoords_;
  bool detailedOutput_;
  std::string backgroundBlock_;
  std::string backgroundSurface_;
  std::string backgroundCutBlock_;
  std::string oversetSurface_;
  std::vector<std::string> oversetBlockVec_;
 OversetUserData()
   : UserData(),
    percentOverlap_(10.0), clipIsoParametricCoords_(false), detailedOutput_(false), backgroundBlock_("na"),
    backgroundSurface_("na"), backgroundCutBlock_("na"), oversetSurface_("na")
    {} 
};
 
struct SymmetryUserData : public UserData {
  SymmetryUserData()
    : UserData()
  {/* nothing yet*/}
};

struct PeriodicUserData : public UserData {

  double searchTolerance_;
  std::string searchMethodName_;

  PeriodicUserData()
    : UserData(),
      searchTolerance_(1.0e-16),
      searchMethodName_("na")
  {}
};

struct NonConformalUserData : public UserData {

  std::string searchMethodName_;
  double expandBoxPercentage_;
  bool clipIsoParametricCoords_;
  double searchTolerance_;

  NonConformalUserData()
    : UserData(),
      searchMethodName_("na"), expandBoxPercentage_(0.0),
    clipIsoParametricCoords_(false), searchTolerance_(1.0e-16)
  {}
};

struct WallBoundaryConditionData : public BoundaryCondition {
  WallBoundaryConditionData(BoundaryConditions& bcs) : BoundaryCondition(bcs){};
  WallUserData userData_;
};

struct InflowBoundaryConditionData : public BoundaryCondition {
  InflowBoundaryConditionData(BoundaryConditions& bcs) : BoundaryCondition(bcs){};
  InflowUserData userData_;
};

struct OpenBoundaryConditionData : public BoundaryCondition {
  OpenBoundaryConditionData(BoundaryConditions& bcs) : BoundaryCondition(bcs){};
  OpenUserData userData_;
};

struct OversetBoundaryConditionData : public BoundaryCondition {
  OversetBoundaryConditionData(BoundaryConditions& bcs) : BoundaryCondition(bcs){};
  OversetUserData userData_;
};

struct ContactBoundaryConditionData : public BoundaryCondition {
  ContactBoundaryConditionData(BoundaryConditions& bcs) : BoundaryCondition(bcs){};
  ContactUserData userData_;
};

struct SymmetryBoundaryConditionData : public BoundaryCondition {
  SymmetryBoundaryConditionData(BoundaryConditions& bcs) : BoundaryCondition(bcs){};
  SymmetryUserData userData_;
};

struct PeriodicBoundaryConditionData : public BoundaryCondition {
  PeriodicBoundaryConditionData(BoundaryConditions& bcs) : BoundaryCondition(bcs){};
  MasterSlave masterSlave_;
  PeriodicUserData userData_;
};

struct NonConformalBoundaryConditionData : public BoundaryCondition {
  NonConformalBoundaryConditionData(BoundaryConditions& bcs) : BoundaryCondition(bcs){};
  MasterSlave masterSlave_;
  NonConformalUserData userData_;
};

struct BoundaryConditionOptions{
  std::string bcSetName_;
  WallBoundaryConditionData wallbc_;
  InflowBoundaryConditionData inflowbc_;
  OpenBoundaryConditionData openbc_;
  OversetBoundaryConditionData oversetbc_;
  ContactBoundaryConditionData contactbc_;
  NonConformalBoundaryConditionData nonConformalbc_;
  SymmetryBoundaryConditionData symmetrybc_;
  PeriodicBoundaryConditionData periodicbc_;
};

struct MeshInput {
  std::string meshName_;
};

// initial conditions
struct ConstantInitialConditionData : public InitialCondition {
  ConstantInitialConditionData(InitialConditions& ics) : InitialCondition(ics) {}
  std::vector<std::string> fieldNames_;
  std::vector<std::vector<double> > data_;
};

struct UserFunctionInitialConditionData : public InitialCondition {
  UserFunctionInitialConditionData(InitialConditions& ics) : InitialCondition(ics) {}
  std::map<std::string, std::string> functionNames_;
  std::map<std::string, std::vector<double> > functionParams_;
};

// now the extraction operators for these types
void operator >> (const YAML::Node& node, Velocity& v);
void operator >> (const YAML::Node& node, Coordinates& x);
void operator >> (const YAML::Node& node, Pressure& p);
void operator >> (const YAML::Node& node, TurbKinEnergy& tke);
void operator >> (const YAML::Node& node, SpecDissRate& sdr);
void operator >> (const YAML::Node& node, Temperature& t);
void operator >> (const YAML::Node& node, MixtureFraction& z);
void operator >> (const YAML::Node& node, MassFraction& yk);
void operator >> (const YAML::Node& node, std::map<std::string,double>& mapName);
void operator >> (const YAML::Node& node, std::map<std::string,std::string>& mapName);
void operator >> (const YAML::Node& node, std::map<std::string,std::vector<std::string> >& mapName);

void operator >> (const YAML::Node& node, WallUserData& wallData);
void operator >> (const YAML::Node& node, InflowUserData& inflowData);
void operator >> (const YAML::Node& node, OpenUserData& openData);
void operator >> (const YAML::Node& node, OversetUserData& bcData);
void operator >> (const YAML::Node& node, ContactUserData& bcData);
void operator >> (const YAML::Node& node, NonConformalUserData& bcData);
void operator >> (const YAML::Node& node, SymmetryUserData& bcData);
void operator >> (const YAML::Node& node, PeriodicUserData& bcData);
void operator >> (const YAML::Node& node, BoundaryConditionOptions& bcOptions);
void operator >> (const YAML::Node& node, WallBoundaryConditionData& wallBC);
void operator >> (const YAML::Node& node, InflowBoundaryConditionData& inflowBC);
void operator >> (const YAML::Node& node, OpenBoundaryConditionData& openBC);
void operator >> (const YAML::Node& node, OversetBoundaryConditionData& oversetBC);
void operator >> (const YAML::Node& node, ContactBoundaryConditionData& contactBC);
void operator >> (const YAML::Node& node, NonConformalBoundaryConditionData& nonConformalBC);
void operator >> (const YAML::Node& node, SymmetryBoundaryConditionData& symmetryBC);
void operator >> (const YAML::Node& node, PeriodicBoundaryConditionData& periodicBC);
void operator >> (const YAML::Node& node, MeshInput& meshInput);
void operator >> (const YAML::Node& node, ConstantInitialConditionData& constIC);
void operator >> (const YAML::Node& node, UserFunctionInitialConditionData& fcnIC);

/// Set @param result if the @param key is present in the @param node, else set it to the given default value
template<typename T>
void get_if_present(const YAML::Node & node, const std::string& key, T& result, const T& default_if_not_present = T())
{
  const YAML::Node *value = node.FindValue(key);
  if (value)
    *value >> result;
  else
    result = default_if_not_present;
}

/// this version doesn't change @param result unless the @param key is present in the @param node
template<typename T>
void get_if_present_no_default(const YAML::Node & node, const std::string& key, T& result)
{
  const YAML::Node *value = node.FindValue(key);
  if (value)
    *value >> result;
}

/// this version requires the @param key to be present
template<typename T>
void get_required(const YAML::Node & node, const std::string& key, T& result)
{
  const YAML::Node *value = node.FindValue(key);
  if (value)
    *value >> result;
  else
    {
      if (!NaluEnv::self().parallel_rank()) {
        std::ostringstream err_msg;
        err_msg << "\n\nError: parsing missing required key: " << key 
                  << " at " << NaluParsingHelper::line_info(node)
                  << " for node= " << std::endl;
        NaluParsingHelper::emit(err_msg, node);
        std::cout << err_msg.str() << std::endl;
      }
      throw std::runtime_error("Error: parsing missing required key: " + key);
    }
}

/// these can be used to check and ensure a type of yaml node is as expected
const YAML::Node *
expect_type(const YAML::Node& node, const std::string& key, YAML::NodeType::value type, bool optional=false);

const YAML::Node *
expect_null(const YAML::Node& node, const std::string& key, bool optional=false);

const YAML::Node *
expect_scalar(const YAML::Node& node, const std::string& key, bool optional=false);

const YAML::Node *
expect_sequence(const YAML::Node& node, const std::string& key, bool optional=false);

const YAML::Node *
expect_map(const YAML::Node& node, const std::string& key, bool optional=false);



} // namespace nalu
} // namespace Sierra

#endif
