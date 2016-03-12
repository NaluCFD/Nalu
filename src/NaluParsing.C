/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

#include <NaluParsing.h>
#include <NaluEnv.h>
#include <Simulation.h>
#include <Enums.h>

#include <stk_util/environment/ReportHandler.hpp>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <string>

namespace sierra{
namespace nalu{

// now the extraction operators for these types
void operator >> (const YAML::Node& node, Velocity& v) {
  node[0] >> v.ux_;
  node[1] >> v.uy_;
  if ( node.size() > 2 )
    node[2] >> v.uz_;
}

void operator >> (const YAML::Node& node, Coordinates& cx) {
  node[0] >> cx.x_;
  node[1] >> cx.y_;
  if ( node.size() > 2 )
    node[2] >> cx.z_;
}

void operator >> (const YAML::Node& node, Pressure& p) {
  node >> p.pressure_;
}

void operator >> (const YAML::Node& node, TurbKinEnergy& tke) {
  node >> tke.turbKinEnergy_;
}

void operator >> (const YAML::Node& node, SpecDissRate& sdr) {
  node >> sdr.specDissRate_;
}

void operator >> (const YAML::Node& node, Temperature& t) {
  node >> t.temperature_;
}

void operator >> (const YAML::Node& node, MixtureFraction& z) {
  node >> z.mixFrac_;
}

void operator >> (const YAML::Node& node, MassFraction& yk) {
  const size_t ykSize = node.size();
  yk.massFraction_.resize(ykSize);
  for ( size_t k = 0; k < ykSize; ++k ) {
    node[k] >> yk.massFraction_[k];
  }
}

void operator >> (const YAML::Node& node, Emissivity& emiss) {
  node >> emiss.emissivity_;
}

void operator >> (const YAML::Node& node, Irradiation& irrad) {
  node >> irrad.irradiation_;
}

void operator >> (const YAML::Node& node, Transmissivity& tmiss) {
  node >> tmiss.transmissivity_;
}

void operator >> (const YAML::Node& node, EnvironmentalT& et) {
  node >> et.environmentalT_;
}

void operator >> (const YAML::Node& node, NormalHeatFlux& q) {
  node >> q.qn_;
}

void operator >> (const YAML::Node& node, ReferenceTemperature& rt) {
  node >> rt.referenceTemperature_;
}

void operator >> (const YAML::Node& node, HeatTransferCoefficient& htc) {
  node >> htc.heatTransferCoefficient_;
}

void operator >> (const YAML::Node& node, RobinCouplingParameter& alpha) {
  node >> alpha.robinCouplingParameter_;
}

void operator >> (const YAML::Node& node, MasterSlave& ms) {
  node[0] >> ms.master_;
  node[1] >> ms.slave_;
}

void operator >> (const YAML::Node& node, WallUserData& wallData) {

  // constant data; all optional
  if ( node.FindValue("velocity" )  ) {
    node["velocity"] >> wallData.u_;
    wallData.bcDataSpecifiedMap_["velocity"] = true;
    wallData.bcDataTypeMap_["velocity"] = CONSTANT_UD;
  }
  if ( node.FindValue("mesh_displacement" )  ) {
    node["mesh_displacement"] >> wallData.dx_;
    wallData.bcDataSpecifiedMap_["mesh_displacement"] = true;
    wallData.bcDataTypeMap_["mesh_displacement"] = CONSTANT_UD;
  }
  if ( node.FindValue("turbulent_ke") ) {
    node["turbulent_ke"] >> wallData.tke_;
    wallData.bcDataSpecifiedMap_["turbulent_ke"] = true;
    wallData.bcDataTypeMap_["turbulent_ke"] = CONSTANT_UD;
  }
  if ( node.FindValue("temperature") ) {
    node["temperature"] >> wallData.temperature_;
    wallData.bcDataSpecifiedMap_["temperature"] = true;
    wallData.bcDataTypeMap_["temperature"] = CONSTANT_UD;
    wallData.tempSpec_ = true;
  }
  if ( node.FindValue("mixture_fraction") ) {
    node["mixture_fraction"] >> wallData.mixFrac_;
    wallData.bcDataSpecifiedMap_["mixture_fraction"] = true;
    wallData.bcDataTypeMap_["mixture_fraction"] = CONSTANT_UD;
  }
  if ( node.FindValue("mass_fraction") ) {
    node["mass_fraction"] >> wallData.massFraction_;
    wallData.bcDataSpecifiedMap_["mass_fraction"] = true;
    wallData.bcDataTypeMap_["mass_fraction"] = CONSTANT_UD;
  }
  if ( node.FindValue("emissivity") ) {
    node["emissivity"] >> wallData.emissivity_; 
    wallData.emissSpec_ = true;
  }
  if ( node.FindValue("transmissivity") ) {
    node["transmissivity"] >> wallData.transmissivity_; 
  }
  if ( node.FindValue("environmental_temperature") ) {
    node["environmental_temperature"] >> wallData.environmentalT_; 
  }
  if ( node.FindValue("adiabatic") ) {
    node["adiabatic"] >> wallData.isAdiabatic_; 
  }
  if ( node.FindValue("interface") ) {
    node["interface"] >> wallData.isInterface_; 
  }
  if ( node.FindValue("reference_temperature") ) {
    node["reference_temperature"] >> wallData.referenceTemperature_;
    wallData.refTempSpec_ = true;
  }
  if ( node.FindValue("heat_transfer_coefficient") ) {
    node["heat_transfer_coefficient"] >> wallData.heatTransferCoefficient_;
    wallData.htcSpec_ = true;
  }
  if ( node.FindValue("irradiation") ) {
    node["irradiation"] >> wallData.irradiation_;
    wallData.irradSpec_ = true;
  }
  if ( node.FindValue("robin_coupling_parameter") ) {
    node["robin_coupling_parameter"] >> wallData.robinCouplingParameter_;
    wallData.robinParameterSpec_ = true;
  }
  if ( node.FindValue("use_wall_function")) {
    node["use_wall_function"] >> wallData.wallFunctionApproach_;
  }
  if ( node.FindValue("pressure") ) {
    node["pressure"] >> wallData.pressure_;
    wallData.bcDataSpecifiedMap_["pressure"] = true;
    wallData.bcDataTypeMap_["pressure"] = CONSTANT_UD;
  }
  if ( node.FindValue("fsi_interface") ) {
    node["fsi_interface"] >> wallData.isFsiInterface_;
  }

  // not appropriate
  if ( node.FindValue("specific_dissipation_rate")) {
    throw std::runtime_error("specific_dissipation rate at walls is provided by a model, not the user");
  }

  // function data
  const bool optional = true;
  const YAML::Node *userFcnNode = expect_map(node, "user_function_name", optional);
  if (NULL != userFcnNode ) {
    for (YAML::Iterator i = userFcnNode->begin(); i != userFcnNode->end(); ++i) {
      const YAML::Node & key = i.first();
      const YAML::Node & value = i.second();
      std::string stringName;
      key >> stringName;
      std::string data;
      value >> data;
      wallData.bcDataSpecifiedMap_[stringName] = true;
      wallData.bcDataTypeMap_[stringName] = FUNCTION_UD;
      wallData.userFunctionMap_[stringName] = data;
    }

    // extract function name and parameters
    if (expect_map( node, "user_function_parameters", true)) {
      node["user_function_parameters"] >> wallData.functionParams_;
    }
  }

  if ( node.FindValue("heat_flux") ) {
    node["heat_flux"] >> wallData.q_;
    wallData.heatFluxSpec_ = true;
  }

}

void operator >> (const YAML::Node& node, InflowUserData& inflowData) {

  // optional
  if ( node.FindValue("velocity" )  ) {
    node["velocity"] >> inflowData.u_;
    inflowData.uSpec_ = true;
  }
  if ( node.FindValue("turbulent_ke") ){
    node["turbulent_ke"] >> inflowData.tke_;
    inflowData.tkeSpec_ = true;
  }
  if ( node.FindValue("specific_dissipation_rate") ){
    node["specific_dissipation_rate"] >> inflowData.sdr_;
    inflowData.sdrSpec_ = true;
  }
  if ( node.FindValue("mixture_fraction") ){
    node["mixture_fraction"] >> inflowData.mixFrac_;
    inflowData.mixFracSpec_ = true;
  }
  if ( node.FindValue("mass_fraction") ){
    node["mass_fraction"] >> inflowData.massFraction_;
    inflowData.massFractionSpec_ = true;
  }
  if ( node.FindValue("temperature") ) {
    node["temperature"] >> inflowData.temperature_;
    inflowData.tempSpec_ = true;
  }

  const bool optional = true;
  const YAML::Node *userFcnNode = expect_map(node, "user_function_name", optional);
  if (NULL != userFcnNode ) {
    for (YAML::Iterator i = userFcnNode->begin(); i != userFcnNode->end(); ++i) {
      const YAML::Node & key = i.first();
      const YAML::Node & value = i.second();
      std::string stringName;
      key >> stringName;
      std::string data;
      value >> data;
      inflowData.bcDataSpecifiedMap_[stringName] = true;
      inflowData.bcDataTypeMap_[stringName] = FUNCTION_UD;
      inflowData.userFunctionMap_[stringName] = data;
    }
    
    // extract function name and parameters
    if (expect_map( node, "user_function_parameters", true)) {
      node["user_function_parameters"] >> inflowData.functionParams_;
    }
  }
  
}
void operator >> (const YAML::Node& node, OversetUserData& oversetData){
  //nothing is optional
  if ( node.FindValue("percent_overlap") ) {
    node["percent_overlap"] >> oversetData.percentOverlap_;
  }
  else {
    throw std::runtime_error("One MUST specify overset overlap percentage");
  }

  if ( node.FindValue("background_block") ) {
    node["background_block"] >> oversetData.backgroundBlock_;
  }
  else {
    throw std::runtime_error("One MUST specify background block");
  }

  if ( node.FindValue("overset_block") ) {
    const YAML::Node &oversetBlock = *node.FindValue("overset_block");
    if (oversetBlock.Type() == YAML::NodeType::Scalar) {
      oversetData.oversetBlockVec_.resize(1);
      oversetBlock >> oversetData.oversetBlockVec_[0];
    }
    else {
      oversetData.oversetBlockVec_.resize(oversetBlock.size());
      for (size_t i=0; i < oversetBlock.size(); ++i) {
        oversetBlock[i] >> oversetData.oversetBlockVec_[i];
      }
    }
  }
  else {
    throw std::runtime_error("One MUST specify overset block(s)");
  }

  if ( node.FindValue("background_cut_block") ) {
    node["background_cut_block"] >> oversetData.backgroundCutBlock_;
  }
  else {
    throw std::runtime_error("One MUST specify background cut block");
  }

  if ( node.FindValue("background_cut_surface") ) {
    node["background_cut_surface"] >> oversetData.backgroundSurface_;
  }
  else {
    throw std::runtime_error("One MUST specify background cut surface");
  }

  if ( node.FindValue("overset_surface") ) {
    node["overset_surface"] >> oversetData.oversetSurface_;
  }
  else {
    throw std::runtime_error("One MUST specify overset surface");
  }

  if ( node.FindValue("clip_isoparametric_coordinates") ) {
     node["clip_isoparametric_coordinates"] >> oversetData.clipIsoParametricCoords_;
  }

  if ( node.FindValue("detailed_output") ) {
     node["detailed_output"] >> oversetData.detailedOutput_;
  }

}

void operator >> (const YAML::Node& node, ContactUserData& contactData) {
  // nothing is optional
  if ( node.FindValue("max_search_radius" )  ) {
    node["max_search_radius"] >> contactData.maxSearchRadius_;
  }
  else {
    throw std::runtime_error("One MUST specify max search radius at contact bcs");
  }
  
  if ( node.FindValue("min_search_radius" )  ) {
    node["min_search_radius"] >> contactData.minSearchRadius_;
  }
  else {
    throw std::runtime_error("One MUST specify min search radius at contact bcs");
  }

  if ( node.FindValue("search_block" )  ) {
    const YAML::Node &searchBlock = *node.FindValue("search_block");
    if (searchBlock.Type() == YAML::NodeType::Scalar) {
      contactData.searchBlock_.resize(1);
      searchBlock >> contactData.searchBlock_[0];
    }
    else {
      contactData.searchBlock_.resize(searchBlock.size());
      for (size_t i=0; i < searchBlock.size(); ++i) {
        searchBlock[i] >> contactData.searchBlock_[i];
      }
    }
  }
  else {
    throw std::runtime_error("One MUST specify search block at contact bcs");
  }

  if ( node.FindValue("extrusion_distance" )  ) {
    node["extrusion_distance"] >> contactData.extrusionDistance_;
    contactData.useExtrusionAlg_ = true;
  }
  else {
    throw std::runtime_error("Specify extrusion distance at contact bcs; simple halo disabled");
  }
  
  if ( node.FindValue("search_method") ) {
    node["search_method"] >> contactData.searchMethodName_;
  }
  
  if ( node.FindValue("expand_box_percentage" )  ) {
    node["expand_box_percentage"] >> contactData.expandBoxPercentage_;
  }

  if ( node.FindValue("clip_isoparametric_coordinates" )  ) {
     node["clip_isoparametric_coordinates"] >> contactData.clipIsoParametricCoords_;
  }

  if ( node.FindValue("hermite_interpolation" )  ) {
    node["hermite_interpolation"] >> contactData.useHermiteInterpolation_;
 }

}

void operator >> (const YAML::Node& node, OpenUserData& openData) {
  // optional
  if ( node.FindValue("velocity" ) ){
    node["velocity"] >> openData.u_;
    openData.uSpec_ = true;
  }
  // optional
  if ( node.FindValue("turbulent_ke") ) {
    node["turbulent_ke"] >> openData.tke_;
    openData.tkeSpec_ = true;
  }
  if ( node.FindValue("specific_dissipation_rate") ){
    node["specific_dissipation_rate"] >> openData.sdr_;
    openData.sdrSpec_ = true;
  }
  if ( node.FindValue("pressure") ) {
    node["pressure"] >> openData.p_;
    openData.pSpec_ = true;
  }
  if ( node.FindValue("mixture_fraction") ){
    node["mixture_fraction"] >> openData.mixFrac_;
    openData.mixFracSpec_ = true;
  }
  if ( node.FindValue("mass_fraction") ){
    node["mass_fraction"] >> openData.massFraction_;
    openData.massFractionSpec_ = true;
  }
  if ( node.FindValue("temperature") ) {
    node["temperature"] >> openData.temperature_;
    openData.tempSpec_ = true;
  }
}

void operator >> (const YAML::Node& node, SymmetryUserData& symmetryData) {
  // nothing as of yet
}

void operator >> (const YAML::Node& node, PeriodicUserData& periodicData) {
  // nothing is optional
  if ( node.FindValue("search_tolerance" )  ) {
    node["search_tolerance"] >> periodicData.searchTolerance_;
  }
  else {
    throw std::runtime_error("One MUST specify search tolerance at periodic bcs");
  }
  if ( node.FindValue("search_method") ) {
    node["search_method"] >> periodicData.searchMethodName_;
  }
}

void operator >> (const YAML::Node& node, NonConformalUserData& nonConformalData) {

  // everything is optional
  if ( node.FindValue("search_method") ) {
    node["search_method"] >> nonConformalData.searchMethodName_;
  }
  if ( node.FindValue("expand_box_percentage" )  ) {
    node["expand_box_percentage"] >> nonConformalData.expandBoxPercentage_;
  }
  if ( node.FindValue("clip_isoparametric_coordinates" )  ) {
     node["clip_isoparametric_coordinates"] >> nonConformalData.clipIsoParametricCoords_;
  }
  if ( node.FindValue("search_tolerance" )  ) {
    node["search_tolerance"] >> nonConformalData.searchTolerance_;
  }
 
}

void operator >> (const YAML::Node& node, BoundaryConditionOptions& bcOptions) {
  node["boundary_conditions"] >> bcOptions.bcSetName_;
  node["wall_boundary_condition"] >> bcOptions.wallbc_;
  node["inflow_boundary_condition"] >> bcOptions.inflowbc_;
  node["open_boundary_condition"] >> bcOptions.openbc_;
  node["overset_boundary_condition"] >> bcOptions.oversetbc_;
  node["contact_boundary_condition"] >> bcOptions.contactbc_;
  node["symmetry_boundary_condition"] >> bcOptions.symmetrybc_;
  node["periodic_boundary_condition"] >> bcOptions.periodicbc_;
  node["non_confomal_boundary_condition"] >> bcOptions.nonConformalbc_;
}

void operator >> (const YAML::Node& node, WallBoundaryConditionData& wallBC) {
  node["wall_boundary_condition"] >> wallBC.bcName_;
  node["target_name"] >> wallBC.targetName_;
  wallBC.theBcType_ = WALL_BC;
  const YAML::Node& wallUserData = node["wall_user_data"];
  wallUserData >> wallBC.userData_;
  // check for typical rogue line command
  const YAML::Node *value = node.FindValue("user_function_name");
  if ( NULL != value )
    throw std::runtime_error("user_function_data is misplaced; it must be under wall_user_data");
}

void operator >> (const YAML::Node& node, InflowBoundaryConditionData& inflowBC) {
  node["inflow_boundary_condition"] >> inflowBC.bcName_;
  node["target_name"] >> inflowBC.targetName_;
  inflowBC.theBcType_ = INFLOW_BC;
  const YAML::Node& inflowUserData = node["inflow_user_data"];
  inflowUserData >> inflowBC.userData_;
  // check for typical rogue line command
  const YAML::Node *value = node.FindValue("user_function_name");
  if ( NULL != value )
    throw std::runtime_error("user_function_data is misplaced; it must be under inflow_user_data");
}

void operator >> (const YAML::Node& node, OpenBoundaryConditionData& openBC) {
  node["open_boundary_condition"] >> openBC.bcName_;
  node["target_name"] >> openBC.targetName_;
  openBC.theBcType_ = OPEN_BC;
  const YAML::Node& openUserData = node["open_user_data"];
  openUserData >> openBC.userData_;
  // check for typical rogue line command
  const YAML::Node *value = node.FindValue("user_function_name");
  if ( NULL != value )
    throw std::runtime_error("user_function_data is misplaced; it must be under open_user_data");
}

void operator >> (const YAML::Node& node, OversetBoundaryConditionData& oversetBC) {
  node["overset_boundary_condition"] >> oversetBC.bcName_;
  oversetBC.theBcType_ = OVERSET_BC;
  const YAML::Node& oversetUserData = node["overset_user_data"];
  oversetUserData >> oversetBC.userData_;
}

void operator >> (const YAML::Node& node, ContactBoundaryConditionData& contactBC) {
  node["contact_boundary_condition"] >> contactBC.bcName_;
  node["target_name"] >> contactBC.targetName_;
  contactBC.theBcType_ = CONTACT_BC;
  const YAML::Node& contactUserData = node["contact_user_data"];
  contactUserData >> contactBC.userData_;
}

void operator >> (const YAML::Node& node, SymmetryBoundaryConditionData& symmetryBC) {
  node["symmetry_boundary_condition"] >> symmetryBC.bcName_;
  node["target_name"] >> symmetryBC.targetName_;
  symmetryBC.theBcType_ = SYMMETRY_BC;
  const YAML::Node& symmetryUserData = node["symmetry_user_data"];
  symmetryUserData >> symmetryBC.userData_;
}

void operator >> (const YAML::Node& node, PeriodicBoundaryConditionData& periodicBC) {
  node["periodic_boundary_condition"] >> periodicBC.bcName_;
  node["target_name"] >> periodicBC.masterSlave_;
  periodicBC.targetName_ = periodicBC.masterSlave_.master_ + "_" + periodicBC.masterSlave_.slave_;
  periodicBC.theBcType_ = PERIODIC_BC;
  const YAML::Node& periodicUserData = node["periodic_user_data"];
  periodicUserData >> periodicBC.userData_;
}

void operator >> (const YAML::Node& node, NonConformalBoundaryConditionData& nonConformalBC) {
  node["non_conformal_boundary_condition"] >> nonConformalBC.bcName_;
  node["target_name"] >> nonConformalBC.masterSlave_;
  nonConformalBC.targetName_ = nonConformalBC.masterSlave_.master_ + "_" + nonConformalBC.masterSlave_.slave_;
  nonConformalBC.theBcType_ = NON_CONFORMAL_BC;
  const YAML::Node& nonConformalUserData = node["non_conformal_user_data"];
  nonConformalUserData >> nonConformalBC.userData_;
}

void operator >> (const YAML::Node& node, MeshInput& meshInput) {
  node["mesh_name"] >> meshInput.meshName_;
}

void operator >> (const YAML::Node& node, ConstantInitialConditionData& constIC)
{
  constIC.theIcType_ = CONSTANT_UD;
  node["constant"] >> constIC.icName_;
  const YAML::Node & targets = node["target_name"];
  if (targets.Type() == YAML::NodeType::Scalar)
  {
    constIC.targetNames_.resize(1);
    targets >> constIC.targetNames_[0];
    NaluEnv::self().naluOutputP0() << "constant IC: name: " << constIC.icName_ << " , target[" << 0 << "] = "
              << constIC.targetNames_[0] << std::endl;
    if (constIC.targetNames_[0].find(',') != std::string::npos)
      throw std::runtime_error("In " + constIC.icName_ +
                               " found ',' in target name - you must enclose in '[...]' for multiple targets");
  }
  else
  {
    constIC.targetNames_.resize(targets.size());
    for (size_t i=0; i < targets.size(); ++i)
    {
      targets[i] >> constIC.targetNames_[i];
      if (constIC.root()->debug())
        NaluEnv::self().naluOutputP0() << "constant IC: name: " << constIC.icName_ << " , target[" << i << "] = "
                  << constIC.targetNames_[i] << std::endl;
    }
  }

  const YAML::Node & value_node = node["value"];
  size_t value_size = value_node.size();
  constIC.fieldNames_.resize(value_size);
  constIC.data_.resize(value_size);
  if (constIC.root()->debug())
  {
    NaluEnv::self().naluOutputP0() << "fieldNames_.size()= " << constIC.fieldNames_.size()
              << " value.size= " << constIC.data_.size() << std::endl;
  }
  size_t jv = 0;
  for (YAML::Iterator i = value_node.begin(); i != value_node.end(); ++i,++jv) {
    const YAML::Node & key   = i.first();
    const YAML::Node & value = i.second();
    key >> constIC.fieldNames_[jv] ;
    size_t nvals = value.size();
    if (nvals)
    {
      constIC.data_[jv].resize(nvals);
      for (size_t iv=0; iv < nvals; ++iv)
      {
        value[iv] >> constIC.data_[jv][iv];
        if (constIC.root()->debug())
        {
          NaluEnv::self().naluOutputP0() << "fieldNames_= " << constIC.fieldNames_[jv] << " value= "
                    << constIC.data_[jv][iv] << std::endl;
        }
      }
    }
    else
    {
      constIC.data_[jv].resize(1);
      value >> constIC.data_[jv][0];
      if (constIC.root()->debug())
      {
        NaluEnv::self().naluOutputP0() << "fieldNames_= " << constIC.fieldNames_[jv] << " value= " << constIC.data_[jv][0] << std::endl;
      }
    }
  }
}

void operator >> (const YAML::Node& node, UserFunctionInitialConditionData& fcnIC)
{
  fcnIC.theIcType_ = FUNCTION_UD;
  node["user_function"] >> fcnIC.icName_;
  const YAML::Node & targets = node["target_name"];
  if (targets.Type() == YAML::NodeType::Scalar) {
    fcnIC.targetNames_.resize(1);
    targets >> fcnIC.targetNames_[0];
  }
  else
  {
    fcnIC.targetNames_.resize(targets.size());
    for (size_t i=0; i < targets.size(); ++i) {
      targets[i] >> fcnIC.targetNames_[i];
    }
  }

  // extract function name and parameters
  if (expect_map( node, "user_function_name", false)) {
    node["user_function_name"] >> fcnIC.functionNames_;
  }

  if (expect_map( node, "user_function_parameters", true)) {
     node["user_function_parameters"] >> fcnIC.functionParams_;
  }

}

void operator >> (const YAML::Node& node, std::map<std::string,double> & mapName)
{
  for (YAML::Iterator i = node.begin(); i != node.end(); ++i)
  {
    const YAML::Node & key = i.first();
    const YAML::Node & value = i.second();
    std::string stringName;
    key >> stringName;
    double data;
    value >> data;
    mapName[stringName] = data;
  }
}

void operator >> (const YAML::Node& node, std::map<std::string,std::string> & mapName)
{
  for (YAML::Iterator i = node.begin(); i != node.end(); ++i)
  {
    const YAML::Node & key = i.first();
    const YAML::Node & value = i.second();
    std::string stringName;
    key >> stringName;
    std::string data;
    value >> data;
    mapName[stringName] = data;
  }
}

void operator >> (const YAML::Node& node, std::map<std::string,std::vector<std::string> >& mapName)
{
  for (YAML::Iterator i = node.begin(); i != node.end(); ++i)
  {
    const YAML::Node & key = i.first();
    const YAML::Node & targets = i.second();
    std::string stringName;
    key >> stringName;
    
    std::vector<std::string> &vecOfStrings = mapName[stringName];
    std::string theName;
    if ( targets.Type() == YAML::NodeType::Scalar ) {
      targets >> theName;
      vecOfStrings.push_back(theName);
    }
    else {
      for (size_t it=0; it < targets.size(); ++it) {
	targets[it] >> theName;
	vecOfStrings.push_back(theName);
      }
    }
  }
}

const YAML::Node *
expect_type(const YAML::Node& node, const std::string& key, YAML::NodeType::value type, bool optional)
{
  static std::string types[] = {"Null", "Scalar", "Sequence", "Map"};
  const YAML::Node *value = node.FindValue(key);
  std::ostringstream err_msg;
  if (!optional && !value)
    {
      if (!NaluEnv::self().parallel_rank()) {
        err_msg << "Error: parsing expected required value " << key << " but it was not found at"
                << NaluParsingHelper::line_info(node)
                << " for Node= " << std::endl;
        NaluParsingHelper::emit(err_msg, node);
        std::cout << err_msg.str() << std::endl;
      }
      throw std::runtime_error("Error: parsing");
    }
  if (value && (value->Type() != type))
    {
      if (!NaluEnv::self().parallel_rank()) {
        err_msg << "Error: parsing expected type " << types[type] << " got type = " << types[value->Type()]
                << " for key= " << key
                << " at " << NaluParsingHelper::line_info(node)
                << " node= " << std::endl;
        NaluParsingHelper::emit(err_msg, node);
        err_msg << "Check indentation of input file.";
        std::cout << err_msg.str() << std::endl;
      }
      throw std::runtime_error("Error: parsing - Check indentation of input file.");
    }
  return value;
}

const YAML::Node *
expect_null(const YAML::Node& node, const std::string& key, bool optional)
{
  return expect_type(node, key, YAML::NodeType::Null, optional);
}
const YAML::Node *
expect_scalar(const YAML::Node& node, const std::string& key, bool optional)
{
  return expect_type(node, key, YAML::NodeType::Scalar, optional);
}
const YAML::Node *
expect_sequence(const YAML::Node& node, const std::string& key, bool optional)
{
  return expect_type(node, key, YAML::NodeType::Sequence, optional);
}
const YAML::Node *
expect_map(const YAML::Node& node, const std::string& key, bool optional)
{
  return expect_type(node, key, YAML::NodeType::Map, optional);
}

} // namespace nalu
} // namespace Sierra
