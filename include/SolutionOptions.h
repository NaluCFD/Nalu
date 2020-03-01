/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SolutionOptions_h
#define SolutionOptions_h

#include <NaluParsing.h>
#include <Enums.h>

// standard c++
#include <string>
#include <map>
#include <utility>
#include <memory>

namespace sierra{
namespace nalu{

class MeshMotionInfo;
struct FixPressureAtNodeInfo;

enum ErrorIndicatorType {
  EIT_NONE                = 0,
  EIT_PSTAB               = 1 << 1,
  EIT_LIMITER             = 1 << 2 ,
  EIT_SIMPLE_BASE         = 1 << 3,
  EIT_SIMPLE_VORTICITY    = EIT_SIMPLE_BASE + (1 << 4),
  EIT_SIMPLE_VORTICITY_DX = EIT_SIMPLE_BASE + (1 << 5),
  EIT_SIMPLE_DUDX2        = EIT_SIMPLE_BASE + (1 << 6)
};


class SolutionOptions
{
public:

  SolutionOptions();
  ~SolutionOptions();

  void load(const YAML::Node & node);
  void initialize_turbulence_constants();

  inline bool has_mesh_motion() const { return meshMotion_; }

  inline bool has_mesh_deformation() const
  {
    return externalMeshDeformation_ | meshDeformation_;
  }

  inline bool does_mesh_move() const
  {
    return has_mesh_motion() | has_mesh_deformation();
  }

  inline std::string get_coordinates_name() const
  {
    return ( (meshMotion_ | meshDeformation_ | externalMeshDeformation_ | initialMeshDisplacement_) 
	     ? "current_coordinates" : "coordinates");    
  }
  
  double get_alpha_factor(const std::string&) const;

  double get_alpha_upw_factor(const std::string&) const;

  double get_upw_factor(const std::string&) const;

  bool primitive_uses_limiter(const std::string&) const;

  bool get_shifted_grad_op(const std::string&) const;
  
  bool get_skew_symmetric(const std::string&) const;

  std::vector<double> get_gravity_vector(const unsigned nDim) const;
 
  double get_turb_model_constant(
    TurbulenceModelConstant turbModelEnum) const;
  
  double get_turb_prandtl(const std::string &dofName) const;

  bool get_noc_usage(const std::string &dofName) const;

  void set_consolidated_bc_solver_alg();

  double hybridDefault_;
  double alphaDefault_;
  double alphaUpwDefault_;
  double upwDefault_;
  double lamScDefault_;
  double turbScDefault_;
  double turbPrDefault_;
  bool nocDefault_;
  bool shiftedGradOpDefault_;
  bool skewSymmetricDefault_;
  std::string tanhFormDefault_;
  double tanhTransDefault_;
  double tanhWidthDefault_;
  double referenceDensity_;
  double referenceTemperature_;
  double thermalExpansionCoeff_;
  double stefanBoltzmann_;
  double includeDivU_;
  bool isTurbulent_;
  TurbulenceModel turbulenceModel_;
  bool meshMotion_;
  bool meshDeformation_;
  bool externalMeshDeformation_;
  bool initialMeshDisplacement_;
  bool errorIndicatorActive_;
  ErrorIndicatorType errorIndicatorType_;
  int errorIndicatorFrequency_;
  bool ncAlgGaussLabatto_;
  bool ncAlgUpwindAdvection_;
  bool ncAlgIncludePstab_;
  bool ncAlgDetailedOutput_;
  bool ncAlgCoincidentNodesErrorCheck_;
  bool ncAlgCurrentNormal_;
  bool ncAlgPngPenalty_;
  bool cvfemShiftMdot_;
  bool cvfemReducedSensPoisson_;
  double inputVariablesRestorationTime_;
  bool inputVariablesInterpolateInTime_;
  double inputVariablesPeriodicTime_;
  bool consistentMMPngDefault_;
  bool useConsolidatedSolverAlg_;
  bool useConsolidatedBcSolverAlg_;
  bool eigenvaluePerturb_;
  double eigenvaluePerturbDelta_;
  int eigenvaluePerturbBiasTowards_;
  double eigenvaluePerturbTurbKe_;
 
  // mdot post processing
  double mdotAlgAccumulation_;
  double mdotAlgInflow_;
  double mdotAlgOpen_;
 
  // turbulence model coeffs
  std::map<TurbulenceModelConstant, double> turbModelConstantMap_;
  
  // numerics related
  std::map<std::string, double> hybridMap_;
  std::map<std::string, double> alphaMap_;
  std::map<std::string, double> alphaUpwMap_;
  std::map<std::string, double> upwMap_;
  std::map<std::string, bool> limiterMap_;
  std::map<std::string, std::string> tanhFormMap_;
  std::map<std::string, double> tanhTransMap_;
  std::map<std::string, double> tanhWidthMap_;
  std::map<std::string, bool> consistentMassMatrixPngMap_;
  std::map<std::string, bool> skewSymmetricMap_;

  // property related
  std::map<std::string, double> lamScMap_;
  std::map<std::string, double> lamPrMap_;
  std::map<std::string, double> turbScMap_;
  std::map<std::string, double> turbPrMap_;

  // source; nodal and fully integrated
  std::map<std::string, std::vector<std::string> > srcTermsMap_;
  std::map<std::string, std::vector<double> > srcTermParamMap_;
  std::map<std::string, std::vector<std::string> > elemSrcTermsMap_;
  std::map<std::string, std::vector<double> > elemSrcTermParamMap_;

  // nodal gradient
  std::map<std::string, std::string> nodalGradMap_;

  // non-orthogonal correction
  std::map<std::string, bool> nocMap_;
  
  // shifting of Laplace operator for the element-based grad_op
  std::map<std::string, bool> shiftedGradOpMap_;
  
  // read any fields from input files
  std::map<std::string, std::string> inputVarFromFileMap_;

  // mesh motion
  std::map<std::string, MeshMotionInfo *> meshMotionInfoMap_;

  // initial displacement
  std::map<std::string, MeshMotionInfo *> initialMeshDisplacementInfoMap_;

  std::vector<double> gravity_;

  // Coriolis source term
  std::vector<double> eastVector_;
  std::vector<double> northVector_;

  //! Flag indicating whether the user has requested pressure referencing
  bool needPressureReference_{false};

  std::unique_ptr<FixPressureAtNodeInfo> fixPressureInfo_;

  std::string name_;

  std::string quadType_;
  
  // allow for rho = f(P)
  bool accousticallyCompressible_;
};

} // namespace nalu
} // namespace Sierra

#endif
