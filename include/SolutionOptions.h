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

namespace sierra{
namespace nalu{

class MeshMotionInfo;

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
  double hybridDefault_;
  double alphaDefault_;
  double alphaUpwDefault_;
  double upwDefault_;
  double lamScDefault_;
  double turbScDefault_;
  double turbPrDefault_;
  bool nocDefault_;
  std::string tanhFormDefault_;
  double tanhTransDefault_;
  double tanhWidthDefault_;
  double referenceDensity_;
  double referenceTemperature_;
  double thermalExpansionCoeff_;
  double stefanBoltzmann_;
  double nearestFaceEntrain_;
  double includeDivU_;
  bool mdotInterpRhoUTogether_;
  bool isTurbulent_;
  TurbulenceModel turbulenceModel_;
  bool meshMotion_;
  bool meshDeformation_;
  bool externalMeshDeformation_;
  bool activateUniformRefinement_;
  bool uniformRefineSaveAfter_;
  std::vector<int> refineAt_;
  bool activateAdaptivity_;
  ErrorIndicatorType errorIndicatorType_;
  int adaptivityFrequency_;
  bool useMarker_;
  double refineFraction_;
  double unrefineFraction_;
  double physicalErrIndCriterion_;
  double physicalErrIndUnrefCriterionMultipler_;
  double maxRefinementNumberOfElementsFraction_;
  bool adapterExtraOutput_;
  bool useAdapter_;
  int maxRefinementLevel_;
  bool ncAlgGaussLabatto_;
  bool ncAlgUpwindAdvection_;
  bool ncAlgIncludePstab_;
  bool ncAlgDetailedOutput_;
  bool ncAlgCurrentNormal_;
  bool cvfemShiftMdot_;
  bool cvfemShiftPoisson_;
  bool cvfemReducedSensPoisson_;
  double inputVariablesRestorationTime_;
  bool inputVariablesInterpolateInTime_;
  bool consistentMMPngDefault_;
  bool useConsolidatedSolverAlg_;
  bool eigenvaluePerturb_;
  double eigenvaluePerturbDelta_;
  int eigenvaluePerturbBiasTowards_;
  double eigenvaluePerturbTurbKe_;
  double earthAngularVelocity_;
  double latitude_;

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

  // read any fields from input files
  std::map<std::string, std::string> inputVarFromFileMap_;

  // mesh motion
  std::map<std::string, MeshMotionInfo *> meshMotionInfoMap_;

  std::vector<double> gravity_;

  // Coriolis source term
  std::vector<double> eastVector_;
  std::vector<double> northVector_;

  std::string name_;

};

} // namespace nalu
} // namespace Sierra

#endif
