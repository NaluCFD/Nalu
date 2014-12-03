/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TimeIntegrator_h
#define TimeIntegrator_h

#include <Enums.h>
// yaml for parsing..
#include <yaml-cpp/yaml.h>

#include <vector>
#include <string>


namespace sierra{
namespace nalu{

class Realm;
class Simulation;

class TimeIntegrator
{
public:

  TimeIntegrator(Simulation& sim);
  ~TimeIntegrator();

  void load(const YAML::Node & node) ;

  void breadboard();

  void initialize();
  Simulation *root();
  Simulation *parent();

  void integrate_realm();
  void provide_mean_norm();
  bool simulation_proceeds();
  Simulation& sim_;

  double totalSimTime_;
  double currentTime_;
  double timeStepFromFile_;
  double timeStepN_;
  double timeStepNm1_;
  double gamma1_;
  double gamma2_;
  double gamma3_;
  int timeStepCount_;
  int maxTimeStepCount_;
  bool secondOrderTimeAccurate_;
  bool adaptiveTimeStep_;
  bool terminateBasedOnTime_;
  int nonlinearIterations_;

  std::string name_;

  std::vector<std::string> realmNamesVec_;

  std::vector<Realm*> realmVec_;

  double get_time_step(
    const NaluState &theState = NALU_STATE_N);
  double get_current_time();
  double get_gamma1();
  double get_gamma2();
  double get_gamma3();
  int get_time_step_count();
  void compute_gamma();
 
};

} // namespace nalu
} // namespace Sierra

#endif
