/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef InputOutputRealm_h
#define InputOutputRealm_h

#include <Realm.h>

// standard c++
#include <map>
#include <string>
#include <vector>
#include <stdint.h>
#include <utility>

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

class Realms;

class InputOutputInfo {

public:
  
 InputOutputInfo() : fieldName_("na"), fieldSize_(1), fieldType_("node_rank") {}
  ~InputOutputInfo() {}

  std::string fieldName_;
  int fieldSize_;
  std::string fieldType_;
  std::vector<std::string> targetNames_;
};
 
class InputOutputRealm : public Realm 
{  
public:

  InputOutputRealm(Realms&, const YAML::Node & node);
  ~InputOutputRealm();

  typedef size_t SizeType;

  // called by others; Realms
  void load(const YAML::Node & node);
  void breadboard() {}
  void initialize();

  // called by others; time integrator; virtual
  void populate_initial_condition() {}
  void populate_boundary_data() {}
  void boundary_data_to_state_data() {}
  void populate_derived_quantities() {}
  void evaluate_properties() {}
  void initial_work() {}
  double compute_adaptive_time_step() { return 1.0e8; }
  void swap_states() {}
  void predict_state() {}   
  void pre_timestep_work() {}  
  void output_banner() {}
  void advance_time_step() {}
  double populate_restart( double &timeStepNm1, int &timeStepCount);
  void populate_external_variables_from_input(const double currentTime);
 
  // internal calls
  void register_io_fields();

  // hold the field information
  std::vector<InputOutputInfo *> inputOutputFieldInfo_;
};

} // namespace nalu
} // namespace Sierra

#endif
