/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <OutputInfo.h>
#include <NaluEnv.h>
#include <NaluParsing.h>

// basic c++
#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// OutputInfo - holder for output information held at realm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
OutputInfo::OutputInfo() 
  : outputDBName_("output.e"),
    outputFreq_(1),
    outputNodeSet_(false),
    serializedIOGroupSize_(0),
    hasOutputBlock_(false),
    hasRestartBlock_(false),
    activateRestart_(false),
    meshAdapted_(false),
    restartTime_(0.0),
    restartDBName_("restart.rst"),
    restartFreq_(500),
    restartMaxDataBaseStepSize_(100000),
    restartNodeSet_(true)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
OutputInfo::~OutputInfo()
{
  // nothing to do
}


//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
OutputInfo::load(
  const YAML::Node & y_node)
{

  const YAML::Node *y_output = y_node.FindValue("output");
  if (y_output)
  {
    // user desires output
    hasOutputBlock_ = true;

    // output data base name
    get_if_present(*y_output, "output_data_base_name", outputDBName_, outputDBName_);
    
    // output frequency
    get_if_present(*y_output, "output_frequency", outputFreq_, outputFreq_);

    // determine if we want nodeset output
    get_if_present(*y_output, "output_node_set", outputNodeSet_, outputNodeSet_);

    // serialize io...
    {
      get_if_present(*y_output, "serialized_io_group_size", serializedIOGroupSize_, serializedIOGroupSize_);
      if (serializedIOGroupSize_) {
        NaluEnv::self().naluOutputP0() << "Info: found non-zero serialized_io_group_size in input file= " << serializedIOGroupSize_ << std::endl;
      }
    }

    const YAML::Node *y_vars = y_output->FindValue("output_variables");
    if (y_vars)
    {
      size_t varSize = y_vars->size();
      for (size_t ioption = 0; ioption < varSize; ++ioption)
      {
        const YAML::Node & y_var = (*y_vars)[ioption];
        std::string fieldName;
        y_var >> fieldName;
        outputFieldNameSet_.insert(fieldName);
      }
    }
  }
  
  // output for restart
  const YAML::Node *y_restart = y_node.FindValue("restart");
  if (y_restart)
  {    
    // some sort of intent to manage a restart event - either clean or restart run
    hasRestartBlock_ = true;

    // restart output data base name
    get_if_present(*y_restart, "restart_data_base_name", restartDBName_, restartDBName_);
    
    // restart output frequency
    get_if_present(*y_restart, "restart_frequency", restartFreq_, restartFreq_);
    
    // determine if we want nodeset restart output
    get_if_present(*y_restart, "restart_node_set", restartNodeSet_, restartNodeSet_);
    
    // max data base size for restart
    get_if_present(*y_restart, "max_data_base_step_size", restartMaxDataBaseStepSize_, restartMaxDataBaseStepSize_);

    // check to see if restart is active for this run
    if ( y_restart->FindValue("restart_time") ) {
      activateRestart_ = true;
      *(y_restart->FindValue("restart_time")) >> restartTime_;
    }
    
    const YAML::Node *y_vars = y_restart->FindValue("restart_variables");
    if (y_vars) {
      NaluEnv::self().naluOutputP0() << "Restart variable specification has been deprecated" << std::endl;
    }
  }

}


} // namespace nalu
} // namespace Sierra
