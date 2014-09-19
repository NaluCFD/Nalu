/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AveragingInfo.h>

#include <NaluParsing.h>

// basic c++
#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AveragingInfo - holder for averaging information held at realm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AveragingInfo::AveragingInfo() 
  : currentTimeFilter_(0.0),
    timeFilterInterval_(1.0e8),
    forcedReset_(false),
    processAveraging_(false)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AveragingInfo::~AveragingInfo()
{
  // nothing to do
}


//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
AveragingInfo::load(
  const YAML::Node & y_node)
{

  // output for results
  const YAML::Node *y_average = y_node.FindValue("turbulence_averaging");
  if (y_average)
  {    
    processAveraging_ = true;
    get_if_present(*y_average, "forced_reset", forcedReset_, forcedReset_);
    get_if_present(*y_average, "time_filter_interval", timeFilterInterval_, timeFilterInterval_);

    // reynolds
    const YAML::Node *y_reynolds = y_average->FindValue("reynolds_averaged_variables");
    if (y_reynolds)
    {
      size_t varSize = y_reynolds->size();
      for (size_t ioption = 0; ioption < varSize; ++ioption)
      {
        const YAML::Node & y_var = (*y_reynolds)[ioption];
        std::string fieldName;
        y_var >> fieldName;
        if ( fieldName != "density" )
          reynoldsFieldNameVec_.push_back(fieldName);
      }
    }

    // favre
    const YAML::Node *y_favre = y_average->FindValue("favre_averaged_variables");
    if (y_favre)
    {
      size_t varSize = y_favre->size();
      for (size_t ioption = 0; ioption < varSize; ++ioption)
      {
        const YAML::Node & y_var = (*y_favre)[ioption];
        std::string fieldName;
        y_var >> fieldName;
        if ( fieldName != "density")
          favreFieldNameVec_.push_back(fieldName);
      }
    }
  }
 
}


} // namespace nalu
} // namespace Sierra
