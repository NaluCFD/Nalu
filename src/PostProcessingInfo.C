/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <PostProcessingInfo.h>
#include <PostProcessingData.h>
#include <NaluParsing.h>

// basic c++
#include <stdexcept>

namespace sierra{
namespace nalu{


//==========================================================================
// Class Definition
//==========================================================================
// PostProcessingInfo - holder for post processing information held at realm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
PostProcessingInfo::PostProcessingInfo() 
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
PostProcessingInfo::~PostProcessingInfo()
{
  std::vector<PostProcessingData *>::iterator ii;
  for( ii=ppDataVec_.begin(); ii!=ppDataVec_.end(); ++ii )
    delete *ii;
}
  
//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
PostProcessingInfo::load(
    const YAML::Node & y_node)
{

  // post processing on?
  const bool optional = true;
  const YAML::Node *y_pp = expect_sequence(y_node, "post_processing", optional);
  if (y_pp) {

    for (size_t itype = 0; itype < y_pp->size(); ++itype) {

      // extract the particular type
      const YAML::Node &y_type = (*y_pp)[itype];
  
      // create the data
      PostProcessingData *ppData = new PostProcessingData();      

      // push the data to vector
      ppDataVec_.push_back(ppData);

      // save the type; surface/node/volume
      y_type["type"] >> ppData->type_;
      
      // physics; traction, etc.
      if ( y_type.FindValue("physics") )
        y_type["physics"] >> ppData->physics_;
      else
        throw std::runtime_error("parser error PostProcessings::load: no physics type specified");
      
      // outfile file
      if ( y_type.FindValue("output_file_name") )
        y_type["output_file_name"] >> ppData->outputFileName_;
      else
        throw std::runtime_error("parser error PostProcessings::load: no output file specified");
      
      // frequency
      if ( y_type.FindValue("frequency"))
        y_type["frequency"] >> ppData->frequency_;

      // misc parameters
      if ( y_type.FindValue("parameters") ) {

        // extract the value(s)
        const YAML::Node & targets = y_type["parameters"];
        if (targets.Type() == YAML::NodeType::Scalar) {
          ppData->parameters_.resize(1);
          targets >> ppData->parameters_[0];
        }
        else {
          ppData->parameters_.resize(targets.size());
          for (size_t i=0; i < targets.size(); ++i) {
            targets[i] >> ppData->parameters_[i];
          }
        }
      }
      
      // extract the target(s)
      const YAML::Node & targets = y_type["target_name"];
      if (targets.Type() == YAML::NodeType::Scalar) {
        ppData->targetNames_.resize(1);
        targets >> ppData->targetNames_[0];
      }
      else {
        ppData->targetNames_.resize(targets.size());
        for (size_t i=0; i < targets.size(); ++i) {
          targets[i] >> ppData->targetNames_[i];
        }
      }
    }
  }
 
}


} // namespace nalu
} // namespace Sierra
