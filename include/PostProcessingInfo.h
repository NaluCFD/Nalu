/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PostProcessingInfo_h
#define PostProcessingInfo_h

#include <NaluParsing.h>

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class PostProcessingData;

class PostProcessingInfo
{

public:
  
  PostProcessingInfo();
  ~PostProcessingInfo();
  
  void load(const YAML::Node & y_node);
  
  std::vector<PostProcessingData *> ppDataVec_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
