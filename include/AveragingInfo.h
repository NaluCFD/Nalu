/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AveragingInfo_h
#define AveragingInfo_h

#include <NaluParsing.h>

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class AveragingInfo
{
public:
  
  AveragingInfo();
  ~AveragingInfo();
  
  void load(const YAML::Node & node);
  
  double currentTimeFilter_; /* providd by restart */
  double timeFilterInterval_; /* user supplied */
  bool forcedReset_; /* allows forhard reset */
  bool processAveraging_;
  
  std::vector<std::string> favreFieldNameVec_;
  std::vector<std::string> reynoldsFieldNameVec_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
