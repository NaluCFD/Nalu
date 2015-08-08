/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef OutputInfo_h
#define OutputInfo_h

#include <NaluParsing.h>

#include <string>
#include <set>

namespace sierra{
namespace nalu{

class OutputInfo
{
public:
  
  OutputInfo();
  ~OutputInfo();
  
  void load(const YAML::Node & node);
  
  std::string outputDBName_;
  int outputFreq_;
  bool outputNodeSet_; 
  int serializedIOGroupSize_;
  bool hasOutputBlock_;
  bool hasRestartBlock_;
  bool activateRestart_;
  bool meshAdapted_;
  double restartTime_;
  std::string restartDBName_;
  int restartFreq_;
  int restartMaxDataBaseStepSize_;
  bool restartNodeSet_;
  std::set<std::string> outputFieldNameSet_;
  std::set<std::string> restartFieldNameSet_;

};

} // namespace nalu
} // namespace Sierra

#endif
