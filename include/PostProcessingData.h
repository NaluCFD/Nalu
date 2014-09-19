/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PostProcessingData_h
#define PostProcessingData_h

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class PostProcessingData {

 public:
  PostProcessingData() : type_("na"), physics_("na"), outputFileName_("na"), frequency_(10){}
  ~PostProcessingData() {}
  
  std::string type_;
  std::string physics_;
  std::string outputFileName_;
  int frequency_;
  std::vector<double> parameters_;
  std::vector<std::string> targetNames_;
};
 
} // namespace nalu
} // namespace Sierra

#endif
