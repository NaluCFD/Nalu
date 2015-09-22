/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MaterialPropertyData_h
#define MaterialPropertyData_h

#include <Enums.h>

#include <string>
#include <vector>
#include <map>

namespace sierra{
namespace nalu{

class MaterialPropertyData {
public:
  
  MaterialPropertyData();
  ~MaterialPropertyData();
  
  MaterialPropertyType type_;
  double constValue_;
  
  // mixture fraction specifics
  double primary_;
  double secondary_;
 
  // table specifics, all single in size, all possibly required to be more general
  std::vector<std::string> indVarName_;
  std::vector<std::string> indVarTableName_;
  std::string auxVarName_;
  std::string tablePropName_;
  std::string tableAuxVarName_;

  // generic property name
  std::string genericPropertyEvaluatorName_;

  // vectors and maps
  std::map<std::string, std::vector<double> > polynomialCoeffsMap_;
  std::map<std::string, std::vector<double> > lowPolynomialCoeffsMap_;
  std::map<std::string, std::vector<double> > highPolynomialCoeffsMap_;
  std::map<std::string, double> cpConstMap_;
  std::map<std::string, double> hfConstMap_;
};

} // namespace nalu
} // namespace Sierra

#endif
