/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PolynomialPropertyEvaluator_h
#define PolynomialPropertyEvaluator_h

#include <property_evaluator/PropertyEvaluator.h>
#include <FieldTypeDef.h>

#include <string>
#include <map>
#include <vector>

namespace stk { 
namespace mesh { 
struct Entity; 
} 
}

namespace sierra{
namespace nalu{

class ReferencePropertyData;

class PolynomialPropertyEvaluator : public PropertyEvaluator
{
public:

  PolynomialPropertyEvaluator(
      const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
      const std::map<std::string, std::vector<double> > &lowPolynomialCoeffsMap,
      const std::map<std::string, std::vector<double> > &highPolynomialCoeffsMap,
      double universalR);
  virtual ~PolynomialPropertyEvaluator();
  
  virtual double execute(
      double *indVarList,
      stk::mesh::Entity node) = 0;
  
  const double universalR_;
  const size_t ykVecSize_;
  const double TlowHigh_;

  std::vector<double> mw_;
  std::vector<std::vector<double> > lowPolynomialCoeffs_;
  std::vector<std::vector<double> > highPolynomialCoeffs_;

};

} // namespace nalu
} // namespace Sierra

#endif
