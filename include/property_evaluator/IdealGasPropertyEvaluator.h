/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef IdealGasPropertyEvaluator_h
#define IdealGasPropertyEvaluator_h

#include <property_evaluator/PropertyEvaluator.h>
#include<FieldTypeDef.h>

#include <vector>

namespace stk { namespace mesh { struct Entity; } }

namespace stk {
namespace mesh {
class MetaData;
}
}

namespace sierra{
namespace nalu{

class IdealGasPrefTYkrefPropertyEvaluator : public PropertyEvaluator
{
public:

  IdealGasPrefTYkrefPropertyEvaluator(
    double pRef,
    double universalR,
    std::vector<std::pair<double, double> > mwMassFracVec);
  virtual ~IdealGasPrefTYkrefPropertyEvaluator();
  
  double execute(
    double *indVarList,
    stk::mesh::Entity node);
  
  const double pRef_;
  const double R_;
  double mwRef_;
};

class IdealGasPrefTYkPropertyEvaluator : public PropertyEvaluator
{
public:

  IdealGasPrefTYkPropertyEvaluator(
      double pRef,
      double universalR,
      std::vector<double>mwVec,
      stk::mesh::MetaData &metaData);
  virtual ~IdealGasPrefTYkPropertyEvaluator();
  
  double execute(
      double *indVarList,
      stk::mesh::Entity node);
  
  double compute_mw(
      const double *yk);
  
  // reference quantities
  const double pRef_;
  const double R_;  

  // field definition and extraction,
  GenericFieldType *massFraction_;

  // reference mw vector; size and declaration
  size_t mwVecSize_;
  std::vector<double> mwVec_;
  
};

class IdealGasPTYkrefPropertyEvaluator : public PropertyEvaluator
{
public:

  IdealGasPTYkrefPropertyEvaluator(
      double universalR,
      std::vector<std::pair<double, double> > mwMassFracVec,
      stk::mesh::MetaData &metaData);
  virtual ~IdealGasPTYkrefPropertyEvaluator();

  double execute(
      double *indVarList,
      stk::mesh::Entity node);

  // reference quantities
  const double R_;

  // molecular weight computed at construction time
  double mwRef_;

  // field definition and extraction,
  ScalarFieldType *pressure_;
};

class IdealGasPrefTrefYkPropertyEvaluator : public PropertyEvaluator
{
public:

  IdealGasPrefTrefYkPropertyEvaluator(
      double pRef,
      double tRef,
      double universalR,
      std::vector<double>mwVec,
      stk::mesh::MetaData &metaData);
  virtual ~IdealGasPrefTrefYkPropertyEvaluator();
  
  double execute(
      double *indVarList,
      stk::mesh::Entity node);
  
  double compute_mw(
      const double *yk);
  
  // reference quantities
  const double pRef_;
  const double tRef_;
  const double R_;  

  // field definition and extraction,
  GenericFieldType *massFraction_;

  // reference mw vector; size and declaration
  size_t mwVecSize_;
  std::vector<double> mwVec_;
  
};

class IdealGasPrefTrefYkrefPropertyEvaluator : public PropertyEvaluator
{
public:

  IdealGasPrefTrefYkrefPropertyEvaluator(
    double pRef,
    double tRef,
    double universalR,
    std::vector<std::pair<double, double> > mwMassFracVec);
  virtual ~IdealGasPrefTrefYkrefPropertyEvaluator();
  
  double execute(
    double *indVarList,
    stk::mesh::Entity node);
  
  double compute_mw(
    const double *yk);
  
  // reference quantities
  const double pRef_;
  const double tRef_;
  const double R_;  
  double mwRef_;
};

} // namespace nalu
} // namespace Sierra

#endif
