/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef BurkeSchumannAlgorithm_h
#define BurkeSchumannAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

#include <string>
#include <map>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;
class ReferencePropertyData;

class BurkeSchumannAlgorithm : public Algorithm
{
public:

  BurkeSchumannAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap);
  virtual ~BurkeSchumannAlgorithm() {}
  virtual void execute();

  const size_t massFracSize_;
  ScalarFieldType *mixtureFraction_;
  GenericFieldType *massFraction_;
  
  double zStoich_;
  const int fuelId_;
  const int oxidizerId_;
  const int carbonDioxideId_;
  const int waterId_;
  const int diluentId_;
  const double M_;
  const double N_;

  std::vector<double> mwVec_;
  std::vector<double> stoichiometryVec_;
  std::vector<double> primaryVec_;
  std::vector<double> secondaryVec_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
