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

namespace stk {
  namespace mesh {
    class FieldBase;
    class Part;
    typedef std::vector<Part*> PartVector;
  }
}

namespace sierra{
namespace nalu{

class AveragingInfo
{
public:

  AveragingInfo();
  ~AveragingInfo();

  // name of this block
  std::string name_;

  // specialty options
  bool computeReynoldsStress_;
  bool computeTke_;
  bool computeFavreStress_;
  bool computeFavreTke_;
  bool computeResolvedStress_{false};
  bool computeSFSStress_{false};
  bool computeVorticity_;
  bool computeQcriterion_;
  bool computeLambdaCI_;
  bool computeMeanResolvedKe_;

  // Temperature stresses
  bool computeTemperatureSFS_{false};
  bool computeTemperatureResolved_{false};
  
  // vector of part names, e.g., block_1, surface_2
  std::vector<std::string> targetNames_;

  // vector of parts
  stk::mesh::PartVector partVec_;

  // vector of favre/reynolds fields
  std::vector<std::string> favreFieldNameVec_;
  std::vector<std::string> reynoldsFieldNameVec_;
  std::vector<std::string> resolvedFieldNameVec_;
  std::vector<std::string> movingAvgFieldNameVec_;


  // vector of pairs of fields
  std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > favreFieldVecPair_;
  std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > reynoldsFieldVecPair_;
  std::vector<std::pair<stk::mesh::FieldBase *, stk::mesh::FieldBase *> > resolvedFieldVecPair_;

  // sizes for each
  std::vector<unsigned> favreFieldSizeVec_;
  std::vector<unsigned> reynoldsFieldSizeVec_;
  std::vector<unsigned> resolvedFieldSizeVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
