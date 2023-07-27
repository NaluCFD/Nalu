/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef Table2dAuxFunction_h
#define Table2dAuxFunction_h

#include <AuxFunction.h>

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class Table2dAuxFunction : public AuxFunction
{
public:

  Table2dAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    std::vector<double> params,
    std::vector<std::string> stringParams);

  virtual ~Table2dAuxFunction();
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
  // helper function to find the lower and upper bounds for the table
  void find_entry(
   double &x, int &low, int &high, const std::vector<double> &table) const;

private:

  // from the user
  int fieldComp_;
  int interpX_;
  int interpY_;
  size_t sizeOfX_;
  size_t sizeOfY_;
  double timeBlending_;
  std::string baseFileName_;

  // internal data structures
  std::vector<double> tableX_;
  std::vector<double> tableY_;
  std::vector<std::vector<double> > A_;
};
 
} // namespace nalu
} // namespace Sierra

#endif
