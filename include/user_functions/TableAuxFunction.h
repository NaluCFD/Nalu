/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef TableAuxFunction_h
#define TableAuxFunction_h

#include <AuxFunction.h>

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class TableAuxFunction : public AuxFunction
{
public:

  TableAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    std::vector<double> params);

  virtual ~TableAuxFunction();
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
  double interp(const std::pair<double,double> &x) const;

private:
  int fieldComp_;
  int interpComp_;
  int tableOffset_;
  std::vector<std::pair<double,double> > table_;
};

} // namespace nalu
} // namespace Sierra

#endif
