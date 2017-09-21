/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef CoriolisSrc_h
#define CoriolisSrc_h

#include <vector>

namespace sierra{
namespace nalu{

class SolutionOptions;

class CoriolisSrc {
public:

  CoriolisSrc(const SolutionOptions& solnOpts);
  virtual ~CoriolisSrc() {}

  int nDim_;
  double earthAngularVelocity_;
  double latitude_;
  double sinphi_;
  double cosphi_;
  double corfac_;
  double Jxy_, Jxz_, Jyz_;
  double pi_;
  std::vector<double> eastVector_;
  std::vector<double> northVector_;
  std::vector<double> upVector_;

};

} // namespace nalu
} // namespace Sierra

#endif
