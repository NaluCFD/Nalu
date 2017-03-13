/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef TensorProductQuadratureRule_h
#define TensorProductQuadratureRule_h

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class TensorProductQuadratureRule
{
public:
  TensorProductQuadratureRule(std::string type, int polyOrder);

  TensorProductQuadratureRule(
    std::string type,
    int numQuad,
    std::vector<double>& scsLocs
  );
  ~TensorProductQuadratureRule() = default;

  const std::vector<double>& abscissae() const { return abscissae_; };
  const std::vector<double>& weights() const { return weights_; };
  const std::vector<double>& scs_loc() const { return scsLoc_; };
  const std::vector<double>& scs_end_loc() const { return scsEndLoc_; };

  int num_quad() const { return numQuad_; }
  bool is_SGL() const { return useSGL_; }

  double abscissa(unsigned j) const { return abscissae_[j]; };
  double weight(unsigned j) const { return weights_[j]; };
  double scs_loc(unsigned j) const { return scsLoc_[j]; };
  double scs_end_loc(unsigned j) const { return scsEndLoc_[j]; };

  double integration_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double integration_point_weight(
    int s1Node, int s2Node, int s3Node,
    int s1Ip, int s2Ip, int s3Ip) const;

  double integration_point_weight(
    int s1Node, int s2Node,
    int s1Ip, int s2Ip) const;

  double integration_point_weight(int s1Node, int s1Ip) const;

  double isoparametric_mapping(
    const double b,
    const double a,
    const double xi) const;

private:
  std::vector<double> scsLoc_;
  std::vector<double> scsEndLoc_;
  int numQuad_;
  std::vector<double> abscissae_;
  std::vector<double> weights_;
  bool useSGL_;
};

} // namespace nalu
} // namespace Sierra

#endif
