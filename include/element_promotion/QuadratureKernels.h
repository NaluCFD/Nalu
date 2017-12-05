/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level nalu      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef QuadratureKernels_h
#define QuadratureKernels_h

#include <Teuchos_BLAS.hpp>

#include <vector>

namespace sierra {
namespace nalu {

  struct ElementDescription;

  class SGLQuadratureOps
  {
  public:
    SGLQuadratureOps(const ElementDescription& elem);

    void volume_2D(
      const double*  nodalValuesTensor,
      double* result);

    void volume_3D(const double*  nodalValue, double* result);

    void surface_2D(
      const double*  integrand,
      double* result,
      int line_offset);

    void surface_3D(
      const double* integrand,
      double* result,
      int face_offset);

    void surfaces_2D(const double* integrand, double* result);
    void surfaces_3D(const double* integrand, double* result);

  private:
    const Teuchos::BLAS<int,double> blas_;
    int nodes1D_;
    int nodesPerElement_;
    int numSurfaces_;
    int nodesPerFace_;
    int size3D_;

    std::vector<double> work2D_;
    std::vector<double> weightTensor_;
    std::vector<double> weightMatrix_;
    double* p_weightTensor_;
    double* p_weightMatrix_;
    double* p_work2D_;
  };

} // namespace nalu
} // namespace Sierra

#endif
