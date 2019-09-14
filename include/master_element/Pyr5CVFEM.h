/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Pyr5CVFEM_h
#define Pyr5CVFEM_h

#include <master_element/MasterElement.h>

#include <AlgTraits.h>

// NGP-based includes
#include "SimdInterface.h"
#include "KokkosInterface.h"

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <array>

namespace stk {
  struct topology;
}

namespace sierra{
namespace nalu{

struct ElementDescription;
class MasterElement;


// Pyramid 5 subcontrol volume
class PyrSCV : public MasterElement
{
public:
  using AlgTraits = AlgTraitsPyr5;


  PyrSCV();
  virtual ~PyrSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType*>& vol);

  void grad_op(
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gradop,
    SharedMemView<DoubleType***>& deriv);

  void shifted_grad_op(
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gradop,
    SharedMemView<DoubleType***>& deriv);

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void pyr_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  void shifted_pyr_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);
};

// Pyramid 5 subcontrol surface
class PyrSCS : public MasterElement
{
public:
  using AlgTraits = AlgTraitsPyr5;

  PyrSCS();
  virtual ~PyrSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType**>& areav);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gradop,
    SharedMemView<DoubleType***>& deriv);

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gradop,
    SharedMemView<DoubleType***>& deriv);

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void gij( 
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gupper,
    SharedMemView<DoubleType***>& glower,
    SharedMemView<DoubleType***>& deriv);

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  const int * scsIpEdgeOrd();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void pyr_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  void shifted_pyr_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  void
  general_shape_fcn(const int numIp, const double* isoParCoord, double* shpfc)
  {
    pyr_shape_fcn(numIp, isoParCoord, shpfc);
  }

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double *error);

  void face_grad_op(
    int face_ordinal,
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gradop) final;

  void shifted_face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void shifted_face_grad_op(
    int face_ordinal,
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gradop) final;

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double *error);

  const int* side_node_ordinals(int sideOrdinal) final;

  double parametric_distance(const std::array<double,3>& x);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

};

} // namespace nalu
} // namespace Sierra

#endif
