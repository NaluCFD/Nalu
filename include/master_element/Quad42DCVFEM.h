/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Quad42DCVFEM_h   
#define Quad42DCVFEM_h   

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

namespace sierra{
namespace nalu{

// 2D Quad 4 subcontrol volume
class Quad42DSCV : public MasterElement
{
public:
  using Traits = AlgTraitsQuad4_2D;
  Quad42DSCV();
  virtual ~Quad42DSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    SharedMemView<DoubleType**> &coords,
    SharedMemView<DoubleType*> &vol);

  void grad_op(
    SharedMemView<DoubleType** >& coords,
    SharedMemView<DoubleType***>& gradop,
    SharedMemView<DoubleType***>& deriv);

  void shifted_grad_op(
    SharedMemView<DoubleType** >& coords,
    SharedMemView<DoubleType***>& gradop,
    SharedMemView<DoubleType***>& deriv);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void quad_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn) ;
};

// 2D Quad 4 subcontrol surface
class Quad42DSCS : public MasterElement
{
public:
  using Traits = AlgTraitsQuad4_2D;
  Quad42DSCS();
  virtual ~Quad42DSCS();

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
    SharedMemView<DoubleType** >& coords,
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
    SharedMemView<DoubleType** >& coords,
    SharedMemView<DoubleType***>& gradop,
    SharedMemView<DoubleType***>& deriv);

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    int face_ordinal,
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gradop) final;

  void face_grad_op(
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

  void shifted_face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij( 
    SharedMemView<DoubleType** >& coords,
    SharedMemView<DoubleType***>& gupper,
    SharedMemView<DoubleType***>& glower,
    SharedMemView<DoubleType***>& deriv);

  void gij(
     const double *coords,
     double *gupperij,
     double *gij,
     double *deriv);

  const int * adjacentNodes();

  const int * scsIpEdgeOrd();

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void quad_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn) ;

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);
  
  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  const int* side_node_ordinals(int sideOrdinal) final;
private :
  void face_grad_op(
    const int face_ordinal,
    const bool shifted,
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType***>& gradop);

};

} // namespace nalu
} // namespace Sierra

#endif
