/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Tri6FEM_h
#define Tri6FEM_h

#include<master_element/MasterElement.h>

namespace sierra{
namespace nalu{

// Tet 10 FEM
class Tri6FEM : public MasterElement
{
public:

  Tri6FEM();
  virtual ~Tri6FEM();
    
  // NGP-ready methods first
  void determinant_fem(
    SharedMemView<DoubleType**>&coords,
    SharedMemView<DoubleType***>&deriv,
    SharedMemView<DoubleType*>&det_j) final;

  void normal_fem(
    SharedMemView<DoubleType**>&coords,
    SharedMemView<DoubleType***>&deriv,
    SharedMemView<DoubleType**>&normal) final;

  void shape_fcn(
    SharedMemView<DoubleType**> &shpfc) final;

  void shape_fcn(
    double *shpfc);

  void tri6_shape_fcn(
    const int  &npts,
    const double *isoParCoord,
    SharedMemView<DoubleType**> &shpfc);

  void tri6_shape_fcn(
    const int  &npts,
    const double *isoParCoord,
    double *shpfc);
};
    
} // namespace nalu
} // namespace Sierra

#endif
