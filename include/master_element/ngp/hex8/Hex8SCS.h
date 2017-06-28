/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Hex8SCS_h
#define Hex8SCS_h

#include "master_element/MasterElement.h"
#include "SimdInterface.h"

#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

// Hex 8 SCS; -1.0/2 : 1.0/2 range
class Hex8SCS : public HexSCS
{
public:

  Hex8SCS();
  virtual ~Hex8SCS();

  void shape_fcn(
    SharedMemView<DoubleType**> &shpfc);

  void shifted_shape_fcn(
    SharedMemView<DoubleType**> &shpfc);

private:

  void hex8_shape_fcn(
    const int  &numIp,
    const double *isoParCoord,
    SharedMemView<DoubleType**> &shpfc);

  void hex8_derivative(
    const int npt, 
    const double* par_coord,
    SharedMemView<DoubleType*8*> &deriv);
};
    
} // namespace nalu
} // namespace Sierra

#endif
