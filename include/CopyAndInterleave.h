/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef CopyAndInterleave_h
#define CopyAndInterleave_h

#include <ElemDataRequests.h>
#include <KokkosInterface.h>
#include <SimdInterface.h>
#include <ScratchViews.h>

namespace sierra{
namespace nalu{

template<typename DTYPE>
void interleave_3D(SharedMemView<DTYPE***>& dview, const SharedMemView<double***>& sview, int simdIndex, int simdLen)
{
    int dim0 = dview.dimension(0);
    int dim1 = dview.dimension(1);
    int dim2 = dview.dimension(2);

//    int len = dim0*dim1*dim2;
//    DTYPE& dval = dview(0,0,0);
//    double* dptr = &dval[simdIndex];
//    double* sptr = &sview(0,0,0);
//    for(int i=0; i<len; ++i) {
//        *dptr = *sptr++;
//        dptr += simdLen;
//    }
    for(int i=0; i<dim0; ++i) {
        for(int j=0; j<dim1; ++j) {
            for(int k=0; k<dim2; ++k) {
                stk::simd::set_data(dview(i,j,k), simdIndex, sview(i,j,k));
            }   
        }   
    }   
}

template<typename DTYPE>
void interleave_1D(SharedMemView<DTYPE*>& dview, const SharedMemView<double*>& sview, int simdIndex, int simdLen)
{
    int dim = dview.dimension(0);
//    DTYPE& dval = dview(0);
//    double* dptr = &dval[simdIndex];
//    double* sptr = &sview(0);

    for(int i=0; i<dim; ++i) {
        stk::simd::set_data(dview(i), simdIndex, sview(i));
    }
//    for(int i=0; i<dim; ++i) {
//       *dptr = *sptr++;
//       dptr += simdLen;
//    }   
}

template<typename DTYPE>
void interleave_2D(SharedMemView<DTYPE**>& dview, const SharedMemView<double**>& sview, int simdIndex)
{
    int dim0 = dview.dimension(0);
    int dim1 = dview.dimension(1);
    for(int i=0; i<dim0; ++i) {
        for(int j=0; j<dim1; ++j) {
            stk::simd::set_data(dview(i,j), simdIndex, sview(i,j));
        }
    }   
}

inline
void interleave_1D(ViewHolder* dest, const ViewHolder* src, int simdIndex, int simdLen)
{
    const SharedMemView<double*>& src_view = static_cast<const ViewT<SharedMemView<double*>>*>(src)->view_;
    SharedMemView<DoubleType*>& dest_view = static_cast<ViewT<SharedMemView<DoubleType*>>*>(dest)->view_;
    interleave_1D(dest_view, src_view, simdIndex, simdLen);
}

inline
void interleave_2D(ViewHolder* dest, const ViewHolder* src, int simdIndex)
{
    const SharedMemView<double**>& src_view = static_cast<const ViewT<SharedMemView<double**>>*>(src)->view_;
    SharedMemView<DoubleType**>& dest_view = static_cast<ViewT<SharedMemView<DoubleType**>>*>(dest)->view_;
    interleave_2D(dest_view, src_view, simdIndex);
}

inline
void interleave_3D(ViewHolder* dest, const ViewHolder* src, int simdIndex, int simdLen)
{
    const SharedMemView<double***>& src_view = static_cast<const ViewT<SharedMemView<double***>>*>(src)->view_;
    SharedMemView<DoubleType***>& dest_view = static_cast<ViewT<SharedMemView<DoubleType***>>*>(dest)->view_;
    interleave_3D(dest_view, src_view, simdIndex, simdLen);
}

inline
void interleave_me_views(MasterElementViews<DoubleType>& dest,
                         const MasterElementViews<double>& src,
                         int simdIndex, int simdLen)
{
  interleave_2D(dest.scs_areav, src.scs_areav, simdIndex);
  interleave_3D(dest.dndx, src.dndx, simdIndex, simdLen);
  interleave_3D(dest.dndx_shifted, src.dndx_shifted, simdIndex, simdLen);
  interleave_3D(dest.dndx_fem, src.dndx_fem, simdIndex, simdLen);
  interleave_1D(dest.deriv, src.deriv, simdIndex, simdLen);
  interleave_1D(dest.deriv_fem, src.deriv_fem, simdIndex, simdLen);
  interleave_1D(dest.det_j, src.det_j, simdIndex, simdLen);
  interleave_1D(dest.det_j_fem, src.det_j_fem, simdIndex, simdLen);
  interleave_1D(dest.scv_volume, src.scv_volume, simdIndex, simdLen);
  interleave_3D(dest.gijUpper, src.gijUpper, simdIndex, simdLen);
  interleave_3D(dest.gijLower, src.gijLower, simdIndex, simdLen);
}

inline
void copy_and_interleave(const std::vector<ScratchViews<double>*>& data,
                         int simdElems,
                         int simdLen,
                         ScratchViews<DoubleType>& simdData,
                         bool copyMEViews = true)
{
    const std::vector<ViewHolder*>& simdFieldViews = simdData.get_field_views();

    for(int simdIndex=0; simdIndex<simdElems; ++simdIndex) {
      const std::vector<ViewHolder*>& fieldViews = data[simdIndex]->get_field_views();
      for(size_t i=0; i<fieldViews.size(); ++i) {
        if (fieldViews[i] != nullptr) {
          switch(fieldViews[i]->dim_) {
          case 1: interleave_1D(simdFieldViews[i], fieldViews[i], simdIndex, simdLen); break;
          case 2: interleave_2D(simdFieldViews[i], fieldViews[i], simdIndex); break;
          case 3: interleave_3D(simdFieldViews[i], fieldViews[i], simdIndex, simdLen); break;
          default: ThrowRequireMsg(fieldViews[i]->dim_ > 0 && fieldViews[i]->dim_ < 4, "ERROR, view dim out of range: "<<fieldViews[i]->dim_);
             break;
          }
        }
      }
    }

    if (copyMEViews)
    {
      for(int simdIndex=0; simdIndex<simdElems; ++simdIndex) {
        if (simdData.has_coord_field(CURRENT_COORDINATES)) {
          interleave_me_views(simdData.get_me_views(CURRENT_COORDINATES), data[simdIndex]->get_me_views(CURRENT_COORDINATES), simdIndex, simdLen);
        }
        if (simdData.has_coord_field(MODEL_COORDINATES)) {
          interleave_me_views(simdData.get_me_views(MODEL_COORDINATES), data[simdIndex]->get_me_views(MODEL_COORDINATES), simdIndex, simdLen);
        }
      }
    }
}

inline
void extract_vector_lane(const SharedMemView<DoubleType*>& simdrhs, int simdIndex, SharedMemView<double*>& rhs)
{
  for(size_t i=0; i<simdrhs.dimension(0); ++i) {
    rhs(i) = stk::simd::get_data(simdrhs(i), simdIndex);
  }
}

inline
void extract_vector_lane(const SharedMemView<DoubleType**>& simdlhs, int simdIndex, SharedMemView<double**>& lhs)
{
  for(size_t i=0; i<simdlhs.dimension(0); ++i) {
    for(size_t j=0; j<simdlhs.dimension(1); ++j) {
      lhs(i,j) = stk::simd::get_data(simdlhs(i,j), simdIndex);
    }
  }
}

} // namespace nalu
} // namespace Sierra

#endif
