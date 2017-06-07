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
void interleave_3D(SharedMemView<DTYPE***>& dview, const SharedMemView<double***> sview, int simdIndex)
{
    int dim0 = dview.dimension(0);
    int dim1 = dview.dimension(1);
    int dim2 = dview.dimension(2);

    for(int i=0; i<dim0; ++i) {
        for(int j=0; j<dim1; ++j) {
            for(int k=0; k<dim2; ++k) {
                stk::simd::set_data(dview(i,j,k), simdIndex, sview(i,j,k));
            }   
        }   
    }   
}

template<typename DTYPE>
void interleave_1D(SharedMemView<DTYPE*>& dview, const SharedMemView<double*>& sview, int simdIndex)
{
    int dim = dview.dimension(0);
    for(int i=0; i<dim; ++i) {
       stk::simd::set_data(dview(i), simdIndex, sview(i));
    }   
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
void interleave_1D(ViewHolder* dest, const ViewHolder* src, int simdIndex)
{
    const SharedMemView<double*>& src_view = static_cast<const ViewT<SharedMemView<double*>>*>(src)->view_;
    SharedMemView<DoubleType*>& dest_view = static_cast<ViewT<SharedMemView<DoubleType*>>*>(dest)->view_;
    interleave_1D(dest_view, src_view, simdIndex);
}

inline
void interleave_2D(ViewHolder* dest, const ViewHolder* src, int simdIndex)
{
    const SharedMemView<double**>& src_view = static_cast<const ViewT<SharedMemView<double**>>*>(src)->view_;
    SharedMemView<DoubleType**>& dest_view = static_cast<ViewT<SharedMemView<DoubleType**>>*>(dest)->view_;
    interleave_2D(dest_view, src_view, simdIndex);
}

inline
void interleave_3D(ViewHolder* dest, const ViewHolder* src, int simdIndex)
{
    const SharedMemView<double***>& src_view = static_cast<const ViewT<SharedMemView<double***>>*>(src)->view_;
    SharedMemView<DoubleType***>& dest_view = static_cast<ViewT<SharedMemView<DoubleType***>>*>(dest)->view_;
    interleave_3D(dest_view, src_view, simdIndex);
}

inline
void interleave_me_views(MasterElementViews<DoubleType>& dest,
                         const MasterElementViews<double>& src,
                         int simdIndex)
{
  interleave_2D(dest.scs_areav, src.scs_areav, simdIndex);
  interleave_3D(dest.dndx, src.dndx, simdIndex);
  interleave_3D(dest.dndx_shifted, src.dndx_shifted, simdIndex);
  interleave_3D(dest.dndx_fem, src.dndx_fem, simdIndex);
  interleave_1D(dest.deriv, src.deriv, simdIndex);
  interleave_1D(dest.deriv_fem, src.deriv_fem, simdIndex);
  interleave_1D(dest.det_j, src.det_j, simdIndex);
  interleave_1D(dest.det_j_fem, src.det_j_fem, simdIndex);
  interleave_1D(dest.scv_volume, src.scv_volume, simdIndex);
  interleave_3D(dest.gijUpper, src.gijUpper, simdIndex);
  interleave_3D(dest.gijLower, src.gijLower, simdIndex);
}

inline
void copy_and_interleave(const std::vector<ScratchViews<double>*>& data,
                         int simdLen,
                         ScratchViews<DoubleType>& simdData)
{
    const std::vector<ViewHolder*>& simdFieldViews = simdData.get_field_views();

    for(int simdIndex=0; simdIndex<simdLen; ++simdIndex) {
      const std::vector<ViewHolder*>& fieldViews = data[simdIndex]->get_field_views();
      for(size_t i=0; i<fieldViews.size(); ++i) {
        if (fieldViews[i] != nullptr) {
          switch(fieldViews[i]->dim_) {
          case 1: interleave_1D(simdFieldViews[i], fieldViews[i], simdIndex); break;
          case 2: interleave_2D(simdFieldViews[i], fieldViews[i], simdIndex); break;
          case 3: interleave_3D(simdFieldViews[i], fieldViews[i], simdIndex); break;
          default: ThrowRequireMsg(fieldViews[i]->dim_ > 0 && fieldViews[i]->dim_ < 4, "ERROR, view dim out of range: "<<fieldViews[i]->dim_);
             break;
          }
        }
      }
    }

    for(int simdIndex=0; simdIndex<simdLen; ++simdIndex) {
      if (simdData.has_coord_field(CURRENT_COORDINATES)) {
        interleave_me_views(simdData.get_me_views(CURRENT_COORDINATES), data[simdIndex]->get_me_views(CURRENT_COORDINATES), simdIndex);
      }
      if (simdData.has_coord_field(MODEL_COORDINATES)) {
        interleave_me_views(simdData.get_me_views(MODEL_COORDINATES), data[simdIndex]->get_me_views(MODEL_COORDINATES), simdIndex);
      }
    }
}

} // namespace nalu
} // namespace Sierra

#endif
