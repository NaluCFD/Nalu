/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MeshDisplacementSimplifiedNeohookeanElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
#include "TimeIntegrator.h"
#include "SolutionOptions.h"

// template and scratch space
#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

template<typename AlgTraits>
MeshDisplacementSimplifiedNeohookeanElemKernel<AlgTraits>::MeshDisplacementSimplifiedNeohookeanElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  VectorFieldType* meshDisplacement,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    lrscv_(sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_)->adjacentNodes())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  meshDisplacement_ = &(meshDisplacement->field_of_state(stk::mesh::StateNP1));
  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());


  mu_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "lame_mu");
  lambda_     = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "lame_lambda");

  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);
  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);

  get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*meshDisplacement_, AlgTraits::nDim_);

  dataPreReqs.add_gathered_nodal_field(*lambda_, 1);
  dataPreReqs.add_gathered_nodal_field(*mu_, 1);

  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);

  dataPreReqs.add_master_element_call(SCS_GRAD_OP, CURRENT_COORDINATES);

}

template<typename AlgTraits>
MeshDisplacementSimplifiedNeohookeanElemKernel<AlgTraits>::~MeshDisplacementSimplifiedNeohookeanElemKernel()
{}

template<typename AlgTraits>
void
MeshDisplacementSimplifiedNeohookeanElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
}

template<typename AlgTraits>
void
MeshDisplacementSimplifiedNeohookeanElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{

  const int nDim2    = AlgTraits::nDim_*AlgTraits::nDim_;
  const int nDim     = AlgTraits::nDim_;
  const int numScsIp = AlgTraits::numScsIp_;
  const int nodesPerElement = AlgTraits::nodesPerElement_;


  NALU_ALIGNED DoubleType kd[nDim2];
  NALU_ALIGNED DoubleType gradTens[nDim2];

  // Define kd
  for ( int i = 0; i < nDim; ++i ) {
    for ( int j = 0; j < nDim; ++j ) {
      kd[nDim*i+j] = 0.0;
      if (i==j)
        kd[nDim*i+j] = 1.0;
    }
  }

  SharedMemView<DoubleType**>& meshDisp = scratchViews.get_scratch_view_2D(*meshDisplacement_);
  SharedMemView<DoubleType*>& mu = scratchViews.get_scratch_view_1D(*mu_);
  SharedMemView<DoubleType*>& lambda = scratchViews.get_scratch_view_1D(*lambda_);

  SharedMemView<DoubleType**>& v_scs_areav = scratchViews.get_me_views(CURRENT_COORDINATES).scs_areav;
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx;

  for ( int ip = 0; ip < numScsIp; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    const int ilNdim = il*nDim;
    const int irNdim = ir*nDim;

    DoubleType  muIp = 0.0;
    DoubleType lambdaIp = 0.0;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      muIp += r*mu(ic);
      lambdaIp += r*lambda(ic);
    }

    // Compute full lagged gradient
    for ( int i = 0; i < nDim2; ++i )
      gradTens[i] = 0;

    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      for ( int i = 0; i < nDim; ++i ) {
        const DoubleType dxi = meshDisp(ic,i);
        for ( int j = 0; j < nDim; ++j ) {
          gradTens[i*nDim+j] += v_dndx(ip,ic,j)*dxi;
        }
      }
    }

    for ( int i = 0; i < nDim; ++i ) {
      for ( int j = 0; j < nDim; ++j ) {
        gradTens[i*nDim+j] += kd[i*nDim+j];
      }
    }

    DoubleType Iconst = 0.0;

    if ( nDim == 3 ) {
      for ( int i = 0; i < nDim; ++i ) 
        Iconst += gradTens[i]*gradTens[i] + gradTens[i+3]*gradTens[i+3] +
          gradTens[i+6]*gradTens[i+6];
    } else {
      for ( int i = 0; i < nDim; ++i ) 
        Iconst += gradTens[i]*gradTens[i] + gradTens[i+2]*gradTens[i+2]; 
    }

    const DoubleType beta = 1.0-1.0/(1.0+Iconst);
    const DoubleType J = nDim == 3 ? (gradTens[0]*(gradTens[4]*gradTens[8]-gradTens[7]*gradTens[5]) -
               gradTens[1]*(gradTens[3]*gradTens[8]-gradTens[6]*gradTens[5]) +
               gradTens[2]*(gradTens[3]*gradTens[7]-gradTens[6]*gradTens[4])) : 
               (gradTens[0]*gradTens[3]-gradTens[2]*gradTens[1]);

    // stress
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      
      const int icNdim = ic*nDim;

      for ( int i = 0; i < nDim; ++i ) {

        const int indexL = ilNdim + i;
        const int indexR = irNdim + i;

        const DoubleType dxi = meshDisp(ic,i);

        DoubleType lhs_riC_i = 0.0; 
        for ( int j = 0; j < nDim; ++j ) {
          const DoubleType axj = v_scs_areav(ip,j);
          const DoubleType dxj = meshDisp(ic,j);

          for ( int k = 0; k < nDim; ++k ) {
            const DoubleType factor = gradTens[nDim*i+k]/J;
            const DoubleType lhsfacDiff_k = -(muIp*beta)*factor*v_dndx(ip,ic,k)*axj;
            
            lhs(indexL, icNdim+j) += lhsfacDiff_k;
            lhs(indexR, icNdim+j) -= lhsfacDiff_k;

            rhs(indexL) -= lhsfacDiff_k*dxj;
            rhs(indexR) += lhsfacDiff_k*dxj;
          }

          const DoubleType lhsfacDiff_j = -((muIp*beta)/J)*v_dndx(ip,ic,j)*axj;
          lhs_riC_i += lhsfacDiff_j;

          const DoubleType lhsfacDiff_i = -((muIp*beta)/J)*v_dndx(ip,ic,i)*axj;
          lhs(indexL,icNdim+j) += lhsfacDiff_i; 
          lhs(indexR,icNdim+j) -= lhsfacDiff_i;
          rhs(indexL) -= lhsfacDiff_i*dxj;
          rhs(indexR) += lhsfacDiff_i*dxj;
 
        }

        lhs(indexL, icNdim+i) += lhs_riC_i;
        lhs(indexR, icNdim+i) -= lhs_riC_i;
        rhs(indexL) -= lhs_riC_i*dxi;
        rhs(indexR) += lhs_riC_i*dxi;

      }
    
    }

    for ( int i = 0; i < nDim; ++i ) {
      const int indexL = ilNdim + i;
      const int indexR = irNdim + i;
      for ( int j = 0; j < nDim; ++j ) {
        const DoubleType axj = v_scs_areav(ip,j);
        rhs(indexL) += lambdaIp*kd[i*nDim+j]*axj*(J-(1.0+0.75*muIp/lambdaIp))+(muIp*beta/J)*kd[i*nDim+j]*axj;
        rhs(indexR) -= lambdaIp*kd[i*nDim+j]*axj*(J-(1.0+0.75*muIp/lambdaIp))+(muIp*beta/J)*kd[i*nDim+j]*axj;
      }
    }

  }

}

INSTANTIATE_KERNEL(MeshDisplacementSimplifiedNeohookeanElemKernel);

}  // nalu
}  // sierra
