/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SCALARELEMDIFFUSIONFUNCTOR_H_
#define SCALARELEMDIFFUSIONFUNCTOR_H_

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <SupplementalAlgorithm.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{


void metricInverse2D(const double * dx_dr, double * dr_dx)
{
  const double Jxr = dx_dr[0]*dx_dr[3]-dx_dr[1]*dx_dr[2];
  const double IJxr = 1/Jxr;

  dr_dx[0] = IJxr * dx_dr[3];
  dr_dx[1] = - IJxr * dx_dr[1];

  dr_dx[2] = - IJxr * dx_dr[2];
  dr_dx[3] = IJxr * dx_dr[0];

}

void metricInverse3D(const double * dx_dr, double * dr_dx)
{
  const double Jxr = dx_dr[0] * (dx_dr[4]*dx_dr[8] - dx_dr[5]*dx_dr[7])
    - dx_dr[1] * (dx_dr[3]*dx_dr[8] - dx_dr[5]*dx_dr[6])
    + dx_dr[2] * (dx_dr[3]*dx_dr[7] - dx_dr[4]*dx_dr[6]);
  const double IJxr = 1/Jxr;

  dr_dx[0] = IJxr * (dx_dr[4]*dx_dr[8] - dx_dr[5]*dx_dr[7]);
  dr_dx[1] = IJxr * (dx_dr[2]*dx_dr[7] - dx_dr[1]*dx_dr[8]);
  dr_dx[2] = IJxr * (dx_dr[1]*dx_dr[5] - dx_dr[2]*dx_dr[4]);

  dr_dx[3] = IJxr * (dx_dr[5]*dx_dr[6] - dx_dr[3]*dx_dr[8]);
  dr_dx[4] = IJxr * (dx_dr[0]*dx_dr[8] - dx_dr[2]*dx_dr[6]);
  dr_dx[5] = IJxr * (dx_dr[2]*dx_dr[3] - dx_dr[0]*dx_dr[5]);

  dr_dx[6] = IJxr * (dx_dr[3]*dx_dr[7] - dx_dr[4]*dx_dr[6]);
  dr_dx[7] = IJxr * (dx_dr[1]*dx_dr[6] - dx_dr[0]*dx_dr[7]);
  dr_dx[8] = IJxr * (dx_dr[0]*dx_dr[4] - dx_dr[1]*dx_dr[3]);

}

class ScalarElemDiffusionFunctor{

protected:
  //Bucket and Element Data
  stk::mesh::Bucket * b_;
  stk::mesh::BulkData & bulk_data_;
  stk::mesh::MetaData & meta_data_;
  MasterElement * meSCS_;

  //InputFields
  ScalarFieldType & scalarQ_;
  ScalarFieldType & diffFluxCoeff_;
  VectorFieldType & coordinates_;

  //Output
  double * p_lhs_;
  double * p_rhs_;
  std::vector<stk::mesh::Entity> * p_connected_nodes_;

  //Parameters
  const int *lrscv;
  const int nDim_;
  int numScsIp_;
  int nodesPerElement_;

public:
  ScalarElemDiffusionFunctor(
      stk::mesh::BulkData & bulk_data,
      stk::mesh::MetaData & meta_data,
      ScalarFieldType & scalarQ,
      ScalarFieldType & diffFluxCoeff,
      VectorFieldType & coordinates,
      int nDim):
      b_(0),
      bulk_data_(bulk_data),
      meta_data_(meta_data),
      meSCS_(0),
      scalarQ_(scalarQ),
      diffFluxCoeff_(diffFluxCoeff),
      coordinates_(coordinates),
      p_lhs_(0),
      p_rhs_(0),
      p_connected_nodes_(0),
      lrscv(0),
      nDim_(nDim),
      numScsIp_(0),
      nodesPerElement_(0)
  {
  }

  virtual void bind_data(stk::mesh::Bucket & b,
    MasterElement & meSCS,
    double * p_lhs,
    double * p_rhs,
    std::vector<stk::mesh::Entity> & connected_nodes)
  {
    b_ = &b;
    meSCS_ = &meSCS;
    p_lhs_ = p_lhs;
    p_rhs_ = p_rhs;
    p_connected_nodes_ = &connected_nodes;
    lrscv = meSCS_->adjacentNodes();
    nodesPerElement_ = meSCS_->nodesPerElement_;
    numScsIp_ = meSCS_->numIntPoints_;
  }

  virtual void release_data()
  {
    // does nothing
  }


  virtual ~ScalarElemDiffusionFunctor(){}

  virtual void operator()(int elem_offset) = 0;

};

class CVFEMScalarElemDiffusionFunctor: public ScalarElemDiffusionFunctor{

  double * p_shape_function_;

public:
  CVFEMScalarElemDiffusionFunctor(
      stk::mesh::BulkData & bulk_data,
      stk::mesh::MetaData & meta_data,
      ScalarFieldType & scalarQ,
      ScalarFieldType & diffFluxCoeff,
      VectorFieldType & coordinates,
      int nDim):
      ScalarElemDiffusionFunctor(bulk_data, meta_data,
        scalarQ, diffFluxCoeff, coordinates, nDim),
      p_shape_function_(0)
  { }

  virtual ~CVFEMScalarElemDiffusionFunctor(){}

  virtual void bind_data(stk::mesh::Bucket & b,
    MasterElement & meSCS,
    double * p_lhs,
    double * p_rhs,
    std::vector<stk::mesh::Entity> & connected_nodes)
  {
    ScalarElemDiffusionFunctor::bind_data(b, meSCS, p_lhs, p_rhs, connected_nodes);

    p_shape_function_ = new double[numScsIp_*nodesPerElement_];
    meSCS_->shape_fcn(&p_shape_function_[0]);
  }

  virtual void release_data()
  {
    delete[] p_shape_function_;
  }

  void operator()(int elem_offset){
    // get elem
    const stk::mesh::Entity elem = (*b_)[elem_offset];

    // get nodes
    stk::mesh::Entity const * node_rels = bulk_data_.begin_nodes(elem);

    // temporary arrays
    double p_scalarQ[nodesPerElement_];
    double p_diffFluxCoeff[nodesPerElement_];
    double p_coordinates[nodesPerElement_*nDim_];
    double p_scs_areav[numScsIp_*nDim_];
    double p_dndx[nDim_*numScsIp_*nodesPerElement_];
    double p_deriv[nDim_*numScsIp_*nodesPerElement_];
    double p_det_j[numScsIp_];

    const int lhsSize = nodesPerElement_*nodesPerElement_;
    const int rhsSize = nodesPerElement_;
    // zero lhs/rhs
    for ( int p = 0; p < lhsSize; ++p )
      p_lhs_[p] = 0.0;
    for ( int p = 0; p < rhsSize; ++p )
      p_rhs_[p] = 0.0;

    for ( int ni = 0; ni < nodesPerElement_; ++ni ) {
      stk::mesh::Entity node = node_rels[ni];

      // set connected nodes
      (*p_connected_nodes_)[ni] = node;

      const double * coords = stk::mesh::field_data(coordinates_, node );

      // gather scalars
      p_scalarQ[ni] = *stk::mesh::field_data(scalarQ_, node );
      p_diffFluxCoeff[ni] = *stk::mesh::field_data(diffFluxCoeff_, node);

      // gather vectors
      const int offSet = ni*nDim_;
      for ( int j=0; j < nDim_; ++j ) {
        p_coordinates[offSet+j] = coords[j];
      }
    }

    // compute geometry
    double scs_error = 0.0;
    meSCS_->determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);
    // compute dndx
    meSCS_->grad_op(1, &p_coordinates[0], &p_dndx[0], &p_deriv[0], &p_det_j[0], &scs_error);

    // start assembly
    for ( int ip = 0; ip < numScsIp_; ++ip ) {

      // left and right nodes for this ip
      const int il = lrscv[2*ip];
      const int ir = lrscv[2*ip+1];

      // corresponding matrix rows
      const int rowL = il*nodesPerElement_;
      const int rowR = ir*nodesPerElement_;

      // save off ip values; offset to Shape Function
      double muIp = 0.0;
      const int offSetSF = ip*nodesPerElement_;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {
        const double r = p_shape_function_[offSetSF+ic];
        muIp += r*p_diffFluxCoeff[ic];
      }

      double qDiff = 0.0;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {

        // diffusion
        double lhsfacDiff = 0.0;
        const int offSetDnDx = nDim_*nodesPerElement_*ip + ic*nDim_;
        for ( int j = 0; j < nDim_; ++j ) {
          lhsfacDiff += -muIp*p_dndx[offSetDnDx+j]*p_scs_areav[ip*nDim_+j];
        }

        qDiff += lhsfacDiff*p_scalarQ[ic];

        // lhs; il then ir
        p_lhs_[rowL+ic] += lhsfacDiff;
        p_lhs_[rowR+ic] -= lhsfacDiff;
      }

      // rhs; il then ir
      p_rhs_[il] -= qDiff;
      p_rhs_[ir] += qDiff;

    }

  }

};

class CollocationScalarElemDiffusionFunctor: public ScalarElemDiffusionFunctor{

  double * p_deriv_;
public:
  CollocationScalarElemDiffusionFunctor(
      stk::mesh::BulkData & bulk_data,
      stk::mesh::MetaData & meta_data,
      ScalarFieldType & scalarQ,
      ScalarFieldType & diffFluxCoeff,
      VectorFieldType & coordinates,
      int nDim):
      ScalarElemDiffusionFunctor(bulk_data, meta_data,
        scalarQ, diffFluxCoeff, coordinates, nDim),
      p_deriv_(0)
  {
  }

  virtual ~CollocationScalarElemDiffusionFunctor()
  {
  }

  virtual void bind_data(stk::mesh::Bucket & b,
    MasterElement & meSCS,
    double * p_lhs,
    double * p_rhs,
    std::vector<stk::mesh::Entity> & connected_nodes)
  {
    ScalarElemDiffusionFunctor::bind_data(b, meSCS, p_lhs, p_rhs, connected_nodes);

    p_deriv_ = new double[nDim_*nodesPerElement_*nodesPerElement_];
    double scs_error = 0.0;
    // compute dndx
    meSCS_->nodal_grad_op(1, p_deriv_, &scs_error);
  }

  virtual void release_data()
  {
    delete p_deriv_;
    p_deriv_ = 0;
  }

  void operator()(int elem_offset){
    // get elem
    const stk::mesh::Entity elem = (*b_)[elem_offset];

    // get nodes
    stk::mesh::Entity const * node_rels = bulk_data_.begin_nodes(elem);

    // temporary arrays
    double p_scalarQ[nodesPerElement_];
    double p_diffFluxCoeff[nodesPerElement_];
    double p_coordinates[nodesPerElement_*nDim_];
    double p_scs_areav[numScsIp_*nDim_];

    const int lhsSize = nodesPerElement_*nodesPerElement_;
    const int rhsSize = nodesPerElement_;
    // zero lhs/rhs
    for ( int p = 0; p < lhsSize; ++p )
      p_lhs_[p] = 0.0;
    for ( int p = 0; p < rhsSize; ++p )
      p_rhs_[p] = 0.0;

    for ( int ni = 0; ni < nodesPerElement_; ++ni ) {
      stk::mesh::Entity node = node_rels[ni];

      // set connected nodes
      (*p_connected_nodes_)[ni] = node;

      const double * coords = stk::mesh::field_data(coordinates_, node );

      // gather scalars
      p_scalarQ[ni] = *stk::mesh::field_data(scalarQ_, node );
      p_diffFluxCoeff[ni] = *stk::mesh::field_data(diffFluxCoeff_, node);

      // gather vectors
      const int offSet = ni*nDim_;
      for ( int j=0; j < nDim_; ++j ) {
        p_coordinates[offSet+j] = coords[j];
      }
    }

    // compute geometry
    double scs_error = 0.0;
    meSCS_->determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

    double x_r[nodesPerElement_*nDim_*nDim_];
    double r_x[nodesPerElement_*nDim_*nDim_];
    for (int i = 0; i < nodesPerElement_*nDim_*nDim_; ++i)
    {
      x_r[i] = 0;
      r_x[i] = 0;
    }

    // compute gradients
    for (int i = 0; i < nodesPerElement_; ++i)
    {
      double * p_x_r = &x_r[i*nDim_*nDim_];

      for (int j = 0; j < nodesPerElement_; ++j)
      {
        const int ioffset = i*nDim_*nodesPerElement_ + j*nDim_;

        for (int gradDir = 0; gradDir < nDim_; ++gradDir)
          for (int coordDir = 0; coordDir < nDim_; ++coordDir)
            p_x_r[gradDir*nDim_+coordDir] += p_deriv_[ioffset+gradDir]*p_coordinates[j*nDim_+coordDir];
      }

      double * p_r_x = &r_x[i*nDim_*nDim_];

      if (nDim_ == 3){
        metricInverse3D(p_x_r, p_r_x);
      }
      else{
        metricInverse2D(p_x_r, p_r_x);
      }

    }

    // start assembly
    for ( int ip = 0; ip < numScsIp_; ++ip ) {

      // left and right nodes for this ip
      const int il = lrscv[2*ip];
      const int ir = lrscv[2*ip+1];

      // corresponding matrix rows
      const int rowL = il*nodesPerElement_;
      const int rowR = ir*nodesPerElement_;

      // diffusion
      double qDiff = 0.0;
      for ( int edgeNode = 0; edgeNode < 2; ++edgeNode ) {
        const int ic = lrscv[2*ip+edgeNode];
        const double r = 0.5*p_diffFluxCoeff[ic];
        for (int iDir = 0; iDir < nDim_; ++iDir)
        {
          const double rfac = r*p_scs_areav[ip*nDim_+iDir];
          const double * p_r_x = &r_x[ic*nDim_*nDim_+iDir*nDim_];
          const int ioffset = ic*nDim_*nodesPerElement_;
          for (int jc = 0; jc < nodesPerElement_; ++jc){
            const int joffset = ioffset+jc*nDim_;
            double lhsfacDiff = 0.0;
            for (int jDir = 0; jDir < nDim_; ++jDir)
              lhsfacDiff -= rfac*p_deriv_[joffset+jDir]*p_r_x[jDir];
            p_lhs_[rowL+jc] += lhsfacDiff;
            p_lhs_[rowR+jc] -= lhsfacDiff;

            qDiff += lhsfacDiff*p_scalarQ[jc];

          }

        }
      }

      // rhs; il then ir
      p_rhs_[il] -= qDiff;
      p_rhs_[ir] += qDiff;

    }

  }

};

}

}

#endif /* SCALARELEMDIFFUSIONFUNCTOR_H_ */
