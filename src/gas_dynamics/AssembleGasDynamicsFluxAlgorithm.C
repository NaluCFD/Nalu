/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <gas_dynamics/AssembleGasDynamicsFluxAlgorithm.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleGasDynamicsFluxAlgorithm - assembles RHS for gas dynamics; AUSM+
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleGasDynamicsFluxAlgorithm::AssembleGasDynamicsFluxAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *density,
  VectorFieldType *momentum,
  VectorFieldType *velocity,
  ScalarFieldType *totalH,
  ScalarFieldType *pressure,
  ScalarFieldType *temperature,
  ScalarFieldType *speedOfSound,
  ScalarFieldType *viscosity,
  ScalarFieldType *thermalCond,
  GenericFieldType *rhsGasDyn)
  : Algorithm(realm, part),
    density_(density),
    momentum_(momentum),
    velocity_(velocity),
    totalH_(totalH),
    pressure_(pressure),
    temperature_(temperature),
    speedOfSound_(speedOfSound),
    viscosity_(viscosity),
    thermalCond_(thermalCond),
    rhsGasDyn_(rhsGasDyn),
    velocityRTM_(NULL),
    edgeAreaVec_(NULL),
    coordinates_(NULL)
{
  // save off mising fields
  stk::mesh::MetaData & metaData = realm_.meta_data();
  if ( realm_.does_mesh_move() )
    velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity");
  edgeAreaVec_ = metaData.get_field<double>(stk::topology::EDGE_RANK, "edge_area_vector");
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleGasDynamicsFluxAlgorithm::execute()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  // sizes
  const int nDim = metaData.spatial_dimension();
  const int cOffset = nDim;
  const int eOffset = nDim + 1;

  // fixed values
  std::vector<double> ws_av(nDim);
  std::vector<double> ws_duidxj(nDim*nDim);
  double tauIp[3][3];

  // constants
  const double oneEighth = 1.0/8.0;
  const double threeSixTeenth = 3.0/16.0;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = metaData.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  //===========================================================
  // assemble edge-based flux operator to the node
  //===========================================================

  stk::mesh::BucketVector const& edge_buckets =
    realm_.get_buckets( stk::topology::EDGE_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets.begin();
        ib != edge_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to edge area vector
    double * av = stk::mesh::field_data(*edgeAreaVec_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      stk::mesh::Entity const * edge_node_rels = b.begin_nodes(k);

      // sanity check on number or nodes
      STK_ThrowAssert( b.num_nodes(k) == 2 );

      // left and right nodes
      stk::mesh::Entity nodeL = edge_node_rels[0];
      stk::mesh::Entity nodeR = edge_node_rels[1];

      // left/right nodes (all constant)
      const double * coordL = stk::mesh::field_data(*coordinates_, nodeL);
      const double * coordR = stk::mesh::field_data(*coordinates_, nodeR);

      // rho
      const double densityL = *stk::mesh::field_data( *density_, nodeL);
      const double densityR = *stk::mesh::field_data( *density_, nodeR);

      // rho*u_i
      const double * momentumL = stk::mesh::field_data( *momentum_, nodeL);
      const double * momentumR = stk::mesh::field_data( *momentum_, nodeR);

      // u_i
      const double * velocityL = stk::mesh::field_data( *velocity_, nodeL);
      const double * velocityR = stk::mesh::field_data( *velocity_, nodeR);

      // u_i - v_i
      const double * vrtmL = stk::mesh::field_data( *velocityRTM_, nodeL);
      const double * vrtmR = stk::mesh::field_data( *velocityRTM_, nodeR);

      // rhoH
      const double totalHL = *stk::mesh::field_data( *totalH_, nodeL);
      const double totalHR = *stk::mesh::field_data( *totalH_, nodeR);

      // p
      const double pressureL = *stk::mesh::field_data( *pressure_, nodeL);
      const double pressureR = *stk::mesh::field_data( *pressure_, nodeR);

      // T
      const double temperatureL = *stk::mesh::field_data( *temperature_, nodeL);
      const double temperatureR = *stk::mesh::field_data( *temperature_, nodeR);

      // c
      const double speedOfSoundL = *stk::mesh::field_data( *speedOfSound_, nodeL);
      const double speedOfSoundR = *stk::mesh::field_data( *speedOfSound_, nodeR);

      // mu
      const double viscL = *stk::mesh::field_data( *viscosity_, nodeL);
      const double viscR = *stk::mesh::field_data( *viscosity_, nodeR);

      // kappa
      const double thermalCondL = *stk::mesh::field_data( *thermalCond_, nodeL);
      const double thermalCondR = *stk::mesh::field_data( *thermalCond_, nodeR);
      
      // left/right nodes (to be assembled)
      double * rhsGasDynL = stk::mesh::field_data( *rhsGasDyn_, nodeL);
      double * rhsGasDynR = stk::mesh::field_data( *rhsGasDyn_, nodeR);

      // compute geometry; save off area vector
      double axdx = 0.0;
      double asq = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        const double axj = av[nDim*k+j];
        const double dxj = coordR[j] - coordL[j];
        asq += axj*axj;
        axdx += axj*dxj;
        ws_av[j] = axj;
      }
      const double aMag = std::sqrt(asq);
      const double inv_axdx = 1.0/axdx;

      // zero tauIp
      for ( int i = 0; i < nDim; ++i )
        for ( int j = 0; j < nDim; ++j )
          tauIp[i][j] = 0.0;

      // form duidxj now omitting over-relaxed procedure of Jasak:
      // dui/dxj = GjUi +[(uiR - uiL) - GlUi*dxl]*Aj/AxDx
      for ( int i = 0; i < nDim; ++i ) {

        // difference between R and L nodes for component i
        const double uidiff = velocityR[i] - velocityL[i];

        // offset into all forms of dudx
        const int nDimI = nDim*i;

        // start sum for NOC contribution; TBD

        // form full tensor dui/dxj without NOC
        for ( int j = 0; j < nDim; ++j ) {
          const int offSetIJ = nDimI+j;
          const double axj = ws_av[j];
          ws_duidxj[offSetIJ] = uidiff*axj*inv_axdx;
        }
      }

      // divU
      double divU = 0.0;
      for ( int j = 0; j < nDim; ++j)
        divU += ws_duidxj[j*nDim+j];

      double machNumberL = 0.0;
      double machNumberR = 0.0;
      for ( int j = 0; j < nDim; ++j) {
        const double nj = ws_av[j]/aMag;
        const double meanSpeedOfSound = 0.5*(speedOfSoundL + speedOfSoundR);
        machNumberL += vrtmL[j]*nj/meanSpeedOfSound;
        machNumberR += vrtmR[j]*nj/meanSpeedOfSound;
      }

      // AUSM quantities
      const double absML = std::abs(machNumberL);
      const double absMR = std::abs(machNumberR);
      const double signML = machNumberL > 0.0 ? 1.0 : -1.0;
      const double signMR = machNumberR > 0.0 ? 1.0 : -1.0;

      // script{M}+Left and script{M}-Right
      const double scriptMpL = ( absML >= 1.0) 
        ? 0.5*(machNumberL + absML)
        : 0.25*std::pow(machNumberL + 1.0, 2.0) 
        + oneEighth*std::pow(machNumberL*machNumberL - 1.0, 2.0);

      const double scriptMmR = ( absMR >= 1.0) 
        ? 0.5*(machNumberR - absMR)
        : -0.25*std::pow(machNumberR - 1.0, 2.0) 
        - oneEighth*std::pow(machNumberR*machNumberR - 1.0, 2.0);

      // script{P}+Left and script{P}-Right
      const double scriptPpL = ( absML >= 1.0 )
        ? 0.5*(1.0 + signML) 
        : 0.25*std::pow(machNumberL + 1.0, 2.0)*(2.0 - machNumberL) 
        + threeSixTeenth*machNumberL*std::pow(machNumberL*machNumberL - 1.0, 2.0);

      const double scriptPmR = ( absMR >= 1.0 )
        ? 0.5*(1.0 - signMR) 
        : 0.25*std::pow(machNumberR - 1.0, 2.0)*(2.0 + machNumberR) 
        - threeSixTeenth*machNumberR*std::pow(machNumberR*machNumberR - 1.0, 2.0);
      
      // left-right m's and p's
      const double mLR = scriptMpL + scriptMmR;
      const double pLR = scriptPpL*pressureL + scriptPmR*pressureR;

      // momentum first
      const double viscIp = 0.5*(viscL + viscR);
      for ( int i = 0; i < nDim; ++i ) {

        // advective
        const double fmL = momentumL[i]*speedOfSoundL;
        const double fmR = momentumR[i]*speedOfSoundR;
        const double fmLR = aMag*(0.5*mLR*(fmL + fmR) - 0.5*std::abs(mLR)*(fmR - fmL)) + pLR*ws_av[i];
        
        // diffusive
        double dfmA = 2.0/3.0*viscIp*divU*ws_av[i];
        tauIp[i][i] = -2.0/3.0*viscIp*divU;
        const int offSetI = nDim*i;
        for ( int j = 0; j < nDim; ++j ) {
          const int offSetTrans = nDim*j+i;
          const double axj = ws_av[j];
          dfmA += -viscIp*(ws_duidxj[offSetI+j] + ws_duidxj[offSetTrans])*axj;
          tauIp[i][j] += viscIp*(ws_duidxj[offSetI+j] + ws_duidxj[offSetTrans]);
        }

        // advection/diffusion assembly
        rhsGasDynL[i] -= (fmLR + dfmA);
        rhsGasDynR[i] += (fmLR + dfmA);
      }

      // continuity 
      const double fcL = densityL*speedOfSoundL;
      const double fcR = densityR*speedOfSoundR;
      const double fcLR = aMag*(0.5*mLR*(fcL + fcR) - 0.5*std::abs(mLR)*(fcR - fcL));
      // no diffusion
      rhsGasDynL[cOffset] -= fcLR;
      rhsGasDynR[cOffset] += fcLR;

      // total energy
      const double feL = totalHL*speedOfSoundL;
      const double feR = totalHR*speedOfSoundR;
      const double feLR = aMag*(0.5*mLR*(feL + feR) - 0.5*std::abs(mLR)*(feR - feL));

      // diffusion, LHS is: d/dxj(qj) - d/dxj(ui*tauij)
      const double thermalCondIp = 0.5*(thermalCondL + thermalCondR);
      double dfeA = -thermalCondIp*asq*inv_axdx*(temperatureR - temperatureL);

      for ( int i = 0; i < nDim; ++i ) {
        const double uiIp = 0.5*(velocityR[i] + velocityL[i]);
        for ( int j = 0; j < nDim; ++j ) {
          dfeA += -uiIp*tauIp[i][j]*ws_av[j];
        }
      }
      
      rhsGasDynL[eOffset] -= (feLR + dfeA);
      rhsGasDynR[eOffset] += (feLR + dfeA);
    }
  }
}

} // namespace nalu
} // namespace Sierra
