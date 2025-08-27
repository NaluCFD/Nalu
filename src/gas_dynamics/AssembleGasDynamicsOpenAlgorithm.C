/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <gas_dynamics/AssembleGasDynamicsOpenAlgorithm.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
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
// AssembleGasDynamicsOpenAlgorithm - assembles RHS for gas dynamics; open
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleGasDynamicsOpenAlgorithm::AssembleGasDynamicsOpenAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *density,
  VectorFieldType *momentum,
  ScalarFieldType *totalH,
  ScalarFieldType *bcPressure,
  ScalarFieldType *bcTemperature,
  ScalarFieldType *cp,
  ScalarFieldType *gamma,
  GenericFieldType *rhsGasDyn)
  : Algorithm(realm, part),
    density_(density),
    momentum_(momentum),
    totalH_(totalH),
    bcPressure_(bcPressure),
    bcTemperature_(bcTemperature),
    cp_(cp),
    gamma_(gamma),
    rhsGasDyn_(rhsGasDyn),
    velocityRTM_(NULL),
    exposedAreaVec_(NULL)
{
  // save off mising fields
  stk::mesh::MetaData & metaData = realm_.meta_data();
  if ( realm_.does_mesh_move() )
    velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity");
  exposedAreaVec_ = metaData.get_field<double>(metaData.side_rank(), "exposed_area_vector");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleGasDynamicsOpenAlgorithm::execute()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  // sizes
  const int nDim = metaData.spatial_dimension();
  const int cOffset = nDim;
  const int eOffset = nDim + 1;

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

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( metaData.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());

    // extract master element specifics
    const int numScsIp = meFC->numIntPoints_;
    const int *ipNodeMap = meFC->ipNodeMap();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, b, k);

      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);

      // start assembly
      for ( int ip = 0; ip < numScsIp; ++ip ) {

        // compute aMag
        double aMag = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          aMag += areaVec[ip*nDim+i]*areaVec[ip*nDim+i];
        }
        aMag = std::sqrt(aMag);

        // nearest node
        const int nn = ipNodeMap[ip];
        
        stk::mesh::Entity nodeNN = face_node_rels[nn];

        // rho
        const double density = *stk::mesh::field_data(*density_, nodeNN);

        // rho*u_i
        const double * momentum = stk::mesh::field_data(*momentum_, nodeNN);

        // u_i - v_i
        const double * velocityRTM = stk::mesh::field_data(*velocityRTM_, nodeNN);

        // rhoH
        const double totalH = *stk::mesh::field_data( *totalH_, nodeNN);

        // p
        const double bcPressure = *stk::mesh::field_data(*bcPressure_, nodeNN);

        // T
        const double bcTemperature = *stk::mesh::field_data(*bcTemperature_, nodeNN);

        // Cp
        const double cp = *stk::mesh::field_data(*cp_, nodeNN);

        // gamma
        const double gamma = *stk::mesh::field_data(*gamma_, nodeNN);

        // pointer to fields to assemble
        double *rhsGasDynNN = stk::mesh::field_data(*rhsGasDyn_, nodeNN);

        // compute uSq and a check for the flow being into or out of the domain
        double uSq = 0.0;
        double inOutCheck = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double nj = areaVec[ip*nDim+j]/aMag;
          inOutCheck += velocityRTM[j]*nj; 
          uSq += velocityRTM[j];
        }
               
        if ( inOutCheck > 0.0 ) {
          // leaving; compute c based on specified static pressure
          const double speedOfSound = std::sqrt(gamma*bcPressure/density);

          double machNumber = 0.0;
          for ( int j = 0; j < nDim; ++j) {
            const double nj = areaVec[ip*nDim+j]/aMag;
            machNumber += velocityRTM[j]*nj/speedOfSound;
          }
          const double absM = std::abs(machNumber);
          const double signM = machNumber > 0.0 ? 1.0 : -1.0;
          
          const double scriptM = ( absM >= 1.0) 
            ? 0.5*(machNumber + absM)
            : 0.25*std::pow(machNumber + 1.0, 2.0) 
            + oneEighth*std::pow(machNumber*machNumber - 1.0, 2.0);
          
          const double scriptP = ( absM >= 1.0 )
            ? 0.5*(1.0 + signM) 
            : 0.25*std::pow(machNumber + 1.0, 2.0)*(2.0 - machNumber) 
            + threeSixTeenth*machNumber*std::pow(machNumber*machNumber - 1.0, 2.0);
          
          for ( int i = 0; i < nDim; ++i ) {
            rhsGasDynNN[i] -= aMag*scriptM*momentum[i]*speedOfSound + scriptP*bcPressure*areaVec[ip*nDim+i];
          }
          
          // continuity
          rhsGasDynNN[cOffset] -= aMag*scriptM*density*speedOfSound;
          
          // energy
          rhsGasDynNN[eOffset] -= aMag*scriptM*totalH*speedOfSound;
        }
        else {
          // entering; compute c based on specified total pressure: pTotal = pStatic + 1/2*rho*u*u
          const double speedOfSound = std::sqrt(gamma*(bcPressure - 0.5*density*uSq)/density);

          double machNumber = 0.0;
          for ( int j = 0; j < nDim; ++j) {
            const double nj = areaVec[ip*nDim+j]/aMag;
            machNumber += velocityRTM[j]*nj/speedOfSound;
          }
          const double absM = std::abs(machNumber);
          const double signM = machNumber > 0.0 ? 1.0 : -1.0;

          // interprete specified temperature as total temperature; Ttot = (1+0.5*Ma*Ma)*T
          const double correctedTemperature = bcTemperature/(1.0 + gamma*machNumber*machNumber);
          const double staticH = cp*correctedTemperature; // reference temperature is assumed to be zero
          const double totalHfromT = density*(staticH + 0.5*uSq);

          const double scriptM = ( absM > 1.0) 
            ? 0.5*(machNumber - absM)
            : -0.25*std::pow(machNumber - 1.0, 2.0) 
            - oneEighth*std::pow(machNumber*machNumber - 1.0, 2.0);
          
          const double scriptP = ( absM >= 1.0 )
            ? 0.5*(1.0 - signM) 
            : 0.25*std::pow(machNumber - 1.0, 2.0)*(2.0 + machNumber) 
            - threeSixTeenth*machNumber*std::pow(machNumber*machNumber - 1.0, 2.0);

          for ( int i = 0; i < nDim; ++i ) {
            rhsGasDynNN[i] -= aMag*scriptM*momentum[i]*speedOfSound 
              + scriptP*(bcPressure - 0.5*density*uSq)*areaVec[ip*nDim+i];
          }
          
          // continuity
          rhsGasDynNN[cOffset] -= aMag*scriptM*density*speedOfSound;
          
          // energy
          rhsGasDynNN[eOffset] -= aMag*scriptM*totalHfromT*speedOfSound;
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
