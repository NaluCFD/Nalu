/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Algorithm.h>
#include <TablePropAlgorithm.h>
#include <FieldTypeDef.h>
#include <PropertyEvaluator.h>
#include <Realm.h>
#include <tabular_props/StateTable.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <vector>

namespace sierra{
namespace nalu{

TablePropAlgorithm::TablePropAlgorithm(
  Realm & realm,
  stk::mesh::Part * part,
  StateTable *stateTable,
  stk::mesh::FieldBase * prop,
  std::string tablePropName,
  std::vector<std::string> &indVarNameVec,
  std::vector<std::string> &indVarTableNameVec,
  const size_t cIndex,
  const stk::mesh::MetaData &meta_data)
  : Algorithm(realm, part),
    prop_(prop),
    tablePropName_(tablePropName),
    indVarTableNameVec_(indVarTableNameVec),
    cIndex_(cIndex),
    indVarSize_(indVarNameVec.size()),
    iMin_(0),
    iMax_(0)
{
  // extract the table; save off iMix, iMax
  table_ = *stateTable->get_entry();
  iMin_ = 0;
  iMax_ = table_.size() -1;

  // extract the independent fields; check if there is one..
  if ( indVarSize_ == 0 )
    throw std::runtime_error("TablePropAlgorithm: independent variable size is zero:");

  indVar_.resize(indVarSize_);
  for ( size_t k = 0; k < indVarSize_; ++k) {
    ScalarFieldType *indVar = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, indVarNameVec[k]);
    if ( NULL == indVar ) {
      throw std::runtime_error("TablePropAlgorithm: independent variable not registered:");
    }
    else {
      indVar_[k] = indVar;
    }
  }

  // resize some work vectors
  workIndVar_.resize(indVarSize_);
  workZ_.resize(indVarSize_);

  // provide some output
  Env::outputP0() << "the Following Table Property name will be extracted: " << tablePropName << std::endl;
  for ( size_t k = 0; k < indVarTableNameVec_.size(); ++k ) {
    Env::outputP0() << "using independent variables: " << indVarTableNameVec_[k] << std::endl;
  }

}

void
TablePropAlgorithm::execute()
{

  // make sure that partVec_ is size one
  ThrowAssert( partVec_.size() == 1 );
  
  stk::mesh::Selector selector = stk::mesh::selectUnion(partVec_);
  
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, selector );

  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // property is always single in size
    double *prop  = (double*) stk::mesh::field_data(*prop_, b);

    // independent variable size can be more than one
    for ( size_t l = 0; l < indVarSize_; ++l) {
      double *indVar  = (double*) stk::mesh::field_data(*indVar_[l], b);
      workIndVar_[l] = indVar;
    }

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // load up the independent work array
      for ( size_t l = 0; l < indVarSize_; ++l) {
        double *z = workIndVar_[l];
        workZ_[l] = z[k];
      }
      prop[k] = interp_property(&workZ_[0]);
    }
  }
}

size_t
TablePropAlgorithm::find_index( 
  const double z,
  size_t iMin,
  size_t iMax)
{  
  // mid point of table
  size_t iMid = iMin + (iMax - iMin )/2;

  // find the mid point balue for z; call it the candidate
  std::vector<double> &tableK = table_[iMid];
  const double candidateZ = tableK[0];
  
  if ( iMid == iMin) {
    return iMid;
  }

  if ( candidateZ > z ) {
    // in the lower
    return find_index(z, iMin, iMid-1);
  }
  else {
    // in the upper
    return find_index(z, iMid+1, iMax);
  }

} 

double
TablePropAlgorithm::interp_property( 
  const double *z)
{
  // find the index
  size_t foundIndex = find_index(z[0], iMin_, iMax_);

  // proceed with the simple linear interpolation
  double interpProp = 0.0;

  // check for high/low bounds
  if ( foundIndex == iMin_ || foundIndex == iMax_)
    interpProp = table_[foundIndex][cIndex_];
  else {
    // not on the outer edge of the table
    if ( table_[foundIndex][0] > z[0] ) {
      // on the high side
      size_t foundM1 = foundIndex-1;
      
      // values of mixFrac
      const double x1 = table_[foundIndex][0];
      const double x0 = table_[foundM1][0];
      
      // property
      double y1 = table_[foundIndex][cIndex_];
      double y0 = table_[foundM1][cIndex_];
      interpProp = y0 + (y1 - y0 )*(z[0]-x0)/(x1-x0);
    }
    else {
      // on the low side
      size_t foundP1 = foundIndex+1;
      
      // values of mixFrac
      const double x0 = table_[foundIndex][0];
      const double x1 = table_[foundP1][0];
      
      // property
      double y0 = table_[foundIndex][cIndex_];
      double y1 = table_[foundP1][cIndex_];
      
      interpProp = y0 + (y1 - y0 )*(z[0]-x0)/(x1-x0);
    }
  }
  return interpProp;
}

} // namespace nalu
} // namespace Sierra
