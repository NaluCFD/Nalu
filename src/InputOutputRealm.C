/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <InputOutputRealm.h>
#include <Realm.h>

// transfer
#include <xfer/Transfer.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// InputOutputRealm - holder of data to be transferred to/from
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
InputOutputRealm::InputOutputRealm(Realms& realms, const YAML::Node & node)
  : Realm(realms, node)
{
  // nothing now
}
  
//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
InputOutputRealm::~InputOutputRealm()
{
  for ( size_t k = 0; k < inputOutputFieldInfo_.size(); ++k ) 
    delete inputOutputFieldInfo_[k];
}
 
//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void 
InputOutputRealm::initialize()
{
  // bar minimum to register fields and to extract from possible mesh file
  register_io_fields();
  ioBroker_->populate_mesh();
  ioBroker_->populate_field_data();
  create_output_mesh();
  input_variables_from_mesh();
}

//--------------------------------------------------------------------------
//-------- register_io_fields ------------------------------------------------------
//--------------------------------------------------------------------------
void 
InputOutputRealm::register_io_fields() {
  // register fields; extract vector of field/part; only nodal for now
  for ( size_t k = 0; k < inputOutputFieldInfo_.size(); ++ k ) {
    const std::string fieldName = inputOutputFieldInfo_[k]->fieldName_;
    const int fieldSize = inputOutputFieldInfo_[k]->fieldSize_;
    const std::string fieldType = inputOutputFieldInfo_[k]->fieldType_;
    std::vector<std::string> targetNames = inputOutputFieldInfo_[k]->targetNames_;
    
    // sanity check on the type
    if ( fieldType != "node_rank" ) {
      throw std::runtime_error("Input/output realm only supports nodal_rank types");
    }
      
    // loop over target parts and declare/put the field
    for ( size_t j = 0; j < targetNames.size(); ++j ) {
      const std::string targetName = targetNames[j];
      stk::mesh::Part *targetPart = metaData_->get_part(targetName);
      if ( NULL == targetPart ) {
        throw std::runtime_error("Sorry, no part name found by the name: " + targetName + " for field: " + fieldName);
      }
      else {
        if ( fieldName == "velocity" ) {
          VectorFieldType *velocity = &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity"));
          stk::mesh::put_field(*velocity, *targetPart, fieldSize);
        }
        else {
          stk::mesh::FieldBase *theField 
            = &(metaData_->declare_field< stk::mesh::Field<double, stk::mesh::SimpleArrayTag> >(stk::topology::NODE_RANK, fieldName));
          stk::mesh::put_field(*theField,*targetPart,fieldSize);
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void 
InputOutputRealm::load(const YAML::Node & node)
{
  // call base class
  Realm::load(node);

  // now proceed with specific line commands to IO Realm
  const YAML::Node *y_field = node.FindValue("field_registration");
  if (y_field) {    
    
    // extract the sequence of types
    const YAML::Node *y_specs = expect_sequence(*y_field, "specifications", false);
    if (y_specs) {
      for (size_t ispec = 0; ispec < y_specs->size(); ++ispec) {
        const YAML::Node &y_spec = (*y_specs)[ispec];
        
        // find the name, size and type
        const YAML::Node *fieldNameNode = y_spec.FindValue("field_name");
        const YAML::Node *fieldSizeNode = y_spec.FindValue("field_size");
        const YAML::Node *fieldTypeNode = y_spec.FindValue("field_type");
        
        if ( NULL == fieldNameNode ) 
          throw std::runtime_error("Sorry, field name must be provided");
        
        if ( NULL == fieldSizeNode ) 
          throw std::runtime_error("Sorry, field size must be provided");
        
        if ( NULL == fieldTypeNode ) 
          throw std::runtime_error("Sorry, field type must be provided");
        
        // new the data
        InputOutputInfo *theInfo = new InputOutputInfo();

        // push data to the info object
        *fieldNameNode >> theInfo->fieldName_;
        *fieldSizeNode >> theInfo->fieldSize_;
        *fieldTypeNode >> theInfo->fieldType_;
        
        const YAML::Node &targets = y_spec["target_name"];
        if (targets.Type() == YAML::NodeType::Scalar) {
          theInfo->targetNames_.resize(1);
          targets >> theInfo->targetNames_[0];
        }
        else {
          theInfo->targetNames_.resize(targets.size());
          for (size_t i=0; i < targets.size(); ++i) {
            targets[i] >> theInfo->targetNames_[i];
          }
        }
        
        inputOutputFieldInfo_.push_back(theInfo);
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- populate_restart ------------------------------------------------
//--------------------------------------------------------------------------
double
InputOutputRealm::populate_restart(
  double &timeStepNm1, int &timeStepCount)
{
  return get_current_time();
}

} // namespace nalu
} // namespace Sierra
