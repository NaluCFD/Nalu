#include <NaluEnv.h>
#include <tabular_props/HDF5Table.h>
#include <tabular_props/Converter.h>
#include <tabular_props/H5IO.h>
#include <tabular_props/BSpline.h>

#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace sierra {
namespace nalu {

//============================================================================
HDF5Table::HDF5Table()
  : fileIO_( ),
    indVarSize_( 0 ),
    dimension_( 0 ),
    valueMin_( 0.0 ),
    valueMax_( 0.0 ),
    spline_(  ),
    clipEventLogSize_( 10 ),
    numClipped_( 0 )
{
}

//============================================================================
HDF5Table::HDF5Table(
   H5IO *fileIO,
   std::string tablePropName,
   std::vector<std::string> &indVarNameVec,
   std::vector<std::string> &indVarTableNameVec)
  : fileIO_( fileIO ),
    tablePropName_(tablePropName),
    indVarTableNameVec_(indVarTableNameVec),
    indVarSize_(indVarNameVec.size()),
    dimension_( 0 ),
    valueMin_( 0.0 ),
    valueMax_( 0.0 ),
    spline_( NULL ),
    clipEventLogSize_( 10 ),
    numClipped_( 0 )
{ 
  // extract the independent fields; check if there is one..
  if ( indVarSize_ == 0 )
    throw std::runtime_error("HDF5Table: independent variable size is zero:");

  //read in table
  read_hdf5_property();
}
//----------------------------------------------------------------------------
HDF5Table::~HDF5Table()
{
  for ( unsigned int i = 0; i < converters_.size(); ++i ) {
    delete converters_[i];
  }
  converters_.clear();
}
//----------------------------------------------------------------------------
void
HDF5Table::add_converter( const Converter * converter )
{
  converters_.push_back( converter );
  update_input_mapping();
}
//----------------------------------------------------------------------------
void
HDF5Table::update_input_mapping()
{
  // Start with a clean, direct-through mapping of the table
  //
  directInputIndex_.resize( inputNames_.size() );
  for ( unsigned int i = 0; i < directInputIndex_.size(); ++i ) {
    directInputIndex_[i] = i;
  }
  convTableIndex_.clear();
  convInputIndex_.clear();

  // Add the converters one-by-one until we have our final mapping
  //
  // Find the converter output variable in the list of inputNames_
  // (tableInputNames is same as inputNames_ at this point)
  const std::vector<std::string> & tableInputNames = input_names();
  for ( unsigned int i = 0; i < converters_.size(); ++i ) {

    // Check that the output of this converter, converters_[i]->name()
    // matches one of the table inputs in tableInputNames = inputNames_
    int idx = findix( tableInputNames, converters_[i]->name() );
    if ( idx < 0 ) {
      std::ostringstream errmsg;
      errmsg
        << "ERROR: Unable to find the Converter output variable '"
          << converters_[i]->name() << "'" << std::endl
        << "       in the list of table inputs:" << std::endl;
      for ( unsigned int j = 0; j < tableInputNames.size(); ++j ) {
        errmsg << "        - " << tableInputNames[j] << std::endl;
      }
      throw std::runtime_error( errmsg.str() );
    }
    // store the idx (based on ordering in indexNames_) of converter output
    // This is later used to provide input for this idx in lookupBuffer_
    // Note that below we add the inputs required by the converter
    convTableIndex_.push_back( idx );

    // Pluck this variable out of the global input variable list and
    // direct input list
    //
    idx = findix( inputNames_, converters_[i]->name() );
    if ( idx < 0 ) {
      std::ostringstream errmsg;
      errmsg
        << "ERROR: You may not hook the output of more than one" << std::endl
        << "       Converter to the same Table variable.  Error" << std::endl
        << "       occurred with the '" << converters_[i]->name() << "'"
        << " Converter." << std::endl;
      throw std::runtime_error( errmsg.str() );
    }
    inputNames_.erase( inputNames_.begin() + idx );
    directInputIndex_.erase( directInputIndex_.begin() + idx );

    // Add any new inputs required by this converter to the global list
    // if they aren't already there
    const std::vector<std::string> & convInputNames =
                                     converters_[i]->input_names();
    for ( unsigned int j = 0; j < convInputNames.size(); ++j ) {
      idx = findix( inputNames_, convInputNames[j] );
      if ( idx < 0 ) {
        inputNames_.push_back( convInputNames[j] );
      }
    }
  }
  // compute the new size of independent variable list with converter inputs
  dimension_ = inputNames_.size();

  // Now that our global input variable list is finalized, grab the
  // indexes into it for each of our converters
  //
  for ( unsigned int i = 0; i < converters_.size(); ++i ) {
    std::vector<unsigned int> convIndex;
    const std::vector<std::string> & convInputNames =
                                     converters_[i]->input_names();
    for ( unsigned int j = 0; j < convInputNames.size(); ++j ) {
      int idx = findix( inputNames_, convInputNames[j] );
      convIndex.push_back( idx );
    }
    // convInputIndex_ has an outer loop over the converters 
    // and inner loop over convInputNames 
    convInputIndex_.push_back( convIndex );
  }

  // Make sure our input buffers for the query() function are allocated
  // to the correct size
  //
  lookupBuffer_.resize( dimension_ );
  lookupBufferChecked_.resize( dimension_ );

  unsigned int maxSize = 0;
  for ( unsigned int i = 0; i < convInputIndex_.size(); ++i ) {
    maxSize = std::max( maxSize, (unsigned int)convInputIndex_[i].size() );
  }
  converterBuf_.resize( maxSize );

  // Now we have finalized the list of inputNames_ as it is expected 
  // by the spline query.  We need to set up a mapping between 
  // indVarTableNameVec_ and inputNames_.
  if ( inputNames_.size() != indVarTableNameVec_.size() ) {
    std::ostringstream errmsg;
    errmsg
      << "ERROR: The number of input variables needed by HDF5Table " << std::endl
      << "is not the same as that defined by indVarTableNameVec_."<< std::endl;
    throw std::runtime_error( errmsg.str() );
  }
  
  indexIndVar_.clear();
  for ( unsigned int i = 0; i < inputNames_.size() ; i++ ) {
    int idx = findix( indVarTableNameVec_, inputNames_[i] );
    if ( idx < 0 ) {
      std::ostringstream errmsg;
      errmsg
        << "ERROR: Unable to find inputNames_ = " << inputNames_[i] << std::endl
	<< "in indVarTableNameVec_"  << std::endl;
      throw std::runtime_error( errmsg.str() );
    }
    indexIndVar_.push_back( idx );
  }


}
//----------------------------------------------------------------------------
// Returns index of second argument in first argument vector
int HDF5Table::findix( const std::vector<std::string> & nameVector,
                   const std::string & name )
{
  std::vector<std::string>::const_iterator iname;
  iname = std::find( nameVector.begin(), nameVector.end(), name );
  if ( iname != nameVector.end() ) {
    return iname - nameVector.begin();
  }
  else {
    return -1;
  }
}
//----------------------------------------------------------------------------
double
HDF5Table::query( const std::vector<double> &inputs ) const
{
  if ( converters_.size() == 0 ) {

    // No converters, so just do a direct lookup in the table
    for ( unsigned int i = 0; i < indexIndVar_.size() ; i++ ) {
      lookupBuffer_[i] = inputs[indexIndVar_[i]];
    }
  }
  else {

    // Set the table input variables that we already know.  These (if any)
    // will be at the beginning of inputs.
    //
    for ( unsigned int i = 0; i < directInputIndex_.size(); ++i ) {
      lookupBuffer_[directInputIndex_[i]] = inputs[i];
    }

    // Process each of our converters, and place their result in the correct
    // spot in the table buffer.
    //
    unsigned int tableIdx = directInputIndex_.size();
    for ( unsigned int i = 0; i < converters_.size(); ++i, ++tableIdx ) {
      for ( unsigned int j = 0; j < convInputIndex_[i].size(); ++j ) {
        converterBuf_[j] = inputs[convInputIndex_[i][j]];
      }
      lookupBuffer_[convTableIndex_[i]] = converters_[i]->query( converterBuf_ );
    }
  }

  bool clipped = false;
  for ( unsigned int i = 0; i < dimension_; ++i ) {
    lookupBufferChecked_[i] = lookupBuffer_[i];
    
    if ( lookupBufferChecked_[i] < inputMin_[i] ) {
      clipped = true;
      lookupBufferChecked_[i] = inputMin_[i];
    }
    
    if ( lookupBufferChecked_[i] > inputMax_[i] ) {
      clipped = true;
      lookupBufferChecked_[i] = inputMax_[i];
    }
    
    // Convert to log scale if required
    if ( inputLogScale_[i] == 1 ) {
      lookupBufferChecked_[i] = std::log( std::max(lookupBufferChecked_[i], 1.e-16) );
    }
  }
  
  if ( clipped ) {
    //
    // Increment the clipping counter and log this set of input coordinates
    // for later diagnostic output
    //
    ++numClipped_;
    if ( clipEventLogSize_ > 0 ) {
      log_clip_event( lookupBuffer_ );
    }
  }
  
  // Perform the query
  return spline_->value( lookupBufferChecked_ );
}
//----------------------------------------------------------------------------
double
HDF5Table::raw_query( const std::vector<double> &inputs ) const
{
  if ( converters_.size() == 0 ) {

    // No converters, so just do a direct lookup in the table.  Note that
    // no input bounds checking is done, so extrapolation is allowed.
    //
    lookupBuffer_ = inputs;
  }
  else {

    // Set the table input variables that we already know.  These (if any)
    // will be at the beginning of inputs.
    //
    for ( unsigned int i = 0; i < directInputIndex_.size(); ++i ) {
      lookupBuffer_[directInputIndex_[i]] = inputs[i];
    }

    // Process each of our converters, and place their result in the correct
    // spot in the table buffer.
    //
    unsigned int tableIdx = directInputIndex_.size();
    for ( unsigned int i = 0; i < converters_.size(); ++i, ++tableIdx ) {
      for ( unsigned int j = 0; j < convInputIndex_[i].size(); ++j ) {
        converterBuf_[j] = inputs[convInputIndex_[i][j]];
      }
      lookupBuffer_[convTableIndex_[i]] = converters_[i]->query( converterBuf_ );
    }
  }

  // Do our final table lookup, and return the result.  Note that no
  // input bounds checking is done, so extrapolation is allowed.
  //  
  return spline_->value( &lookupBuffer_[0] );

}
//--------------------------------------------------------------------
void
HDF5Table::set_clipping_log_size( unsigned int size ) 
{
  clipEventLogSize_ = size ;
}
//--------------------------------------------------------------------
unsigned int
HDF5Table::num_clipping_events() const
{
  return numClipped_;
}
//--------------------------------------------------------------------
const ClipEventLog &
HDF5Table::clipping_event_log() const
{
  return clipEventLog_;
}
//--------------------------------------------------------------------
const std::vector<double> &
HDF5Table::clipping_event_min_bounds() const
{
  return inputMin_;
}
//--------------------------------------------------------------------
const std::vector<double> &
HDF5Table::clipping_event_max_bounds() const
{
  return inputMax_;
}
//--------------------------------------------------------------------
void
HDF5Table::clear_clipping_log() 
{
  numClipped_ = 0;
  clipEventLog_.clear();
}
//----------------------------------------------------------------------------
void
HDF5Table::log_clip_event( const std::vector<double> & values ) const
{
  double sev = 1.0;
  for ( unsigned int i = 0; i < inputMin_.size(); ++i ) {
    if ( values[i] < inputMin_[i] ) {
      // Below lower limit
      sev *= 1.0 + (inputMin_[i] - values[i]) / (inputMax_[i] - inputMin_[i]);
    }
    else  if ( values[i] > inputMax_[i] ) {
      // Above upper limit
      sev *= 1.0 + (values[i] - inputMax_[i]) / (inputMax_[i] - inputMin_[i]);
    }
  }
  ClipEvent event;
  event.severity = sev;
  event.values = values;
  clipEventLog_.insert( event );

  // Remove the event with smallest severity (the last event in the set)
  // if the above insertion pushed us over the log size limit.
  if ( clipEventLog_.size() > clipEventLogSize_ ) {
    clipEventLog_.erase( --(clipEventLog_.end()) );
  }
}
//--------------------------------------------------------------------
bool
HDF5Table::has_attribute( const std::string & name ) const
{
  Attributes::const_iterator iatr;
  iatr = attributes_.find( name );
  return ( iatr != attributes_.end() );
}
//--------------------------------------------------------------------
const std::string &
HDF5Table::attribute( const std::string & name ) const
{
  Attributes::const_iterator iatr;
  iatr = attributes_.find( name );
  if ( iatr != attributes_.end() ) {
    return (*iatr).second;
  }
  else {
    static std::string empty;
    return empty;
  }
}
//--------------------------------------------------------------------
void
HDF5Table::read_hdf5_property( )
{

  // Won't need H5IO to be passed, instead get it in the constructor at the file level
  H5IO propIO = fileIO_->open_group( tablePropName_ );
  //
  // Read the data describing this property configuration
  //
  propIO.read_attribute( "Name", name_ );
  propIO.read_attribute( "Dimension", dimension_ );
  propIO.read_attribute( "InputNames", inputNames_ );

  // Read the Converters (if any)
  //
  unsigned int nConverters = 0;
  propIO.read_attribute( "NConverters", nConverters );

  for ( unsigned int i = 0; i < nConverters; ++i ) {
    std::ostringstream label;
    label << "Converter_" << i;

    H5IO converterIO = propIO.open_group( label.str() );
    std::string converterType;
    converterIO.read_attribute( "ConverterType", converterType );

    ConverterFactory converterFactory;
    Converter * converter = converterFactory.create( converterType );

    converter->read_hdf5( converterIO );
    converters_.push_back( converter );
  }

  // Read the central Table for the dependent variable
  H5IO tableIO = propIO.open_group( "Table" );
  read_hdf5_table( tableIO );
  //
  // Wire up all of the inputs and outputs between the Table and Converters
  //
  update_input_mapping();

}
//============================================================================
//--------------------------------------------------------------------
void
HDF5Table::read_hdf5_table( H5IO & io )
{

    //
    // Read the data describing this state variable
    //
    io.read_attribute( "Name", name_ );
    io.read_attribute( "Dimension", dimension_ );
    io.read_attribute( "InputNames", inputNames_ );
    io.read_attribute( "InputLogScale", inputLogScale_ );
    io.read_attribute( "InputMin", inputMin_ );
    io.read_attribute( "InputMax", inputMax_ );
    io.read_attribute( "MeshMin", meshMin_ );
    io.read_attribute( "MeshMax", meshMax_ );
    io.read_attribute( "ValueMin", valueMin_ );
    io.read_attribute( "ValueMax", valueMax_ );
    
    //
    // Read any arbitrary metadata.
    //
    H5IO attrIO = io.open_group( "Attributes" );
    unsigned int numAttrs = attrIO.num_attributes();
    for ( unsigned int i = 0; i < numAttrs; ++i ) {
      std::string name;
      std::string value;
      attrIO.read_attribute( i, name, value );
      attributes_[name] = value;
    }
    
    //
    // Read the independent variable mesh points
    //
    for ( unsigned int i = 0; i < dimension_; ++i ) {
      std::ostringstream label;
      label << "Mesh_" << i;
      std::vector<double> mesh;
      if ( io.file_version() >= 2 ) {
	// Written as an attribute to improve performance
	io.read_attribute( label.str(), mesh );
      }
      else {
	io.read_dataset( label.str(), mesh );
      }
      mesh_.push_back( mesh );
    }
    
    //
    // Construct a skeletal spline, and then read it in
    //
    const bool allowClipping = false; // The Table is responsible for clipping
    
    switch( dimension_ ){
    case 1:
      spline_ = new BSpline1D( allowClipping );
      break;
    case 2:
      spline_ = new BSpline2D( allowClipping );
      break;
    case 3:
      spline_ = new BSpline3D( allowClipping );
      break;
    case 4:
      spline_ = new BSpline4D( allowClipping );
      break;
    case 5:
      spline_ = new BSpline5D( allowClipping );
      break;
    default:
      std::ostringstream errmsg;
      errmsg << "ERROR: unsupported dimension for BSpline creation!";
      throw std::runtime_error( errmsg.str() );
    }
    
    H5IO splineIO = io.open_group( "BSpline" );
    spline_->read_hdf5( splineIO );
    
    lookupBuffer_.resize( dimension_ );
    lookupBufferChecked_.resize( dimension_ );

}
//============================================================================

} // end nalu namespace
} // end sierra namespace

