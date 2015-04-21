#include <tabular_props/HDF5FilePtr.h>
#include <tabular_props/H5IO.h>

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
namespace sierra {
namespace nalu {

//=============================================================================
HDF5FilePtr::HDF5FilePtr( const std::string & fileName )
  : fileName_( fileName )
{
  // Increment the exported file version when making structural changes
  // to the HDF5 file layout, and modify importing code to handle all
  // supported versions.  Version key:
  //
  // 1  = Versions prior to version tagging
  // 2  = All arrays written as attributes instead of datasets to work around
  //      O(N^2) scaling behavior in write performance for very large numbers
  //      of small (less than about 2 kB) data blocks.
  // 3  = Unnecessary "SLFM" prefix removed from converters.  Backward
  //      compatibility added for reading old files.
  // 4  = Multiple mixture fraction support for Converters
  //
  exportFileVersion_ = 4;
  fileIO_ = new H5IO();
  //  fileIO_->open_file( fileName_ );
  read_hdf5();
}
//-----------------------------------------------------------------------------
HDF5FilePtr::~HDF5FilePtr()
{
  fileIO_->close_file();
  delete fileIO_;
}
//-----------------------------------------------------------------------------
bool
HDF5FilePtr::has_entry( const std::string & name ) const
{
  for (unsigned int i = 0; i < propertyNames_.size(); i++ ){
    if (propertyNames_[i] == name ) {
      return true;
    } 
  } 
  return false;
}
//-----------------------------------------------------------------------------
std::vector<std::string>
HDF5FilePtr::property_names() const
{
  return propertyNames_;
}
//-----------------------------------------------------------------------------
void
HDF5FilePtr::read_hdf5()

{
  fileIO_->open_file( fileName_ );
  fileIO_->read_attribute( "PropertyNames", propertyNames_ );
}
//--------------------------------------------------------------------
H5IO* HDF5FilePtr::get_H5IO()
{
  return fileIO_;
} 

} // end nalu namespace
} // end sierra namespace
