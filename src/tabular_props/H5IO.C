#include <tabular_props/H5IO.h>

#include <string.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

using std::endl;
using std::ostringstream;

namespace sierra {
namespace nalu {

//============================================================================
H5IO::H5IO()
  : file_( -1 ),
    group_( -1 ),
    fileVersion_( 0 )
{ }
//----------------------------------------------------------------------------
H5IO::~H5IO()
{ }
//----------------------------------------------------------------------------
void
H5IO::create_file( const std::string & name, int version )
{
  if ( file_ >= 0 ) {
    close_file();
  }

  // Increase the small data block size and metadata block size from the
  // default 2 kB.  This delays the onset of the O(N^2) write slowdown
  // experienced when writing large amounts of small data in large numbers
  // of groups.  Note that this will be the new minimum file size.
  //
  const hsize_t sd_block_size = 2097152;  // Small data block size
  const hsize_t md_block_size = 2097152;  // Metadata block size

  hid_t fapl_id = H5Pcreate( H5P_FILE_ACCESS );
  H5Pset_small_data_block_size( fapl_id, sd_block_size );
  H5Pset_meta_block_size( fapl_id, md_block_size );

  file_ = H5Fcreate( name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id );

  if ( file_ < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not create new HDF5 file for output: '" << name
           << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
  fileName_ = name;

  // Set the root group
  groupName_ = "/";

  // Tag the file with a version number
  fileVersion_ = version;
  write_attribute( "FileVersion", version );
}
//----------------------------------------------------------------------------
void
H5IO::open_file( const std::string & name )
{
  if ( file_ >= 0 ) {
    close_file();
  }

  // Open the file
  file_ = H5Fopen( name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  if ( file_ < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not open HDF5 file for input: '" << name
           << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
  fileName_ = name;

  // Set the root group
  groupName_ = "/";

  // Read the file version number, if present.  If not present (which
  // indicates files from before we started the tagging), then set a
  // version of 1.
  if ( has_attribute( "FileVersion" ) ) {
    read_attribute( "FileVersion", fileVersion_ );
  }
  else {
    fileVersion_ = 1;
  }
}
//----------------------------------------------------------------------------
void
H5IO::close_file()
{
  if ( file_ >= 0 ) {

    herr_t err = H5Fclose( file_ );

    if ( err < 0 ) {
      ostringstream errmsg;
      errmsg << "ERROR: Could not close HDF5 file: '" << fileName_
             << "'" << endl;
      throw std::runtime_error( errmsg.str() );
    }

    group_ = -1;
    file_ = -1;
    groupName_.clear();
    fileName_.clear();

  }
}
//----------------------------------------------------------------------------
H5IO
H5IO::create_group( const std::string & name )
{
  if ( file_ < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Cannot create HDF5 group '" << name << "'" << endl
           << "       before opening file." << endl;
    throw std::runtime_error( errmsg.str() );
  }

  // Create a new H5IO object that is a copy of this object
  H5IO newIO( *this );
  newIO.h5io_create_group( name );
  return newIO;
}
//----------------------------------------------------------------------------
H5IO
H5IO::open_group( const std::string & name )
{
  if ( file_ < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Cannot open HDF5 group '" << name << "'" << endl
           << "       before opening file." << endl;
    throw std::runtime_error( errmsg.str() );
  }

  // Create a new H5IO object that is a copy of this object
  H5IO newIO( *this );
  newIO.h5io_open_group( name );
  return newIO;
}
//----------------------------------------------------------------------------
unsigned int
H5IO::num_attributes()
{
  // Return the number of group attributes
  h5io_open_group();
  int n_attrs = H5Aget_num_attrs( group_ );
  h5io_close_group();

  if ( n_attrs < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not query number of attributes in group '"
           << groupName_ << "'" << endl
           << "       in file '" << fileName_ << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
  return n_attrs;
}
//----------------------------------------------------------------------------
void
H5IO::h5io_create_group( const std::string & name )
{
  hid_t newGroup = -1;
  std::string newGroupName;
  newGroupName = (name.at(0) == '/') ? name                     // Absolute
                                     : groupName_ + "/" + name; // Relative
  newGroup = H5Gcreate( file_, newGroupName.c_str(), 0 , H5P_DEFAULT, H5P_DEFAULT);

  if ( newGroup < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not create HDF5 group named '" << name << "'" <<endl
           << "       in file '" << fileName_ << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  group_ = newGroup;
  groupName_ = newGroupName;
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::h5io_open_group( const std::string & name )
{
  hid_t newGroup = -1;
  std::string newGroupName;
  newGroupName = (name.at(0) == '/') ? name                     // Absolute
                                     : groupName_ + "/" + name; // Relative
  newGroup = H5Gopen( file_, newGroupName.c_str() , H5P_DEFAULT );

  if ( newGroup < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not open HDF5 group named '" << name << "'" << endl
           << "       in file '" << fileName_ << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  group_ = newGroup;
  groupName_ = newGroupName;
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::h5io_open_group()
{
  group_ = H5Gopen( file_, groupName_.c_str() , H5P_DEFAULT );

  if ( group_ < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not open HDF5 group named '" << groupName_ << "'"
             << endl
           << "       in file '" << fileName_ << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
}
//----------------------------------------------------------------------------
void
H5IO::h5io_close_group()
{
  herr_t err = H5Gclose( group_ );
  group_ = -1;

  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not close HDF5 group named '" << groupName_ << "'"
             << endl
           << "       in file '" << fileName_ << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
}
//----------------------------------------------------------------------------
hid_t
H5IO::h5io_create_scalar()
{
  hid_t space_id = H5Screate( H5S_SCALAR );
  if ( space_id < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not create scalar dataspace in HDF5 group" << endl
           << "       '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
  return space_id;
}
//----------------------------------------------------------------------------
hid_t
H5IO::h5io_create_1D_array( unsigned int size )
{
  hsize_t hsize = size;
  hid_t space_id = H5Screate_simple( 1, &hsize, NULL );
  if ( space_id < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not create 1D array dataspace in HDF5 group" << endl
           << "       '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
  return space_id;
}

//----------------------------------------------------------------------------
hid_t
H5IO::h5io_create_attribute( const std::string & name, hid_t type, hid_t space )
{
  hid_t attr_id = H5Acreate( group_, name.c_str(), type, space, H5P_DEFAULT , H5P_DEFAULT );
  if ( attr_id < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not create attribute '" << name << "' in" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
  return attr_id;
}
//----------------------------------------------------------------------------
hid_t
H5IO::h5io_create_dataset( const std::string & name, hid_t type, hid_t space )
{
  hid_t data_id = H5Dcreate( group_, name.c_str(), type, space, H5P_DEFAULT , H5P_DEFAULT , H5P_DEFAULT );
  if ( data_id < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not create attribute '" << name << "' in" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }
  return data_id;
}

//----------------------------------------------------------------------------
void
H5IO::write_attribute( const std::string & name, int value )
{
  h5io_open_group();
  hid_t space_id = h5io_create_scalar();
  hid_t attr_id = h5io_create_attribute( name, H5T_NATIVE_INT, space_id );

  herr_t err = H5Awrite( attr_id, H5T_NATIVE_INT, &value );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write attribute '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Aclose( attr_id );
  H5Sclose( space_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::write_attribute( const std::string & name, unsigned int value )
{
  h5io_open_group();
  hid_t space_id = h5io_create_scalar();
  hid_t attr_id = h5io_create_attribute( name, H5T_NATIVE_UINT, space_id );

  herr_t err = H5Awrite( attr_id, H5T_NATIVE_UINT, &value );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write attribute '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Aclose( attr_id );
  H5Sclose( space_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::write_attribute( const std::string & name, double value )
{
  h5io_open_group();
  hid_t space_id = h5io_create_scalar();
  hid_t attr_id = h5io_create_attribute( name, H5T_NATIVE_DOUBLE, space_id );

  herr_t err = H5Awrite( attr_id, H5T_NATIVE_DOUBLE, &value );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write attribute '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Aclose( attr_id );
  H5Sclose( space_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::write_attribute( const std::string & name, const std::string & value )
{
  h5io_open_group();
  hid_t space_id = h5io_create_scalar();
  hid_t type_id = H5Tcopy( H5T_C_S1 );
  H5Tset_size( type_id, value.size() + 1 );
  hid_t attr_id = h5io_create_attribute( name, type_id, space_id );

  herr_t err = H5Awrite( attr_id, type_id, value.c_str() );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write attribute '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Aclose( attr_id );
  H5Tclose( type_id );
  H5Sclose( space_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::write_attribute( const std::string & name,
                       const std::vector<int> & value )
{
  h5io_open_group();
  hid_t space_id  = h5io_create_1D_array( value.size() );
  hid_t attr_id = h5io_create_attribute( name, H5T_NATIVE_INT, space_id );

  herr_t err = H5Awrite( attr_id, H5T_NATIVE_INT, &value[0] );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write attribute '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Aclose( attr_id );
  H5Sclose( space_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::write_attribute( const std::string & name,
                       const std::vector<unsigned int> & value )
{
  h5io_open_group();
  hid_t space_id  = h5io_create_1D_array( value.size() );
  hid_t attr_id = h5io_create_attribute( name, H5T_NATIVE_UINT, space_id );

  herr_t err = H5Awrite( attr_id, H5T_NATIVE_UINT, &value[0] );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write attribute '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Aclose( attr_id );
  H5Sclose( space_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::write_attribute( const std::string & name,
                       const std::vector<double> & value )
{
  h5io_open_group();
  hid_t space_id  = h5io_create_1D_array( value.size() );
  hid_t attr_id = h5io_create_attribute( name, H5T_NATIVE_DOUBLE, space_id );

  herr_t err = H5Awrite( attr_id, H5T_NATIVE_DOUBLE, &value[0] );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write attribute '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Aclose( attr_id );
  H5Sclose( space_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::write_attribute( const std::string & name,
                       const std::vector<std::string> & value )
{
  h5io_open_group();
  char ** buf = new char*[value.size()];
  for ( unsigned int i = 0; i < value.size(); ++i ) {
    buf[i] = new char[value[i].length() + 1];
    strcpy( buf[i], value[i].c_str() );
  }

  hid_t space_id  = h5io_create_1D_array( value.size() );
  hid_t type_id = H5Tcopy( H5T_C_S1 );
  H5Tset_size( type_id, H5T_VARIABLE );
  hid_t attr_id = h5io_create_attribute( name, type_id, space_id );

  herr_t err = H5Awrite( attr_id, type_id, buf );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write attribute '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Aclose( attr_id );
  H5Tclose( type_id );
  H5Sclose( space_id );

  for ( unsigned int i = 0; i < value.size(); ++i ) {
    delete [] buf[i];
  }
  delete [] buf;
  h5io_close_group();
}
//----------------------------------------------------------------------------
bool
H5IO::has_attribute( const std::string & name )
{
  const size_t NAMESIZE = 128;
  char nameBuf[NAMESIZE];
  bool foundAttribute = false;

  unsigned int numAttrs = num_attributes();

  h5io_open_group();
  for ( unsigned int i = 0; i < numAttrs; ++i ) {
    hid_t attr_id = H5Aopen_idx( group_, i );
    H5Aget_name( attr_id, NAMESIZE-1, nameBuf );
    H5Aclose( attr_id );
    if ( name == nameBuf ) {
      foundAttribute = true;
      break;
    }
  }
  h5io_close_group();

  return foundAttribute;
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( const std::string & name, int & value )
{
  h5io_open_group();
  hid_t attr_id = H5Aopen_name( group_, name.c_str() );
  H5Aread( attr_id, H5T_NATIVE_INT, &value );
  H5Aclose( attr_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( const std::string & name, unsigned int & value )
{
  h5io_open_group();
  hid_t attr_id = H5Aopen_name( group_, name.c_str() );
  H5Aread( attr_id, H5T_NATIVE_UINT, &value );
  H5Aclose( attr_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( const std::string & name, double & value )
{
  h5io_open_group();
  hid_t attr_id = H5Aopen_name( group_, name.c_str() );
  H5Aread( attr_id, H5T_NATIVE_DOUBLE, &value );
  H5Aclose( attr_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( const std::string & name, std::string & value )
{
  h5io_open_group();
  hid_t attr_id = H5Aopen_name( group_, name.c_str() );
  hid_t type_id = H5Aget_type( attr_id );
  hsize_t size = H5Tget_size( type_id );
  char * buf = new char[size];
  H5Aread( attr_id, type_id, buf );
  value = buf;
  delete [] buf;
  H5Tclose( type_id );
  H5Aclose( attr_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( unsigned int index, std::string & name,
                      std::string & value )
{
  const size_t NAMESIZE = 128;
  h5io_open_group();
  hid_t attr_id = H5Aopen_idx( group_, index );

  char nameBuf[NAMESIZE];
  H5Aget_name( attr_id, NAMESIZE-1, nameBuf );
  name = nameBuf;

  hid_t type_id = H5Aget_type( attr_id );
  hsize_t size = H5Tget_size( type_id );
  char * buf = new char[size];
  H5Aread( attr_id, type_id, buf );
  value = buf;
  delete [] buf;

  H5Tclose( type_id );
  H5Aclose( attr_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( const std::string & name,
                      std::vector<int> & value )
{
  h5io_open_group();
  hid_t attr_id = H5Aopen_name( group_, name.c_str() );
  hid_t space_id = H5Aget_space( attr_id );
  //int dimensions = H5Sget_simple_extent_ndims( space_id );  // Check == 1
  hsize_t size, maxsize;
  H5Sget_simple_extent_dims( space_id, &size, &maxsize );

  value.resize( size, 0 );
  H5Aread( attr_id, H5T_NATIVE_INT, &value[0] );

  H5Sclose( space_id );
  H5Aclose( attr_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( const std::string & name,
                      std::vector<unsigned int> & value )
{
  h5io_open_group();
  hid_t attr_id = H5Aopen_name( group_, name.c_str() );
  hid_t space_id = H5Aget_space( attr_id );
  //int dimensions = H5Sget_simple_extent_ndims( space_id );  // Check == 1
  hsize_t size, maxsize;
  H5Sget_simple_extent_dims( space_id, &size, &maxsize );

  value.resize( size, 0 );
  H5Aread( attr_id, H5T_NATIVE_UINT, &value[0] );

  H5Sclose( space_id );
  H5Aclose( attr_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( const std::string & name,
                      std::vector<double> & value )
{
  h5io_open_group();
  hid_t attr_id = H5Aopen_name( group_, name.c_str() );
  hid_t space_id = H5Aget_space( attr_id );
  //int dimensions = H5Sget_simple_extent_ndims( space_id );  // Check == 1
  hsize_t size, maxsize;
  H5Sget_simple_extent_dims( space_id, &size, &maxsize );

  value.resize( size, 0.0 );
  H5Aread( attr_id, H5T_NATIVE_DOUBLE, &value[0] );

  H5Sclose( space_id );
  H5Aclose( attr_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_attribute( const std::string & name,
                      std::vector<std::string> & value )
{
  h5io_open_group();
  hid_t attr_id = H5Aopen_name( group_, name.c_str() );
  hid_t type_id = H5Aget_type( attr_id );
  hid_t space_id = H5Aget_space( attr_id );
  //int dimensions = H5Sget_simple_extent_ndims( space_id );  // Check == 1
  hsize_t size, maxsize;
  H5Sget_simple_extent_dims( space_id, &size, &maxsize );

  char ** buf = new char*[size];
  H5Aread( attr_id, type_id, buf );

  value.resize( size, "" );
  for ( unsigned int i = 0; i < size; ++i ) {
    value[i] = buf[i];
  }

  H5Dvlen_reclaim( type_id, space_id, H5P_DEFAULT, buf );
  delete [] buf;
  H5Sclose( space_id );
  H5Tclose( type_id );
  H5Aclose( attr_id );
  h5io_close_group();
}

//----------------------------------------------------------------------------
void
H5IO::write_dataset( const std::string & name,
                     const std::vector<double> & value )
{
  h5io_open_group();
  hid_t space_id  = h5io_create_1D_array( value.size() );

  // create a dataset and write it
  hid_t data_id = h5io_create_dataset( name, H5T_NATIVE_DOUBLE, space_id );
  herr_t err = H5Dwrite( data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, &value[0] );
  if ( err < 0 ) {
    ostringstream errmsg;
    errmsg << "ERROR: Could not write dataset '" << name << "' to" << endl
           << "       HDF5 group '" << groupName_ << "' in file '" << fileName_
             << "'" << endl;
    throw std::runtime_error( errmsg.str() );
  }

  H5Dclose( data_id );
  H5Sclose( space_id );
  h5io_close_group();
}
//----------------------------------------------------------------------------
void
H5IO::read_dataset( const std::string & name, std::vector<double> & value )
{
  h5io_open_group();
  hid_t data_id = H5Dopen( group_, name.c_str(), H5P_DEFAULT );
  int size = H5Dget_storage_size( data_id ) / sizeof(double);
  value.resize( size, 0.0 );
  H5Dread( data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &value[0] );
  H5Dclose( data_id );
  h5io_close_group();
}

//----------------------------------------------------------------------------

} // end nalu namespace
} // end sierra namespace


