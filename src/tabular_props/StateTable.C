/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <tabular_props/StateTable.h>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// StateTable - holds the table
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
StateTable::StateTable(
  const std::string propertyTableName)
{
  // read in table; obvious room for improvement as we want a single copy of this guy
  std::ifstream input(propertyTableName.c_str());
  if ( input ) {
    std::string lineData;
    while ( getline(input, lineData) ) {
      double d;
      std::vector<double> row;
      std::stringstream lineStream(lineData);
      while ( lineStream >> d)
        row.push_back(d);
      table_.push_back(row);
    }
  }
  else {
    throw std::runtime_error("Unable to open the user provided table:");
  }
}

//--------------------------------------------------------------------------
//-------- destructor-------------------------------------------------------
//--------------------------------------------------------------------------
StateTable::~StateTable() {
  // nothing to do...
}

//--------------------------------------------------------------------------
//-------- get_entry -------------------------------------------------------
//--------------------------------------------------------------------------
std::vector<std::vector<double> > *
StateTable::get_entry(
  /*const std::string propertyName*/)
{
  return &table_;
}

} // namespace nalu
} // namespace Sierra
