/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef StateTable_h
#define StateTable_h

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class StateTable
{
public:

  StateTable(
    const std::string propertyTableName );
  ~StateTable();

  std::vector<std::vector<double > > table_;

  std::vector<std::vector<double> > *get_entry();

};

} // namespace nalu
} // namespace Sierra

#endif
