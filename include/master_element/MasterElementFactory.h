/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MasterElementFactory_h
#define MasterElementFactory_h

#include <string>
#include <map>
#include <memory>

namespace stk { struct topology; }

namespace sierra{
namespace nalu{
  class MasterElement;

  struct MasterElementRepo
  {
  public:
    static MasterElement*
    get_surface_master_element(
      const stk::topology& theTopo);

    static MasterElement*
    get_volume_master_element(
      const stk::topology& theTopo);

    static MasterElement*
    get_fem_master_element(
      const stk::topology& theTopo);

    static void clear();
  private:
    MasterElementRepo() = default;
    // allow support of all three types of master elements in a given simulation
    static std::map<stk::topology, std::unique_ptr<MasterElement>> surfaceMeMap_;
    static std::map<stk::topology, std::unique_ptr<MasterElement>> volumeMeMap_;
    static std::map<stk::topology, std::unique_ptr<MasterElement>> femMeMap_;
  };

} // namespace nalu
} // namespace sierra

#endif
