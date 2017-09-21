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
      const stk::topology& theTopo,
      int dimension = 0,
      std::string quadType = "GaussLegendre");

    static MasterElement*
    get_volume_master_element(
      const stk::topology& theTopo,
      int dimension = 0,
      std::string quadType = "GaussLegendre");

    static void clear();
  private:
    MasterElementRepo() = default;
    static std::map<stk::topology, std::unique_ptr<MasterElement>> surfaceMeMap_;
    static std::map<stk::topology, std::unique_ptr<MasterElement>> volumeMeMap_;
  };

} // namespace nalu
} // namespace sierra

#endif
