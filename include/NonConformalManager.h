/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NonConformalManager_h
#define NonConformalManager_h

//==============================================================================
// Includes and forwards
//==============================================================================

// stk
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <vector>
#include <map>

namespace sierra {
namespace nalu {

class DgInfo;
class Realm;
class NonConformalInfo;

//=============================================================================
// Class Definition
//=============================================================================
// NonConformalManager
//=============================================================================
/**
 * * @par Description:
 * - class to manage all NonConformalInfo objects.
 *
 * @par Design Considerations:
 * -
 */
//=============================================================================
class NonConformalManager {

 public:

  // constructor and destructor
  NonConformalManager(
    Realm & realm,
    const bool ncAlgDetailedOutput );

  ~NonConformalManager();

  void initialize();
  void manage_ghosting();

  Realm &realm_;
  const bool ncAlgDetailedOutput_;

  /* ghosting for all surface:block pair */
  stk::mesh::Ghosting *nonConformalGhosting_;

  uint64_t needToGhostCount_;
 
  stk::mesh::EntityProcVec elemsToGhost_;
  std::vector<NonConformalInfo *> nonConformalInfoVec_;

};

} // end sierra namespace
} // end Acon namespace

#endif
