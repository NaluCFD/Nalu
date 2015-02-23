/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContactManager_h
#define ContactManager_h

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

class Realm;
class ContactInfo;

//=============================================================================
// Class Definition
//=============================================================================
// ContactManager
//=============================================================================
/**
 * * @par Description:
 * - class to manage all ContactInfo objects.
 *
 * @par Design Considerations:
 * -
 */
//=============================================================================
class ContactManager {

 public:

  // constructor and destructor
  ContactManager(
    Realm & realm);

  ~ContactManager();

  void initialize();
  void manage_ghosting();

  Realm &realm_;

  /* ghosting for all surface:block pair */
  stk::mesh::Ghosting *contactGhosting_;

  uint64_t needToGhostCount_;
  bool provideDetailedOutput_;

  stk::mesh::EntityProcVec elemsToGhost_;
  std::vector<ContactInfo *> contactInfoVec_;

};

} // end sierra namespace
} // end Acon namespace

#endif
