/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Realm.h>
#include <Realms.h>
#include <TimeIntegrator.h>
#include <Simulation.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <NaluParsing.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Realms - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
  
//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Realms::~Realms()
{
  for (size_t ir = 0; ir < this->size(); ++ir)
    delete (*this)[ir];
}

void 
Realms::load(const YAML::Node & node) 
{
  const YAML::Node *realms = node.FindValue("realms");
  if (realms)
  {
    for ( size_t irealm = 0; irealm < realms->size(); ++irealm )
    {
      const YAML::Node & realm_node = (*realms)[irealm];
      Realm *realm = new Realm(*this);
      realm->load(realm_node);
      this->push_back(realm);
    }
  }
  else
    throw std::runtime_error("parser error Realms::load");
}
  
void 
Realms::breadboard()
{
  for ( size_t irealm = 0; irealm < this->size(); ++irealm ) {
    (*this)[irealm]->breadboard();
  }
}

void 
Realms::initialize()
{
#if 0
  for ( size_t irealm = 0; irealm < this->size(); ++irealm )
    {
      (*this)[irealm]->initialize();
    }
#endif
}

Simulation *Realms::root() { return parent()->root(); }
Simulation *Realms::parent() { return &simulation_; }

} // namespace nalu
} // namespace Sierra
