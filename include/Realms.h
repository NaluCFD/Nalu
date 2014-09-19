/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Realms_h
#define Realms_h

#include <Enums.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <Realm.h>
#include <NaluParsing.h>

#include <map>
#include <string>
#include <algorithm>
#include <vector>

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

typedef std::vector<Realm *> RealmVector;

class Simulation;

class Realms : public RealmVector {
public:
  Realms(Simulation& sim) : simulation_(sim) {}

  ~Realms();

  void load(const YAML::Node & node) ;
  void breadboard();
  void initialize();
  Simulation *root();
  Simulation *parent();

  // helper
  struct IsString {
    IsString(std::string& str) : str_(str) {}
    std::string& str_;
    bool operator()(Realm *realm) { return realm->name_ == str_; }
  };

  Realm *find_realm(std::string realm_name)
  {
    RealmVector::iterator realm_iter = std::find_if(begin(), end(), IsString(realm_name));
    if (realm_iter != end())
      return *realm_iter;
    else
      return 0;
  }

  Simulation &simulation_;

};

} // namespace nalu
} // namespace Sierra

#endif
