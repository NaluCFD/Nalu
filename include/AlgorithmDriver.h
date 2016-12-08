/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AlgorithmDriver_h
#define AlgorithmDriver_h

#include<Enums.h>

#include<map>

namespace sierra{
namespace nalu{

class Realm;
class Algorithm;

class AlgorithmDriver
{
public:

  AlgorithmDriver(
    Realm &realm);
  virtual ~AlgorithmDriver();

  virtual void pre_work(){};
  virtual void execute();
  virtual void post_work(){};

  Realm &realm_;
  std::map<AlgorithmType, Algorithm *> algMap_;
};

} // namespace nalu
} // namespace Sierra

#endif
