/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ErrorIndicatorAlgorithmDriver_h
#define ErrorIndicatorAlgorithmDriver_h

#include <AlgorithmDriver.h>
#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

#if defined (NALU_USES_PERCEPT)

class ErrorIndicatorAlgorithmDriver : public AlgorithmDriver
{
public:

  ErrorIndicatorAlgorithmDriver(
    const Realm &realm);
  ~ErrorIndicatorAlgorithmDriver();

  void pre_work();
  void post_work();

public:
  GenericFieldType *errorIndicator_;
  percept::RefineFieldType *refineField_;
  percept::RefineFieldType *refineFieldOrig_;
  percept::RefineLevelType *refineLevelField_;

  double maxErrorIndicator_;
};

#else

class ErrorIndicatorAlgorithmDriver : public AlgorithmDriver
{
public:
  
 ErrorIndicatorAlgorithmDriver( const Realm &realm) : AlgorithmDriver(realm) {}
  ~ErrorIndicatorAlgorithmDriver() {}

  void pre_work() {}
  void post_work() {}
};

#endif
} // namespace nalu
} // namespace Sierra

#endif
