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

#if defined (NALU_USES_PERCEPT)
#include <percept/FieldTypes.hpp>
#endif

namespace sierra{
namespace nalu{

class Realm;

class ErrorIndicatorAlgorithmDriver : public AlgorithmDriver
{
public:

  ErrorIndicatorAlgorithmDriver(
    Realm &realm);
  ~ErrorIndicatorAlgorithmDriver();

  void pre_work();
  void post_work();

public:

#if defined (NALU_USES_PERCEPT)
  GenericFieldType *errorIndicator_;
  percept::RefineFieldType *refineField_;
  percept::RefineFieldType *refineFieldOrig_;
  percept::RefineLevelType *refineLevelField_;
#else
  GenericFieldType *errorIndicator_;
  ScalarFieldType *refineField_;
  ScalarFieldType *refineFieldOrig_;
  ScalarFieldType *refineLevelField_;
#endif

  double maxErrorIndicator_;
};

} // namespace nalu
} // namespace Sierra

#endif

