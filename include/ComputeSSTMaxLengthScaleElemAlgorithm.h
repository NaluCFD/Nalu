/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeSSTMaxLengthScaleElemAlgorithm_h
#define ComputeSSTMaxLengthScaleElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class ComputeSSTMaxLengthScaleElemAlgorithm : public Algorithm
{
public:
  ComputeSSTMaxLengthScaleElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~ComputeSSTMaxLengthScaleElemAlgorithm() {}

  virtual void execute();
  
  VectorFieldType *coordinates_;
  ScalarFieldType *maxLengthScale_;
};

} // namespace nalu
} // namespace Sierra

#endif
