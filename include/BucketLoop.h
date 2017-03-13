/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level NaluUnit      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef BucketLoop_h
#define BucketLoop_h

#include <stk_mesh/base/Bucket.hpp>

namespace sierra { namespace nalu {

  template<class LOOP_BODY>
  void bucket_loop(const stk::mesh::BucketVector& buckets, LOOP_BODY inner_loop_body)
  {
    for (const stk::mesh::Bucket* bptr : buckets)    {
      const stk::mesh::Bucket& bkt = *bptr;
      for(size_t j=0; j<bkt.size(); ++j)  {
        inner_loop_body(bkt[j]);
      }
    }
  }

}}

#endif
