/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SIMDINTERFACE_H
#define SIMDINTERFACE_H

/** \file
 *  \brief STK SIMD Interface
 *
 *  Nalu wrapper to provide SIMD capability for vectorizing sierra::nalu::Kernel
 *  algorithms.
 */

#include "stk_simd/Simd.hpp"

#include <vector>

namespace sierra {
namespace nalu {

typedef stk::simd::Double SimdDouble;

typedef SimdDouble DoubleType;

template<typename T>
using AlignedVector = std::vector<T, non_std::AlignedAllocator<T, 64>>;

using ScalarAlignedVector = AlignedVector<DoubleType>;

}  // nalu
}  // sierra

typedef sierra::nalu::DoubleType DoubleType;
#endif /* SIMDINTERFACE_H */
