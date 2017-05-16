/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef BuildTemplates_h
#define BuildTemplates_h

#include <AlgTraits.h>

namespace sierra{
namespace nalu{

#define INSTANTIATE_KERNEL(ClassName)                             \
template class ClassName<AlgTraitsHex8>;                          \
template class ClassName<AlgTraitsHex27>;                         \
template class ClassName<AlgTraitsTet4>;                          \
template class ClassName<AlgTraitsPyr5>;                          \
template class ClassName<AlgTraitsWed6>;                          \
template class ClassName<AlgTraitsQuad4_2D>;                      \
template class ClassName<AlgTraitsQuad9_2D>;                      \
template class ClassName<AlgTraitsTri3_2D>;                       \

} // namespace nalu
} // namespace Sierra

#endif
