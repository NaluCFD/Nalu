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

#define INSTANTIATE_KERNEL_3D(ClassName)                          \
template class ClassName<AlgTraitsHex8>;                          \
template class ClassName<AlgTraitsHex27>;                         \
template class ClassName<AlgTraitsTet4>;                          \
template class ClassName<AlgTraitsPyr5>;                          \
template class ClassName<AlgTraitsWed6>;                          \

#define INSTANTIATE_KERNEL_2D(ClassName)                          \
template class ClassName<AlgTraitsQuad4_2D>;                      \
template class ClassName<AlgTraitsQuad9_2D>;                      \
template class ClassName<AlgTraitsTri3_2D>;                       \

#define INSTANTIATE_KERNEL_3D_HO(ClassName)                       \
template class ClassName<AlgTraitsHexGL<2>>;                      \
template class ClassName<AlgTraitsHexGL<3>>;                      \
template class ClassName<AlgTraitsHexGL<4>>;                      \
template class ClassName<AlgTraitsHexGL<5>>;                      \

#define INSTANTIATE_KERNEL_2D_HO(ClassName)                       \
template class ClassName<AlgTraitsQuadGL<2>>;                     \
template class ClassName<AlgTraitsQuadGL<3>>;                     \
template class ClassName<AlgTraitsQuadGL<4>>;                     \
template class ClassName<AlgTraitsQuadGL<5>>;                     \

#define INSTANTIATE_KERNEL(ClassName)           \
  INSTANTIATE_KERNEL_3D(ClassName)              \
  INSTANTIATE_KERNEL_2D(ClassName)              \
  INSTANTIATE_KERNEL_3D_HO(ClassName)           \
  INSTANTIATE_KERNEL_2D_HO(ClassName)           \


} // namespace nalu
} // namespace Sierra

#endif
