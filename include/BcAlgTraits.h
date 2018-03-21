/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef BcAlgTraits_h
#define BcAlgTraits_h

#include <stk_topology/topology.hpp>

namespace sierra {
namespace nalu {

struct BcAlgTraitsHex8Quad4 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerElement_ = 8;
  static constexpr int nodesPerFace_ = 4;
  static constexpr int numFaceIp_ = 4;
  static constexpr int numScsIp_ = 12;
  static constexpr int numScvIp_ = 8;
  static constexpr stk::topology::topology_t faceTopo_ = stk::topology::QUAD_4;
  static constexpr stk::topology::topology_t elemTopo_ = stk::topology::HEX_8;
};

struct BcAlgTraitsQuad4 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerFace_ = 4;
  static constexpr int numFaceIp_ = 4;
  static constexpr stk::topology::topology_t faceTopo_ = stk::topology::QUAD_4;
};

struct BcAlgTraitsQuad9 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerFace_ = 9;
  static constexpr int numFaceIp_ = 36;
  static constexpr stk::topology::topology_t faceTopo_ = stk::topology::QUAD_9;
};

struct BcAlgTraitsTri3 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerFace_ = 3;
  static constexpr int numFaceIp_ = 3;
  static constexpr stk::topology::topology_t faceTopo_ = stk::topology::TRI_3;
};

struct BcAlgTraitsLine2 {
  static constexpr int nDim_ = 2;
  static constexpr int nodesPerFace_ = 2;
  static constexpr int numFaceIp_ = 2;
  static constexpr stk::topology::topology_t faceTopo_ = stk::topology::LINE_2;
};

struct BcAlgTraitsLine3 {
  static constexpr int nDim_ = 2;
  static constexpr int nodesPerFace_ = 3;
  static constexpr int numFaceIp_ = 6;
  static constexpr int numGp_ = 2;
  static constexpr stk::topology::topology_t faceTopo_ = stk::topology::LINE_3;
};


} // namespace nalu
} // namespace Sierra

#endif
