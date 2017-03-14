/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AlgTraits_h
#define AlgTraits_h

#include <stk_topology/topology.hpp>

namespace sierra {
namespace nalu {

// limited supported now (P=1 3D elements)
struct AlgTraitsHex8 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerElement_ = 8;
  static constexpr int numScsIp_ = 12;
  static constexpr int numScvIp_ = 8;
  static constexpr int numGp_ = 8; // for FEM
  static constexpr stk::topology::topology_t topo_ = stk::topology::HEX_8;
};

struct AlgTraitsHex27 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerElement_ = 27;
  static constexpr int numScsIp_ = 216;
  static constexpr int numScvIp_ = 216;
  static constexpr int numGp_ = 27; // for FEM (not supported)
  static constexpr stk::topology::topology_t topo_ = stk::topology::HEX_27;
};

struct AlgTraitsTet4 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerElement_ = 4;
  static constexpr int numScsIp_ = 6;
  static constexpr int numScvIp_ = 4;
  static constexpr int numGp_ = 4; // for FEM (not supported)
  static constexpr stk::topology::topology_t topo_ = stk::topology::TET_4;
};

struct AlgTraitsPyr5 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerElement_ = 5;
  static constexpr int numScsIp_ = 8;
  static constexpr int numScvIp_ = 5;
  static constexpr int numGp_ = 5; // for FEM (not supported)
  static constexpr stk::topology::topology_t topo_ = stk::topology::PYRAMID_5;
};

struct AlgTraitsWed6 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerElement_ = 6;
  static constexpr int numScsIp_ = 9;
  static constexpr int numScvIp_ = 6;
  static constexpr int numGp_ = 6; // for FEM (not supported)
  static constexpr stk::topology::topology_t topo_ = stk::topology::WEDGE_6;
};

struct AlgTraitsQuad4_2D {
  static constexpr int nDim_ = 2;
  static constexpr int nodesPerElement_ = 4;
  static constexpr int numScsIp_ = 4;
  static constexpr int numScvIp_ = 4;
  static constexpr int numGp_ = 4; // for FEM (not supported)
  static constexpr stk::topology::topology_t topo_ = stk::topology::QUAD_4_2D;
};

struct AlgTraitsTri3_2D {
  static constexpr int nDim_ = 2;
  static constexpr int nodesPerElement_ = 3;
  static constexpr int numScsIp_ = 3;
  static constexpr int numScvIp_ = 3;
  static constexpr int numGp_ = 3; // for FEM (not supported)
  static constexpr stk::topology::topology_t topo_ = stk::topology::TRI_3_2D;
};

} // namespace nalu
} // namespace Sierra

#endif
