#ifndef STKHELPERS_H
#define STKHELPERS_H

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_topology/topology.hpp>

#include <Realm.h>

namespace sierra {
namespace nalu {

inline
void populate_ghost_comm_procs(const stk::mesh::BulkData& bulk_data, stk::mesh::Ghosting& ghosting,
                            std::vector<int>& ghostCommProcs)
{
    ghostCommProcs.clear();

    std::vector<stk::mesh::EntityProc> sendList;
    ghosting.send_list(sendList);

    for(const stk::mesh::EntityProc& entProc : sendList) {
        stk::util::insert_keep_sorted_and_unique(entProc.second, ghostCommProcs);
    }

    std::vector<stk::mesh::EntityKey> recvList;
    ghosting.receive_list(recvList);

    for(const stk::mesh::EntityKey& key : recvList) {
        stk::mesh::Entity entity = bulk_data.get_entity(key);
        stk::util::insert_keep_sorted_and_unique(bulk_data.parallel_owner_rank(entity), ghostCommProcs);
    }
}

inline
stk::topology get_elem_topo(const Realm& realm, const stk::mesh::Part& surfacePart)
{
  std::vector<const stk::mesh::Part*> blockParts = realm.meta_data().get_blocks_touching_surface(&surfacePart);

  ThrowRequireMsg(blockParts.size() >= 1, "Error, expected at least 1 block for surface "<<surfacePart.name());

  stk::topology elemTopo = blockParts[0]->topology();
  if (blockParts.size() > 1) {
    for(size_t i=1; i<blockParts.size(); ++i) {
      ThrowRequireMsg(blockParts[i]->topology() == elemTopo,
                 "Error, found blocks of different topology connected to surface '"
                  <<surfacePart.name()<<"', "<<elemTopo<<" and "<<blockParts[i]->topology());
    }
  }

  ThrowRequireMsg(elemTopo != stk::topology::INVALID_TOPOLOGY,
                  "Error, didn't find valid topology block for surface "<<surfacePart.name());
  return elemTopo;
}

} // namespace nalu
} // namespace sierra

#endif /* STKHELPERS_H */
