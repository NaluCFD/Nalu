#ifndef STKHELPERS_H
#define STKHELPERS_H

#include <stk_mesh/base/Ghosting.hpp>
#include <stk_util/util/SortAndUnique.hpp>

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

} // namespace nalu
} // namespace sierra

#endif /* STKHELPERS_H */
