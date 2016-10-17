#include <gtest/gtest.h>
#include <limits>
#include <ostream>
#include <sstream>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_topology/topology.hpp>

#include <NonConformalManager.h>

#include "UnitTestUtils.h"

namespace {

static const std::string ghostingName = "myGhosting";

stk::mesh::Ghosting* find_ghosting_by_name(const std::string& name, const stk::mesh::BulkData& bulk)
{
    const std::vector<stk::mesh::Ghosting*>& ghostings = bulk.ghostings();
    for(stk::mesh::Ghosting* ghost_ptr : ghostings) {
        if (ghost_ptr->name() == name) {
            return ghost_ptr;
        }
    }
    return (stk::mesh::Ghosting*)nullptr;
}

void ghostElems(stk::mesh::BulkData& bulk, const stk::mesh::EntityProcVec& elemsToGhost)
{
    bulk.modification_begin();

    stk::mesh::Ghosting* myGhosting = find_ghosting_by_name(ghostingName, bulk);
    if (myGhosting == nullptr) {
        myGhosting = &bulk.create_ghosting(ghostingName);
    }

    bulk.change_ghosting(*myGhosting, elemsToGhost);

    bulk.modification_end();
}

void convert_elem_ids_to_elemprocs(const stk::mesh::BulkData& bulk, int proc, const std::vector<stk::mesh::EntityId>& elemIds,
                               stk::mesh::EntityProcVec& elemsToGhost)
{
    elemsToGhost.clear();
    for(stk::mesh::EntityId elemId : elemIds) {
        elemsToGhost.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK, elemId), proc));
    }
}

bool compare_vectors(const std::string& prefix, const stk::mesh::BulkData& bulk,
                     const stk::mesh::EntityProcVec& gold, const stk::mesh::EntityProcVec& vec)
{
    if (gold.size() != vec.size()) {
        std::cerr<<"P"<<bulk.parallel_rank()<<" "<<prefix<<" gold size="<<gold.size()<<", size="<<vec.size()<<std::endl;
        for(const stk::mesh::EntityProc& entityProc : vec) {
            std::cerr<<"P"<<bulk.parallel_rank()<<" "<<prefix<<" val: "<<bulk.entity_key(entityProc.first)<<",proc="<<entityProc.second<<std::endl;
        }
        return false;
    }

    bool returnVal = true;
    std::ostringstream os;
    for(size_t i=0; i<gold.size(); ++i) {
        if (gold[i] != vec[i]) {
            os<<prefix<<" "<<i<<"-th value mismatch: gold="<<bulk.entity_key(gold[i].first)<<",proc="<<gold[i].second<<", val="<<bulk.entity_key(vec[i].first)<<",proc="<<vec[i].second<<std::endl;
            returnVal = false;
        }
    }
    std::cerr<<os.str();
    return returnVal;
}

bool compare_vectors(const std::string& prefix, int myProc, const std::vector<stk::mesh::EntityKey>& gold, const std::vector<stk::mesh::EntityKey>& vec)
{
    if (gold.size() != vec.size()) {
        std::ostringstream os;
        os<<prefix<<" gold size="<<gold.size()<<", size="<<vec.size()<<std::endl;
        for(const stk::mesh::EntityKey& key : vec) {
            os<<"P"<<myProc<<" "<<prefix<<" val: "<<key<<std::endl;
        }
        std::cerr<<os.str();
        return false;
    }

    bool returnVal = true;
    std::ostringstream os;
    for(size_t i=0; i<gold.size(); ++i) {
        if (gold[i] != vec[i]) {
            os<<prefix<<" "<<i<<"-th value mismatch: gold="<<gold[i]<<", val="<<vec[i]<<std::endl;
            returnVal = false;
        }
    }
    std::cerr<<os.str();
    return returnVal;
}

void fill_send_ghosts(const stk::mesh::BulkData& bulk, stk::mesh::EntityProcVec& sendGhosts)
{
    stk::mesh::Ghosting* myGhosting = find_ghosting_by_name(ghostingName, bulk);
    if (myGhosting != nullptr) {
        myGhosting->send_list(sendGhosts);
    }
}

void verify_ghosting_lists(const std::string& prefix, stk::mesh::BulkData& bulk,
                           const std::vector<std::vector<stk::mesh::EntityId> >& elemIdsToGhostPerSendProc,
                           const std::vector<std::vector<stk::mesh::EntityId> >& newElemIdsToGhostPerSendProc,
                           const std::vector<stk::mesh::EntityKey>& expected_recvGhostsToRemove)
{
    int myProc = bulk.parallel_rank(), otherProc = 1 - bulk.parallel_rank();
    if (myProc == 0) {
        std::cout<<prefix<<std::endl;
    }

    stk::mesh::EntityProcVec elemsToGhost;
    convert_elem_ids_to_elemprocs(bulk, otherProc, elemIdsToGhostPerSendProc[myProc], elemsToGhost);

    stk::mesh::EntityProcVec expected_elemsToGhost;
    convert_elem_ids_to_elemprocs(bulk, otherProc, newElemIdsToGhostPerSendProc[myProc], expected_elemsToGhost);

    std::vector<stk::mesh::EntityKey> recvGhostsToRemove;
    stk::mesh::EntityProcVec currentSendGhosts;
    fill_send_ghosts(bulk, currentSendGhosts);

    sierra::nalu::NonConformalManager::computePreciseGhostingLists(bulk, elemsToGhost, currentSendGhosts, recvGhostsToRemove);

    EXPECT_TRUE(compare_vectors(prefix+" elemsToGhost", bulk, expected_elemsToGhost, elemsToGhost));
    EXPECT_TRUE(compare_vectors(prefix+" recvGhostsToRemove", bulk.parallel_rank(), expected_recvGhostsToRemove, recvGhostsToRemove));

    ghostElems(bulk, elemsToGhost);
}

void print_entities_and_downward_connections(const stk::mesh::BulkData& bulk, const stk::mesh::EntityProcVec& entityProcs)
{
    std::ostringstream os;
    int p = bulk.parallel_rank();
    for(const stk::mesh::EntityProc& entityProc : entityProcs) {
        os <<"P"<<p<<" "<<bulk.entity_key(entityProc.first)<<", -> "<<std::endl;
        stk::mesh::EntityRank erank = bulk.entity_rank(entityProc.first);
        for(stk::mesh::EntityRank irank = stk::topology::NODE_RANK; irank < erank; ++irank) {
            unsigned num = bulk.num_connectivity(entityProc.first, irank);
            const stk::mesh::Entity* ents = bulk.begin(entityProc.first, irank);
            os<<"    {";
            for(unsigned i=0; i<num; ++i) {
                os << bulk.entity_key(ents[i])<<" ";
            }
            os<<"}"<<std::endl;
        }
    }
    std::cerr<<os.str();
}

using stk::mesh::EntityKey;

void check_computePreciseGhostingLists(stk::mesh::BulkData& bulk)
{
    int myProc = bulk.parallel_rank();
    std::vector<std::vector<stk::mesh::EntityKey> > expected_recvGhostsToRemovePerProc = { {}, {} };

    verify_ghosting_lists("initial", bulk, {{2}, {3}}, {{2}, {3}}, expected_recvGhostsToRemovePerProc[myProc]);

    verify_ghosting_lists("add ghost", bulk, {{2, 1}, {3, 4}}, {{1}, {4}}, expected_recvGhostsToRemovePerProc[myProc]);

    expected_recvGhostsToRemovePerProc =
    {
      {EntityKey(stk::topology::NODE_RANK,17), EntityKey(stk::topology::NODE_RANK,18),
       EntityKey(stk::topology::NODE_RANK,19), EntityKey(stk::topology::NODE_RANK,20),
       EntityKey(stk::topology::ELEM_RANK,4)},
      { EntityKey(stk::topology::NODE_RANK,1), EntityKey(stk::topology::NODE_RANK,2),
       EntityKey(stk::topology::NODE_RANK,3), EntityKey(stk::topology::NODE_RANK,4),
       EntityKey(stk::topology::ELEM_RANK,1)}
    };

    verify_ghosting_lists("remove ghost", bulk, {{2}, {3}}, {{}, {}}, expected_recvGhostsToRemovePerProc[myProc]);
}

}//namespace

TEST(NonConformalGhosting, computePreciseGhostingLists)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(comm) == 2) {

        unsigned spatialDimension = 3;
        stk::mesh::MetaData meta(spatialDimension);
        stk::mesh::BulkData bulk(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);

        unit_test_utils::fill_mesh_hex8(bulk, "generated:1x1x4");

        check_computePreciseGhostingLists(bulk);
    }
}

