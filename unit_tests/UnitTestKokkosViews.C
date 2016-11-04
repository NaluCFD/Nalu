#include <gtest/gtest.h>
#include <limits>

#include "UnitTestUtils.h"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <Kokkos_Core.hpp>

#include <master_element/MasterElement.h>

namespace {

typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;

class Hex8Mesh : public ::testing::Test
{
protected:
    Hex8Mesh()
    : comm(MPI_COMM_WORLD), spatialDimension(3),
      meta(spatialDimension), bulk(meta, comm),
      topo(stk::topology::HEX_8),
      elemCentroidField(&meta.declare_field<VectorFieldType>(stk::topology::ELEM_RANK, "elemCentroid")),
      nodalPressureField(&meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "nodalPressure")),
      dpdxIpField(&meta.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdxIp"))
    {
      stk::mesh::put_field(*elemCentroidField, meta.universal_part(), spatialDimension, (double*)nullptr);
      double one = 1.0;
      stk::mesh::put_field(*nodalPressureField, meta.universal_part(), 1, &one);
      stk::mesh::put_field(*dpdxIpField, meta.universal_part(), spatialDimension, (double*)nullptr);
    }

    ~Hex8Mesh() {}

    void fill_mesh(const std::string& meshSpec = "generated:100x100x100")
    {
      unit_test_utils::fill_hex8_mesh(meshSpec, bulk);
    }

    stk::ParallelMachine comm;
    unsigned spatialDimension;
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk;
    stk::topology topo;
    VectorFieldType* elemCentroidField;
    ScalarFieldType* nodalPressureField;
    VectorFieldType* dpdxIpField;
};

template<class LOOP_BODY>
void bucket_loop(const stk::mesh::BucketVector& buckets, LOOP_BODY inner_loop_body)
{
    for(const stk::mesh::Bucket* bptr : buckets)
    {
        const stk::mesh::Bucket& bkt = *bptr;
        for(size_t j=0; j<bkt.size(); ++j)
        {
            inner_loop_body(bkt[j]);
        }
    }
}

template<class LOOP_BODY>
void kokkos_bucket_loop(const stk::mesh::BucketVector& buckets, LOOP_BODY inner_loop_body)
{
    Kokkos::parallel_for(buckets.size(), [&](const size_t& i)
    {
        const stk::mesh::Bucket& bkt = *buckets[i];
        for(size_t j=0; j<bkt.size(); ++j)
        {
            inner_loop_body(bkt[j]);
        }
    });
}

template<class LOOP_BODY>
void kokkos_thread_team_bucket_loop(const stk::mesh::BucketVector& buckets, LOOP_BODY inner_loop_body)
{
    typedef Kokkos::Schedule<Kokkos::Dynamic> DynamicScheduleType;
    typedef typename Kokkos::TeamPolicy<typename Kokkos::DefaultExecutionSpace, DynamicScheduleType>::member_type TeamHandleType;

    Kokkos::parallel_for(Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(buckets.size(), Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
    {
        const stk::mesh::Bucket& bkt = *buckets[team.league_rank()];
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bkt.size()), [&](const size_t& j)
        {
            inner_loop_body(bkt[j]);
        });
    });
}

TEST_F(Hex8Mesh, indexing_raw_arrays)
{
    fill_mesh();
    const VectorFieldType* coordField = static_cast<const VectorFieldType*>(meta.coordinate_field());
    EXPECT_TRUE(coordField != nullptr);

    sierra::nalu::HexSCS meSCS;
    double scs_error = 0.0;

    const int nodesPerElem = topo.num_nodes();
    const int numScsIp = meSCS.numIntPoints_;
    const int nDim = spatialDimension;

    std::vector<double> elemNodeCoords(nodesPerElem*nDim, 0);
    std::vector<double> elemNodePressures(nodesPerElem, 0);

    double* p_elemNodeCoords = elemNodeCoords.data();
    double* p_elemNodePressures = elemNodePressures.data();

    std::vector<double> scs_areav(numScsIp*nDim, 0);
    std::vector<double> dndx(nDim*numScsIp*nodesPerElem, 0);
    std::vector<double> deriv(nDim*numScsIp*nodesPerElem, 0);
    std::vector<double> det_j(numScsIp, 0);

    double* p_scs_areav = scs_areav.data();
    double* p_dndx = dndx.data();
    double* p_deriv = deriv.data();
    double* p_det_j = det_j.data();

    std::vector<double*> p_dpdxIp(nodesPerElem, nullptr);

    const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());

    kokkos_thread_team_bucket_loop(elemBuckets, [&](stk::mesh::Entity elem)
    {
        const stk::mesh::Entity* elemNodes = bulk.begin_nodes(elem);

        for(int n=0; n<nodesPerElem; ++n) {
            const double* nodeCoords = stk::mesh::field_data(*coordField, elemNodes[n]);
            p_dpdxIp[n] = stk::mesh::field_data(*dpdxIpField, elemNodes[n]);

            const int nodeOffset = n*nDim;
            for(int d=0; d<nDim; ++d) {
                p_elemNodeCoords[nodeOffset+d] = nodeCoords[d];
            }
            const double* nodePressure = stk::mesh::field_data(*nodalPressureField, elemNodes[n]);
            p_elemNodePressures[n] = nodePressure[0];
        }

        meSCS.determinant(1, p_elemNodeCoords, p_scs_areav, &scs_error);
        meSCS.grad_op(1, p_elemNodeCoords, p_dndx, p_deriv, p_det_j, &scs_error);

        for (int ip = 0; ip < numScsIp; ++ip ) {
            const int ipOffset = nDim*nodesPerElem*ip;
            for ( int ic = 0; ic < nodesPerElem; ++ic) {
                const int offSetDnDx = ipOffset + ic*nDim;
                for ( int j = 0; j < nDim; ++j ) {
                    p_dpdxIp[ic][j] += p_dndx[offSetDnDx+j]*p_elemNodePressures[ic]*p_scs_areav[ip*nDim+j];
                }
            }
        }
    });
}

TEST_F(Hex8Mesh, indexing_views)
{
    fill_mesh();
    const VectorFieldType* coordField = static_cast<const VectorFieldType*>(meta.coordinate_field());
    EXPECT_TRUE(coordField != nullptr);

    sierra::nalu::HexSCS meSCS;
    double scs_error = 0.0;

    const int nodesPerElem = topo.num_nodes();
    const int numScsIp = meSCS.numIntPoints_;
    const int nDim = spatialDimension;

    Kokkos::View<double**> elemNodeCoords("coords", nodesPerElem, nDim);
    Kokkos::View<double*> elemNodePressures("pressures", nodesPerElem);

    Kokkos::View<double**> scs_areav("areav", numScsIp, nDim);
    Kokkos::View<double***> dndx("dndx", nDim, nodesPerElem, numScsIp);
    std::vector<double> deriv(nDim*numScsIp*nodesPerElem, 0);
    std::vector<double> det_j(numScsIp, 0);

    double* p_deriv = deriv.data();
    double* p_det_j = det_j.data();

    std::vector<double*> p_dpdxIp(nodesPerElem, nullptr);

    const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());

    kokkos_thread_team_bucket_loop(elemBuckets, [&](stk::mesh::Entity elem)
    {
        const stk::mesh::Entity* elemNodes = bulk.begin_nodes(elem);

        for(int n=0; n<nodesPerElem; ++n) {
            const double* nodeCoords = stk::mesh::field_data(*coordField, elemNodes[n]);
            p_dpdxIp[n] = stk::mesh::field_data(*dpdxIpField, elemNodes[n]);

            for(int d=0; d<nDim; ++d) {
                elemNodeCoords(n,d) = nodeCoords[d];
            }
            const double* nodePressure = stk::mesh::field_data(*nodalPressureField, elemNodes[n]);
            elemNodePressures[n] = nodePressure[0];
        }

        meSCS.determinant(1, &elemNodeCoords(0,0), &scs_areav(0,0), &scs_error);
        meSCS.grad_op(1, &elemNodeCoords(0,0), &dndx(0,0,0), p_deriv, p_det_j, &scs_error);

        for (int ip = 0; ip < numScsIp; ++ip ) {
            for ( int ic = 0; ic < nodesPerElem; ++ic) {
                for ( int j = 0; j < nDim; ++j ) {
                    p_dpdxIp[ic][j] += dndx(j,ic,ip)*elemNodePressures(ic)*scs_areav(ip,j);
                }
            }
        }
    });
}

}

