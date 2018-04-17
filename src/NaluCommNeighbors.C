// ***** TEMPORARY COPY OF STK FILE ********
// Will be delete after stk is updated with performance improvements in
// stk_util/parallel/CommBufferV.hpp and stk/parallel/CommNeighbors.[hpp,cpp]
// ***** TEMPORARY COPY OF STK FILE ********
//
// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stdlib.h>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>

#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <NaluCommNeighbors.hpp>

namespace nalu_stk {

//-----------------------------------------------------------------------

#if defined( STK_HAS_MPI )

void CommNeighbors::rank_error( const char * method , int p ) const
{
  std::ostringstream os ;
  os << "stk::CommNeighbors::" << method
     << "(" << p << ") ERROR: Not in [0:" << m_size << ")" ;
  throw std::range_error( os.str() );
}

//----------------------------------------------------------------------

CommNeighbors::~CommNeighbors()
{
  m_comm = stk::parallel_machine_null();
  m_size = 0 ;
  m_rank = 0 ;
}

stk::ParallelMachine CommNeighbors::setup_neighbor_comm(stk::ParallelMachine fullComm,
                                                        const std::vector<int>& sendProcs,
                                                        const std::vector<int>& recvProcs)
{
  MPI_Info info;
  MPI_Info_create(&info);
  int reorder = 0;
  const int* weights = (int*)MPI_UNWEIGHTED;
  stk::ParallelMachine neighborComm;
  MPI_Dist_graph_create_adjacent(fullComm,
              recvProcs.size(), recvProcs.data(), weights,
              sendProcs.size(), sendProcs.data(), weights,
              info, reorder, &neighborComm);
  MPI_Info_free(&info);
  return neighborComm;
}

CommNeighbors::CommNeighbors( stk::ParallelMachine comm, const std::vector<int>& neighbor_procs)
  : m_comm( comm ),
    m_size( stk::parallel_machine_size( comm ) ),
    m_rank( stk::parallel_machine_rank( comm ) ),
    m_send(),
    m_recv(),
    m_send_data(),
    m_recv_data(),
    m_send_procs(neighbor_procs),
    m_recv_procs(neighbor_procs)
{
  m_send.resize(m_size);
  m_recv.resize(m_size);
  stk::util::sort_and_unique(m_send_procs);
  stk::util::sort_and_unique(m_recv_procs);

  m_comm = setup_neighbor_comm(comm, m_send_procs, m_recv_procs);
}

CommNeighbors::CommNeighbors( stk::ParallelMachine comm, const std::vector<int>& send_procs, const std::vector<int>& recv_procs)
  : m_comm( comm ),
    m_size( stk::parallel_machine_size( comm ) ),
    m_rank( stk::parallel_machine_rank( comm ) ),
    m_send(),
    m_recv(),
    m_send_data(),
    m_recv_data(),
    m_send_procs(send_procs),
    m_recv_procs(recv_procs)
{
  m_send.resize(m_size);
  m_recv.resize(m_size);
  stk::util::sort_and_unique(m_send_procs);
  stk::util::sort_and_unique(m_recv_procs);

  m_comm = setup_neighbor_comm(comm, m_send_procs, m_recv_procs);
}

void setup_buffers(const std::vector<int>& procs,
                   const std::vector<CommBufferV>& data,
                   std::vector<int>& counts,
                   std::vector<int>& displs,
                   std::vector<unsigned char>& buf)
{
  counts.resize(procs.size());
  displs.resize(procs.size());

  int totalBytes = 0;
  for(size_t i=0; i<procs.size(); ++i) {
    int p = procs[i];
    counts[i] = data[p].size_in_bytes();
    displs[i] = totalBytes;
    totalBytes += counts[i];
  }

  buf.resize(totalBytes, 0);
  buf.clear();

  for(size_t i=0; i<procs.size(); ++i) {
    int p = procs[i];
    const unsigned char* rawbuf = data[p].raw_buffer();
    size_t len = counts[i];
    buf.insert(buf.end(), rawbuf, rawbuf+len);
  }
}

void store_recvd_data(const std::vector<unsigned char>& recvBuf,
                      const std::vector<int>& recvCounts,
                      const std::vector<int>& recvDispls,
                      const std::vector<int>& recvProcs,
                      std::vector<CommBufferV>& m_recv)
{
  for(size_t i=0; i<recvProcs.size(); ++i) {
    int p = recvProcs[i];
    const unsigned char* buf = recvBuf.data() + recvDispls[i];
    int len = recvCounts[i];
    m_recv[p].resize(len);
    std::memcpy(m_recv[p].raw_buffer(), buf, len);
  }
}

void CommNeighbors::perform_neighbor_communication(MPI_Comm neighborComm,
                                                   const std::vector<unsigned char>& sendBuf,
                                                   const std::vector<int>& sendCounts,
                                                   const std::vector<int>& sendDispls,
                                                         std::vector<unsigned char>& recvBuf,
                                                         std::vector<int>& recvCounts,
                                                         std::vector<int>& recvDispls)
{
  ThrowAssertMsg(sendCounts.size()==m_send_procs.size(), "Error, sendCounts should be same size as m_send_procs.");
  ThrowAssertMsg(sendDispls.size()==m_send_procs.size(), "Error, sendDispls should be same size as m_send_procs.");
  ThrowAssertMsg(recvCounts.size()==m_recv_procs.size(), "Error, recvCounts should be same size as m_recv_procs.");
  ThrowAssertMsg(recvDispls.size()==m_recv_procs.size(), "Error, recvDispls should be same size as m_recv_procs.");

  MPI_Neighbor_alltoall(sendCounts.data(), 1, MPI_INT,
                        recvCounts.data(), 1, MPI_INT, neighborComm);

  int totalRecv = 0;
  for(size_t i=0; i<recvCounts.size(); ++i) {
    recvDispls[i] = totalRecv;
    totalRecv += recvCounts[i];
  }
  recvBuf.resize(totalRecv);

  MPI_Neighbor_alltoallv(
      sendBuf.data(), sendCounts.data(), sendDispls.data(), MPI_BYTE,
      recvBuf.data(), recvCounts.data(), recvDispls.data(), MPI_BYTE, neighborComm);
}

void CommNeighbors::communicate()
{
  if (m_size == 1) {
    int len = m_send[0].size_in_bytes();
    const unsigned char* buf = m_send[0].raw_buffer();
    m_recv[0].resize(len);
    std::memcpy(m_recv[0].raw_buffer(), buf, len);
    return;
  }

  std::vector<int> sendCounts, recvCounts(m_recv_procs.size(), 0);
  std::vector<int> sendDispls, recvDispls(m_recv_procs.size(), 0);
  std::vector<unsigned char> sendBuf, recvBuf;

  setup_buffers(m_send_procs, m_send, sendCounts, sendDispls, sendBuf);

  perform_neighbor_communication(m_comm, sendBuf, sendCounts, sendDispls,
                                         recvBuf, recvCounts, recvDispls);

  store_recvd_data(recvBuf, recvCounts, recvDispls, m_recv_procs, m_recv);
}

#endif

}

