/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <NaluEnv.h>

#include <mpi.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

#include <stk_util/environment/WallTime.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// NaluEnv - manage parallel and parallel output in Nalu
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
NaluEnv::NaluEnv()
  : parallelCommunicator_(MPI_COMM_WORLD),
    pSize_(-1),
    pRank_(-1),
    stdoutStream_(std::cout.rdbuf()),
    naluLogStream_(&std::cout), // std::cout redirects to log file
    naluParallelStream_(new std::ostream(&naluParallelStreamBuffer_)),
    parallelLog_(false),
    baseName_("nalu")
{
  // initialize
  MPI_Comm_size(parallelCommunicator_, &pSize_);
  MPI_Comm_rank(parallelCommunicator_, &pRank_);
}

//--------------------------------------------------------------------------
//-------- self ------------------------------------------------------------
//--------------------------------------------------------------------------
NaluEnv &NaluEnv::self()
{
  static NaluEnv s;
  return s;
}

//--------------------------------------------------------------------------
//-------- naluOutputP0 ----------------------------------------------------
//--------------------------------------------------------------------------
std::ostream &
NaluEnv::naluOutputP0()
{
  return *naluLogStream_;
}

//--------------------------------------------------------------------------
//-------- naluOutput ------------------------------------------------------
//--------------------------------------------------------------------------
std::ostream &
NaluEnv::naluOutput()
{
  return *naluParallelStream_;
}

//--------------------------------------------------------------------------
//-------- parallel_size ---------------------------------------------------
//--------------------------------------------------------------------------
int
NaluEnv::parallel_size()
{
  return pSize_;
}

//--------------------------------------------------------------------------
//-------- parallel_rank ---------------------------------------------------
//--------------------------------------------------------------------------
int
NaluEnv::parallel_rank()
{
  return pRank_;
}

//--------------------------------------------------------------------------
//-------- parallel_comm ---------------------------------------------------
//--------------------------------------------------------------------------
MPI_Comm
NaluEnv::parallel_comm()
{
  return parallelCommunicator_;
}

//--------------------------------------------------------------------------
//-------- set_log_file_stream ---------------------------------------------
//--------------------------------------------------------------------------
void
NaluEnv::set_log_file_stream(std::string naluLogName, bool pprint)
{
  if ( pRank_ == 0 ) {
    naluStreamBuffer_.open(naluLogName.c_str(), std::ios::out);
    naluLogStream_->rdbuf(&naluStreamBuffer_);
  }
  else {
    naluLogStream_->rdbuf(&naluEmptyStreamBuffer_);
  }

  // default to an empty stream buffer for parallel unless pprint is set
  parallelLog_ = pprint;
  if (parallelLog_) {
    int numPlaces = static_cast<int>(std::log10(pSize_-1)+1);

    std::stringstream paddedRank;
    paddedRank << std::setw(numPlaces) << std::setfill('0') << parallel_rank();

    // inputname.log -> inputname.log.16.02 for the rank 2 proc of a 16 proc job
    std::string parallelLogName =
        naluLogName + "." + std::to_string(pSize_) + "." + paddedRank.str();
    naluParallelStreamBuffer_.open(parallelLogName.c_str(), std::ios::out);
    naluParallelStream_->rdbuf(&naluParallelStreamBuffer_);
  }
  else {
    naluParallelStream_->rdbuf(stdoutStream_);
  }
}

//--------------------------------------------------------------------------
//-------- close_log_file_stream -------------------------------------------
//--------------------------------------------------------------------------
void
NaluEnv::close_log_file_stream()
{
  if ( pRank_ == 0 ) {
    naluStreamBuffer_.close();
  }
  if (parallelLog_) {
    naluParallelStreamBuffer_.close();
  }

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
NaluEnv::~NaluEnv()
{
  close_log_file_stream();
  delete naluParallelStream_;
}

//--------------------------------------------------------------------------
//-------- nalu_time -------------------------------------------------------
//--------------------------------------------------------------------------
double
NaluEnv::nalu_time()
{
  return stk::wall_time();
}

//--------------------------------------------------------------------------
//-------- set_base_name ---------------------------------------------------
//--------------------------------------------------------------------------
void
NaluEnv::set_base_name(std::string baseName)
{
  baseName_ = baseName;
}

//--------------------------------------------------------------------------
//-------- get_base_name ---------------------------------------------------
//--------------------------------------------------------------------------
std::string
NaluEnv::get_base_name()
{
  return baseName_;
}

} // namespace nalu
} // namespace Sierra
