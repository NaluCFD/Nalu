/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Simulation.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <NaluParsing.h>
#include <Realms.h>
#include <xfer/Transfers.h>
#include <TimeIntegrator.h>
#include <LinearSolvers.h>
#include <UnitTests.h>

#include <stk_util/environment/Env.hpp>

#include <Ioss_SerializeIO.h>

// basic c++
#include <iostream>
#include <map>
#include <math.h>
#include <utility>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Simulation - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
bool Simulation::debug_ = false;

Simulation::Simulation(const YAML::Node& root_node) :
    m_root_node(root_node),
    timeIntegrator_(NULL),
    realms_(NULL),
    transfers_(NULL),
    linearSolvers_(NULL),
    unitTests_(NULL),
    serializedIOGroupSize_(0)
{}

Simulation::~Simulation() {
  delete realms_;
  delete transfers_;
  delete timeIntegrator_;
  delete linearSolvers_;
  if (unitTests_) delete unitTests_;
}

// Timers
//static
stk::diag::TimerSet &
Simulation::rootTimerSet()
{
  static stk::diag::TimerSet s_timerSet(sierra::Diag::TIMER_ALL);

  return s_timerSet;
}

//static
stk::diag::Timer& Simulation::rootTimer()
{
  static stk::diag::Timer s_timer = stk::diag::createRootTimer("Nalu", rootTimerSet());

  return s_timer;
}

//static
stk::diag::Timer& Simulation::outputTimer()
{
  static stk::diag::Timer s_timer("Output", rootTimer());
  return s_timer;
}


void Simulation::load(const YAML::Node & node)
{

  high_level_banner();

  if (0 != node.FindValue("UnitTests"))
    {
      sierra::Env::outputP0() << "\n\n Running Unit Tests \n\n" << std::endl;
      sierra::nalu::UnitTests unit_tests(*this);
      unit_tests.load(node);
    }

  // load the linear solver configs
  linearSolvers_ = new LinearSolvers(*this);
  linearSolvers_->load(node);

  // create the realms
  realms_ = new Realms(*this);
  realms_->load(node);

  // create the time integrator
  sierra::Env::outputP0() << std::endl;
  sierra::Env::outputP0() << "Time Integrator Review:  " << std::endl;
  sierra::Env::outputP0() << "=========================" << std::endl;
  timeIntegrator_ = new TimeIntegrator(*this);
  timeIntegrator_->load(node);

  // create the transfers; mesh is already loaded in realm
  sierra::Env::outputP0() << std::endl;
  sierra::Env::outputP0() << "Transfer Review:         " << std::endl;
  sierra::Env::outputP0() << "=========================" << std::endl;
  transfers_ = new Transfers(*this);
  transfers_->load(node);

}

void Simulation::setSerializedIOGroupSize(int siogs)
{
  if (siogs)
    {
      if (siogs < 0 || siogs > Env::parallel_size() || Env::parallel_size() % siogs != 0)
        {
          Env::outputP0() << "Error: Job requested serialized_io_group_size of " << siogs
                          << " which is incompatible with MPI size= " << Env::parallel_size()
                          << "... shutting down." << std::endl;
          Env::abort();
        }
      serializedIOGroupSize_ = siogs;
      Ioss::SerializeIO::setGroupFactor(siogs);
    }
}

void Simulation::breadboard()
{
  realms_->breadboard();
  timeIntegrator_->breadboard();
  transfers_->breadboard();
}
void Simulation::initialize()
{
  realms_->initialize();
  timeIntegrator_->initialize();
  transfers_->initialize();
}
void Simulation::run()
{
  Env::outputP0() << std::endl;
  Env::outputP0() << "*******************************************************" << std::endl;
  Env::outputP0() << "Simulation Shall Commence: number of processors= " << sierra::Env::parallel_size() << std::endl;
  Env::outputP0() << "*******************************************************" << std::endl;

  if (unitTests_)
    {
      unitTests_->run();
      if (runOnlyUnitTests_)
        return;
    }
  timeIntegrator_->integrate_realm();
}

void Simulation::high_level_banner() {

  Env::outputP0() << std::endl;
  Env::outputP0() << "=================================================================" << std::endl;
  Env::outputP0() << "                            Nalu:                                " << std::endl;
  Env::outputP0() << "      A low Mach number, turbulent reacting flow code with       " << std::endl;
  Env::outputP0() << "            coupling to PMR and object response (CHT)            " << std::endl;
  Env::outputP0() << "=================================================================" << std::endl;
  Env::outputP0() << std::endl;
  Env::outputP0() << " TPLs: Boost (BSL), STK (TBA), Trilinos (BSD) and YAML_cpp (MIT) " << std::endl;
  Env::outputP0() << std::endl;
  Env::outputP0() << "     Sandia National Laboratories, Albuquerque, New Mexico       " << std::endl;
  Env::outputP0() << "-----------------------------------------------------------------" << std::endl;
  Env::outputP0() << "       Notice: This computer software was prepared by            " << std::endl;
  Env::outputP0() << "       Sandia Corporation, hereinafter the Contractor            " << std::endl;
  Env::outputP0() << "       under Contract DE-AC04-94AL85000 with the                 " << std::endl;
  Env::outputP0() << "       Department of Energy (DOE).  All rights in the            " << std::endl;
  Env::outputP0() << "       computer software are reserved by DOE on behalf           " << std::endl;
  Env::outputP0() << "       of the United States Government and the                   " << std::endl;
  Env::outputP0() << "       Contractor as provided in the Contract. You are           " << std::endl;
  Env::outputP0() << "       authorized to use this computer software for              " << std::endl;
  Env::outputP0() << "       Governmental purposes but it is not to be                 " << std::endl;
  Env::outputP0() << "       released or distributed to the public. NEITHER            " << std::endl;
  Env::outputP0() << "       THE U.S.  GOVERNMENT NOR THE CONTRACTOR MAKES             " << std::endl;
  Env::outputP0() << "       ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY          " << std::endl;
  Env::outputP0() << "       LIABILITY FOR THE USE OF THIS SOFTWARE.                   " << std::endl;
  Env::outputP0() << "-----------------------------------------------------------------" << std::endl;
  Env::outputP0() << std::endl;

}
} // namespace nalu
} // namespace Sierra
