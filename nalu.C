/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <Slib_Startup.h>
#include <stk_util/environment/Env.hpp>
#include <stk_util/diag/PrintTimer.hpp>

// nalu
#include <NaluParsing.h>
#include <Simulation.h>

// util
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

// boost for input params
#include <boost/program_options.hpp>

// yaml for parsing..
#include <yaml-cpp/yaml.h>

#include <iostream>
#include <iomanip>
#include <stdexcept>

static std::string human_bytes_double(double bytes)
{
  const double K = 1024;
  const double M = K*1024;
  const double G = M*1024;

  std::ostringstream out;
  if (bytes < K) {
    out << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << bytes << " K";
  } else if (bytes < G) {
    bytes /= M;
    out << bytes << " M";
  } else {
    bytes /= G;
    out << bytes << " G";
  }
  return out.str();
}


int main( int argc, char ** argv )
{
  stk::diag::setEnabledTimerMetricsMask(stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME);

  sierra::nalu::Simulation::rootTimer().start();

  // start initial time
  double start_time = stk::cpu_time();

  std::string inputFileName;
  bool debug = false;
  int serializedIOGroupSize = 0;

  // Add my command line options to the option descriptions.
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h","produce help message")
    ("input-deck,i", boost::program_options::value<std::string>(&inputFileName)->default_value("myInput.i"),
     "Analysis input file")
    ("serialized-io-group-size,s",
     boost::program_options::value<int>(&serializedIOGroupSize)->default_value(0),
     "Specifies the number of processors which can concurrently perform I/O. Specifying zero disables serialization.")
    ("debug,D", "debug print on");
  //    ("debug,D", boost::program_options::value<bool>(&debug)->default_value(false),"debug print on");

  stk::get_options_description().add(desc);

  sierra::Env::Startup startup__(&argc, &argv, "nalu", __DATE__ " " __TIME__); //, opts);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, stk::get_options_description()), vm);

  boost::program_options::notify(vm);

  if (vm.count("version"))
    {
      if (!sierra::Env::parallel_rank())
        std::cerr << "version = Nalu1.0" << std::endl;
      return 0;
    }

  if (vm.count("debug"))
    {
      debug = true;
    }

  // deal with YAML; for now only one doc
  std::ifstream fin(inputFileName.c_str());
  if (!fin.good())
    {
      if (!sierra::Env::parallel_rank())
        std::cerr << "Input file is not specified or does not exist: user specified (or default) name= " << inputFileName << std::endl;
      return 0;
    }


  YAML::Parser parser(fin);
  YAML::Node doc;

  try {
    parser.GetNextDocument(doc);
    if (debug)
      {
        //sierra::nalu::NaluParsingHelper::traverse(std::cout, doc);
        if (!sierra::Env::parallel_rank())
          sierra::nalu::NaluParsingHelper::emit(std::cout, doc);
      }
  }
  catch (YAML::ParserException &e) {
    std::cout << e.what() << std::endl;
  }

  sierra::nalu::Simulation sim(doc);
  if (serializedIOGroupSize)
    {
      sierra::Env::outputP0() << "Info: found non-zero serialized_io_group_size on command-line= "
                              << serializedIOGroupSize << " (takes precedence over input file value)."
                              << std::endl;
      sim.setSerializedIOGroupSize(serializedIOGroupSize);
    }
  sim.debug_ = debug;
  sim.load(doc);
  sim.breadboard();
  sim.initialize();
  sim.run();

  // stop timer
  const double stop_time = stk::cpu_time();
  const double total_time = stop_time - start_time;
  const char* timer_name = "Total Time";

  // parallel reduce overall times
  double g_sum, g_min, g_max;
  stk::all_reduce_min(sierra::Env::parallel_comm(), &total_time, &g_min, 1);
  stk::all_reduce_max(sierra::Env::parallel_comm(), &total_time, &g_max, 1);
  stk::all_reduce_sum(sierra::Env::parallel_comm(), &total_time, &g_sum, 1);
  const int nprocs = sierra::Env::parallel_size();

  // output total time
  sierra::Env::outputP0() << "Timing for Simulation: nprocs= " << nprocs << std::endl;
  sierra::Env::outputP0() << "           main() --  " << " \tavg: " << g_sum/double(nprocs)
                          << " \tmin: " << g_min << " \tmax: " << g_max << std::endl;

  // output memory usage
  {
    size_t now, hwm;
    stk::get_memory_usage(now, hwm);
    // min, max, sum
    size_t global_now[3] = {now,now,now};
    size_t global_hwm[3] = {hwm,hwm,hwm};

    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceSum<1>( &global_now[2] ) );
    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceMin<1>( &global_now[0] ) );
    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceMax<1>( &global_now[1] ) );

    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceSum<1>( &global_hwm[2] ) );
    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceMin<1>( &global_hwm[0] ) );
    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceMax<1>( &global_hwm[1] ) );

    sierra::Env::outputP0() << "Memory Overview: " << std::endl;

    sierra::Env::outputP0() << "nalu memory: total (over all cores) current/high-water mark= "
                            << std::setw(15) << human_bytes_double(global_now[2])
                            << std::setw(15) << human_bytes_double(global_hwm[2])
                            << std::endl;

    sierra::Env::outputP0() << "nalu memory:   min (over all cores) current/high-water mark= "
                            << std::setw(15) << human_bytes_double(global_now[0])
                            << std::setw(15) << human_bytes_double(global_hwm[0])
                            << std::endl;

    sierra::Env::outputP0() << "nalu memory:   max (over all cores) current/high-water mark= "
                            << std::setw(15) << human_bytes_double(global_now[1])
                            << std::setw(15) << human_bytes_double(global_hwm[1])
                            << std::endl;
  }

  sierra::nalu::Simulation::rootTimer().stop();

  //output timings consistent w/ rest of Sierra
  stk::diag::Timer & sierra_timer = sierra::nalu::Simulation::rootTimer();
  const double elapsed_time = sierra_timer.getMetric<stk::diag::CPUTime>().getAccumulatedLap(false);
  stk::diag::Timer & mesh_output_timer = sierra::nalu::Simulation::outputTimer();
  double mesh_output_time = mesh_output_timer.getMetric<stk::diag::CPUTime>().getAccumulatedLap(false);
  double time_without_output = elapsed_time-mesh_output_time;

  stk::parallel_print_time_without_output_and_hwm(sierra::Env::parallel_comm(), time_without_output, sierra::Env::outputP0());

  if (!sierra::Env::parallel_rank())
    stk::print_timers_and_memory(&timer_name, &total_time, 1 /*num timers*/);

  stk::diag::printTimersTable(sierra::Env::outputP0(), sierra::nalu::Simulation::rootTimer(),
                              stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME,
                              false, sierra::Env::parallel_comm());

  return 0;
}
