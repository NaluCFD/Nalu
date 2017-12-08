
#ifndef ABLFORCINGALGORITHM_H
#define ABLFORCINGALGORITHM_H

#include "NaluParsing.h"
#include "FieldTypeDef.h"

#include "stk_mesh/base/Selector.hpp"

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <unordered_set>

namespace sierra {
namespace nalu {

class Realm;
class Transfer;
class Transfers;

/**
 * \brief ABL Forcing Source terms for Momentum and Temperature equations
 *
 * This class parses the user inputs and provides a common interface to the
 * momentum and temperature ABL forcing source term implementations within Nalu.
 * The ABL forcing capability is turned on by the presence of a sub-section
 * titled `abl_forcing` within the Realm section of the Nalu input file.
 *
 * ```
 *   abl_forcing:
 *     search_method: stk_octree
 *     search_tolerance: 0.0001
 *     search_expansion_factor: 1.5
 *
 *     from_target_part: [block_101]
 *
 *     momentum:
 *       type: computed
 *       relaxation_factor: 1.0
 *       heights: [80.0]
 *       target_part_format: "abl_height_%.1f"
 *       velocity_x:
 *         - [0.0, 10.0]                 # [Time0, vxH0, ... , vxHN]
 *         - [100000.0, 10.0]            # [TimeN, vxH0, ... , vxHN]
 *
 *       velocity_y:
 *         - [0.0, 0.0]
 *         - [10000.0, 0.0]
 *
 *       velocity_z:
 *         - [0.0, 0.0]
 *         - [10000.0, 0.0]
 *
 *     temperature:
 *       type: computed
 *       relaxation_factor: 1.0
 *       heights: [80.0]
 *       target_part_format: "abl_height_%.1f"
 *       temperature:
 *         - [0.0, 300.0]
 *         - [10000.0, 300.0]
 * ```
 *
 * There are two optional sub-sections in `abl_forcing`: "momentum" and
 * "temperature".
 */
class ABLForcingAlgorithm
{
public:
  template <typename T>
  using Array2D = std::vector<std::vector<T>>;

  /**
   * Types of ABL forcing available
   */
  enum ABLForcingTypes {
    OFF = 0,              //!< No ABL forcing applied
    USER_DEFINED = 1,     //!< Source terms provided by user
    COMPUTED = 2,         //!< Forcing computed by code given target profiles
    NUM_ABL_FORCING_TYPES //!< Guard
  };

  ABLForcingAlgorithm(Realm&, const YAML::Node&);

  ~ABLForcingAlgorithm();

  //! Parse input file for user options and initialize
  void load(const YAML::Node&);

  //! Setup ABL forcing (steps before mesh creation)
  void setup();

  //! Initialize ABL forcing (steps after mesh creation)
  void initialize();

  //! Execute field transfers, compute planar averaging, and determine source
  //! terms at desired levels.
  void execute();

  //! Evaluate the ABL forcing source contribution at a node
  void eval_momentum_source(
    const double,        //!< Height of the node from terrain
    std::vector<double>& //!< Source vector to be populated
    );

  //! Evaluate the ABL forcing source contribution (temperature)
  void eval_temperature_source(
    const double, //!< Height of the node from terrain
    double&       //!< Temperature source term to be populated
    );

  //! Inactive selector representing union of all the parts
  inline stk::mesh::Selector& inactive_selector() { return inactiveSelector_; }

  inline bool momentumForcingOn() { return (momSrcType_ != OFF); }

  inline bool temperatureForcingOn() { return (tempSrcType_ != OFF); }

  inline bool ablForcingOn()
  {
    return (momentumForcingOn() || temperatureForcingOn());
  }

private:
  ABLForcingAlgorithm();
  ABLForcingAlgorithm(const ABLForcingAlgorithm&);

  //! Utility function to parse momentum forcing options from input file.
  void load_momentum_info(const YAML::Node&);

  //! Helper method to parse temperature forcing options from input file.
  void load_temperature_info(const YAML::Node&);

  //! Create 2-D interpolation lookup tables from YAML data structures
  void create_interp_arrays(
    const std::vector<double>::size_type,
    const Array2D<double>&,
    std::vector<double>&,
    Array2D<double>&);

  //! Helper method that determines the part corresponding to a desired
  //! vertical level and ensures that part exists in the mesh database.
  void determine_part_names(
    std::vector<double>&, std::vector<std::string>&, bool, std::string&);

  //! Register velocity and temperature fields on the appropriate parts based
  //! on user input.
  void register_fields();

  //! Create transfer that handles mapping of velocity and temperature from
  //! fluidRealm to the planar nodesets.
  void create_transfers();

  void populate_transfer_data(std::string, const std::vector<std::string>&);

  //! Compute mean velocity and estimate source term for a given timestep
  void compute_momentum_sources();

  //! Compute average planar temperature and estimate source term
  void compute_temperature_sources();

  //! Helper method to compute the average velocity on a z-planes
  void calc_mean_velocity();

  void calc_mean_temperature();

  //! Reference to Realm
  Realm& realm_;

  //! Momentum Forcing Source Type
  ABLForcingTypes momSrcType_;

  //! Temperature Forcing Source Type
  ABLForcingTypes tempSrcType_;

  //! Relaxation factor for momentum sources
  double alphaMomentum_;

  //! Relaxation factor for temperature sources
  double alphaTemperature_;

  //! Heights where velocity information is provided
  std::vector<double> velHeights_; // Array of shape [num_Uheights]

  //! Heights where temperature information is provided
  std::vector<double> tempHeights_; // Array of shape [num_Theights]

  //! Times where velocity information is available
  std::vector<double> velXTimes_; // Arrays of shape [num_Utimes]
  std::vector<double> velYTimes_;
  std::vector<double> velZTimes_;

  //! Times where temperature information is available
  std::vector<double> tempTimes_; // Array of shape [num_Ttimes]

  // The following arrays are shaped [num_UHeights, num_Utimes]
  Array2D<double> velX_;
  Array2D<double> velY_;
  Array2D<double> velZ_;
  // The temperature array is shaped [num_Theights, num_Ttimes]
  Array2D<double> temp_;

  //! Planar average velocity calculated on the surface [num_UHeights, 3]
  Array2D<double> UmeanCalc_;

  //! Planar average density calculated on the surface [num_UHeights]
  std::vector<double> rhoMeanCalc_;

  //! U source as a function of height [3,num_UHeights]
  Array2D<double> USource_;

  //! Planar average temperature calculated on the surface [num_THeights]
  std::vector<double> TmeanCalc_;

  //! T source as a function of height [num_THeights]
  std::vector<double> TSource_;

  //! stk::Transfer search methods
  std::string searchMethod_;
  //! stk::Transfer search tolerance
  double searchTolerance_;
  //! stk::Transfer search expansion factor
  double searchExpansionFactor_;

  //! Domains where velocity/temperature are averaged
  std::vector<std::string> fromTargetNames_;

  //! Part names
  std::vector<std::string> velPartNames_;
  std::vector<std::string> tempPartNames_;
  std::unordered_set<std::string> allPartNames_;

  stk::mesh::PartVector allParts_;
  stk::mesh::Selector inactiveSelector_;

  Transfers* transfers_;

  //! Flag indicating whether to generate part names list for velocity field
  bool velGenPartList_;

  //! Flag indicating whether to generate part names list for temperature field
  bool tempGenPartList_;

  //! Format string specifier for generating velocity parts list
  std::string velPartFmt_;

  //! Format string specifier for generating temperature parts list
  std::string tempPartFmt_;

  //! Write frequency for source term output
  int outputFreq_;

  //! Format string specifier indicating the file name for output. The
  //! specification takes one `%s` specifier that is used to populate Ux, Uy,
  //! Uz, T. Default is "abl_sources_%s.dat"
  std::string outFileFmt_;
};
}
}

#endif /* ABLFORCINGALGORITHM_H */
