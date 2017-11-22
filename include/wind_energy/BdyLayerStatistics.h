/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef BDYLAYERSTATISTICS_H
#define BDYLAYERSTATISTICS_H

#include "NaluParsing.h"
#include "FieldTypeDef.h"

#include "stk_mesh/base/Part.hpp"

#include <memory>

namespace sierra {
namespace nalu {

class Realm;
class TurbulenceAveragingPostProcessing;
class AveragingInfo;
class BdyHeightAlgorithm;

/** Boundary layer statistics post-processing utility
 *
 *  A post-processing utility to compute statistics of flow quantities that are
 *  temporally and spatially averaged. Primarily used to evaluate atmospheric
 *  boundary layer characteristics during precursor simulations, this utility
 *  can also be used for channel flows.
 *
 *  The temporal averaging is perfomed via
 *  sierra::nalu::TurbulenceAveragingPostProcessing class.
 */
class BdyLayerStatistics
{
public:
  BdyLayerStatistics(
    Realm&,
    const YAML::Node&);

  virtual ~BdyLayerStatistics();

  /** Check for parts and initialize indexing field
   */
  void setup();

  /** Determine the number of height levels and node association with height levels
   */
  void initialize();

  /** Compute statistics of interest and output them if necessary
   */
  void execute();

  /** Return the spatial average of the instantaneous velocity field at a given height
   *
   *  The method interpolates the spatially averaged velocity from available
   *  heights to return the average velocity at an arbitrary height.
   *
   *  @param[in] height The height where velocity is desired
   *  @param[out] velArray A pointer to an array of nDim which contains the components of velocity
   */
  void velocity(double, double*);

  /** Return the spatially and temporally averaged velocity at a given height
   *
   *  The method interpolates the spatially averaged velocity from available
   *  heights to return the average velocity at an arbitrary height.
   *
   *  @param[in] height The height where velocity is desired
   *  @param[out] velArray A pointer to an array of nDim which contains the components of velocity
   */
  void time_averaged_velocity(double, double*);

private:
  BdyLayerStatistics() = delete;
  BdyLayerStatistics(const BdyLayerStatistics&) = delete;

  //! Process the user inputs and initialize class data
  void load(const YAML::Node&);

  //! Initialize necessary parameters in sierra::nalu::TurbulenceAveragingPostProcessing
  void setup_turbulence_averaging(const double);

  //! Process the velocity data and compute averages
  void compute_velocity_stats();

  //! Process the temperature field and compute averages
  void compute_temperature_stats();

  //! Output averaged velocity and stress profiles as a function of height
  void output_velocity_averages();

  //! Output averaged temperature profiles as a function of height
  void output_temperature_averages();

  /** Helper method to perform interpolations across data with multiple components
   *
   *  @param[in] nComp The number of components: scalar=1, vector=3, tensor=6 and so on
   *  @param[in] varArray Reference to the array that contains averaged data at all heights
   *  @param[in] height The height where data is to be interpolated
   *  @param[out] interpArry Pointer to an array of size nComp where interpolated data is populated
   */
  void interpolate_variable(
    int,
    std::vector<double>&,
    double,
    double*);

  //! Reference to Realm object
  Realm& realm_;

  //! Spatially averaged instantaneous velocity at desired heights [nHeights, nDim]
  std::vector<double> velAvg_;

  //! Spatially and temporally averaged velocity at desired heights [nHeights, nDim]
  std::vector<double> velBarAvg_;

  //! Spatially and temporally averaged resolved stress field at desired heights [nHeights, nDim * 2]
  std::vector<double> uiujAvg_;

  //! Spatially and temporally averaged SFS field at desired heights [nHeights, nDim * 2]
  std::vector<double> sfsAvg_;

  //! Spatially averaged instantaneous temperature field [nHeights]
  std::vector<double> thetaAvg_;

  //! Spatially and temporally averaged temperature field [nHeights]
  std::vector<double> thetaBarAvg_;

  //! Spatially and temporally averaged Temperature SFS field
  std::vector<double> thetaSFSAvg_;

  std::vector<double> thetaUjAvg_;

  //! Total nodal volume at each height level used for volumetric averaging
  std::vector<double> sumVol_;

  //! Height from the wall
  std::vector<double> heights_;

  //! Part names for post-processing
  std::vector<std::string> partNames_;

  //! Parts of the fluid mesh where velocity/temperature averaging is performed
  stk::mesh::PartVector fluidParts_;

  //! Dimensionality of the mesh
  int nDim_{3};

  //! Output frequency for averaged statistics
  int outputFrequency_{10};

  //! Output frequency for time histories
  int timeHistOutFrequency_{10};

  //! Height index field
  ScalarIntFieldType* heightIndex_;

  std::unique_ptr<BdyHeightAlgorithm> bdyHeightAlg_;

  //! Calculate temperature statistics
  bool calcTemperatureStats_{true};

  //! Flag indicating whether initialization must be performed
  bool doInit_{true};
};

}  // nalu
}  // sierra


#endif /* BDYLAYERSTATISTICS_H */
