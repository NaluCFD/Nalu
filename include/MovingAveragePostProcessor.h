/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MovingAveragePostProcessor_h
#define MovingAveragePostProcessor_h

#include <FieldTypeDef.h>

#include <map>
#include <vector>

namespace stk { namespace mesh { class BulkData; } }

namespace sierra{
namespace nalu{

class TimeIntegrator;



class ExponentialMovingAverager
{
public:
  ExponentialMovingAverager(double timeScale = 0.0, bool isInit = false, double alpha = -1.);

  double compute_updated_average(double oldAvg, double newVal);
  void compute_and_set_alpha(double dt);
  void init_state(bool init);
private:
  double timeScale_;
  bool isInit_;
  double alpha_;
};


class MovingAveragePostProcessor
{
public:

  // Field naming rule
  static std::string filtered_field_name(std::string unfilteredFieldName)
  {
    // "_ma" for symmetry with the turbulence averaging naming rule
    return (unfilteredFieldName + "_ma");
  }

  MovingAveragePostProcessor(
    stk::mesh::BulkData& bulk,
    TimeIntegrator& timeIntegrator,
    bool init = true);

  void execute();

  void add_fields(std::vector<std::string> fieldName);
  void set_time_scale(std::string fieldName, double timeScale);
 
  std::map<stk::mesh::FieldBase*, stk::mesh::FieldBase*>& get_field_map()
  {
    return fieldMap_;
  }

private:
  stk::mesh::BulkData& bulk_;
  TimeIntegrator& timeIntegrator_;
  bool isRestarted_;
  std::map<stk::mesh::FieldBase*, stk::mesh::FieldBase*> fieldMap_;
  std::map<std::string, ExponentialMovingAverager> averagers_;
};

} // namespace nalu
} // namespace Sierra

#endif
