/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MovingAveragePostProcessor_h
#define MovingAveragePostProcessor_h

#include <Algorithm.h>
#include <FieldTypeDef.h>

namespace stk { namespace mesh { class BulkData; } }

namespace sierra{
namespace nalu{

class TimeIntegrator;

class ExponentialMovingAverager
{
public:
  ExponentialMovingAverager(double timeScale);

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
    double timeScale);

  void execute();

  void add_fields(std::vector<std::string> fieldName);
  void add_parts_for_all_fields(stk::mesh::PartVector parts);
  void add_parts_for_field(std::string name, stk::mesh::PartVector parts);

  std::map<stk::mesh::FieldBase*, stk::mesh::FieldBase*>& get_field_map()
  {
    return fieldMap_;
  }

private:
  stk::mesh::BulkData& bulk_;
  TimeIntegrator& timeIntegrator_;

  std::map<stk::mesh::FieldBase*, stk::mesh::FieldBase*> fieldMap_;
  ExponentialMovingAverager averager_;

  std::map<std::string, stk::mesh::PartVector> partVecs_;

};

} // namespace nalu
} // namespace Sierra

#endif
