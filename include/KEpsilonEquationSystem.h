/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef KEpsilonEquationSystem_h
#define KEpsilonEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

namespace stk{
struct topology;
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class EquationSystems;
class AlgorithmDriver;
class TurbKineticEnergyEquationSystem;
class TurbDissipationEquationSystem;

class KEpsilonEquationSystem : public EquationSystem {

public:

  KEpsilonEquationSystem(
    EquationSystems& equationSystems,
    const bool outputClippingDiag);
  virtual ~KEpsilonEquationSystem();
  
  virtual void initialize();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void solve_and_update();

  void initial_work();
  void update_and_clip();

  const bool outputClippingDiag_;

  TurbKineticEnergyEquationSystem *tkeEqSys_;
  TurbDissipationEquationSystem *epsEqSys_;

  ScalarFieldType *tke_;
  ScalarFieldType *eps_;

  bool isInit_;     
};

} // namespace nalu
} // namespace Sierra

#endif
