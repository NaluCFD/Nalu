/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef UNITTESTLINEARSYSTEM_H
#define UNITTESTLINEARSYSTEM_H

#include "LinearSystem.h"
#include "EquationSystem.h"

namespace unit_test_utils {

class TestLinearSystem : public sierra::nalu::LinearSystem
{
public:

 TestLinearSystem( sierra::nalu::Realm &realm, const unsigned numDof, sierra::nalu::EquationSystem *eqSys)
   : sierra::nalu::LinearSystem(realm, numDof, eqSys, nullptr), numSumIntoCalls_(0),
     ignoredLhs_(false)
  {}

  virtual ~TestLinearSystem() {}

  // Graph/Matrix Construction
  virtual void buildNodeGraph(const stk::mesh::PartVector & parts) {}
  virtual void buildFaceToNodeGraph(const stk::mesh::PartVector & parts) {}
  virtual void buildEdgeToNodeGraph(const stk::mesh::PartVector & parts) {}
  virtual void buildElemToNodeGraph(const stk::mesh::PartVector & parts) {}
  virtual void buildReducedElemToNodeGraph(const stk::mesh::PartVector & parts) {}
  virtual void buildFaceElemToNodeGraph(const stk::mesh::PartVector & parts) {}
  virtual void buildNonConformalNodeGraph(const stk::mesh::PartVector & parts) {}
  virtual void buildOversetNodeGraph(const stk::mesh::PartVector & parts) {}
  virtual void finalizeLinearSystem() {}

  // Matrix Assembly
  virtual void zeroSystem() {}
  virtual void zeroRhs() {}

  virtual void sumInto(
      unsigned numEntities,
      const stk::mesh::Entity* entities,
      const sierra::nalu::SharedMemView<const double*> & rhs,
      const sierra::nalu::SharedMemView<const double**> & lhs,
      const sierra::nalu::SharedMemView<int*> & localIds,
      const sierra::nalu::SharedMemView<int*> & sortPermutation,
      const char * trace_tag,
      bool ignoreLhs = false)
  {
    if (numSumIntoCalls_ == 0) {
      rhs_ = Kokkos::View<double*>("rhs_",rhs.dimension(0));
      for(size_t i=0; i<rhs.dimension(0); ++i) {
        rhs_(i) = rhs(i);
      }
      lhs_ = Kokkos::View<double**>("lhs_",lhs.dimension(0), lhs.dimension(1));
      if (!ignoreLhs) {
        for(size_t i=0; i<lhs.dimension(0); ++i) {
          for(size_t j=0; j<lhs.dimension(1); ++j) {
            lhs_(i,j) = lhs(i,j);
          }
        }
      }
    }
    if (ignoreLhs) {
      ignoredLhs_ = true;
    }
    Kokkos::atomic_add(&numSumIntoCalls_, 1u);
  }

  virtual void sumInto(
    const std::vector<stk::mesh::Entity> & sym_meshobj,
    std::vector<int> &scratchIds,
    std::vector<double> &scratchVals,
    const std::vector<double> & rhs,
    const std::vector<double> & lhs,
    const char *trace_tag=0,
    bool ignoreLhs = false)
  {
    if (ignoreLhs) {
      ignoredLhs_ = true;
    }
  }

  virtual void applyDirichletBCs(
    stk::mesh::FieldBase * solutionField,
    stk::mesh::FieldBase * bcValuesField,
    const stk::mesh::PartVector & parts,
    const unsigned beginPos,
    const unsigned endPos)
  {}

  virtual void prepareConstraints(
    const unsigned beginPos,
    const unsigned endPos)
  {}

  // Solve
  virtual int solve(stk::mesh::FieldBase * linearSolutionField) { return -1; }
  virtual void loadComplete() {}

  virtual void writeToFile(const char * filename, bool useOwned=true) {}
  virtual void writeSolutionToFile(const char * filename, bool useOwned=true) {}

  virtual void resetRows(
    std::vector<stk::mesh::Entity> nodeList,
    const unsigned beginPos,
    const unsigned endPos) {}

  unsigned numSumIntoCalls_;
  bool ignoredLhs_;
  Kokkos::View<double**> lhs_;
  Kokkos::View<double*> rhs_;

protected:
  virtual void beginLinearSystemConstruction() {}
  virtual void checkError(
    const int err_code,
    const char * msg) {}
};

}

#endif /* UNITTESTLINEARSYSTEM_H */

