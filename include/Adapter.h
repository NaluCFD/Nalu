/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_PERCEPT)

#ifndef Adapter_h
#define Adapter_h

#include <stk_mesh/base/Selector.hpp>

namespace percept {
  class PerceptMesh;
  class AdaptedMeshVerifier;
  struct ElementRefinePredicate;
  template<class RefinePredicate> class TransitionElementAdapter;
  class UniformRefinerPatternBase;
  class UniformRefiner;
}

namespace sierra{
namespace nalu{

class Realm;

enum AdapterInsruction {
  ADAPT_NOTHING = 0,
  ADAPT_REFINE = 1 << 1,
  ADAPT_UNREFINE = 1 << 2
};

class Adapter
{
public:

  Adapter(
          const Realm &realm);
  ~Adapter();

  void do_adapt(int what_to_do);
  void do_uniform_refine();

  const Realm& realm_;
  percept::UniformRefinerPatternBase *uniformRefinementPattern_;
  percept::UniformRefinerPatternBase *refinementPattern_;
  percept::PerceptMesh *perceptMesh_;
  percept::UniformRefiner *uniformBreaker_;
  percept::TransitionElementAdapter<percept::ElementRefinePredicate> *breaker_;
  percept::ElementRefinePredicate *elementRefinePredicate_;
  stk::mesh::Selector *selector_;
  // mesh verifier
  percept::AdaptedMeshVerifier * adaptedMeshVerifier_;

 private:
  void setNaluGlobalId();
};

} // namespace nalu
} // namespace Sierra

#endif

#endif
