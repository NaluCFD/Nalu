/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PointInfo_h
#define PointInfo_h

//==============================================================================
// Includes and forwards
//==============================================================================

#include <master_element/MasterElement.h>

// stk
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Ghosting.hpp>

// stk_search
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <vector>
#include <map>

namespace sierra {
namespace nalu {

typedef stk::search::IdentProc<uint64_t,int>  uint64IdentProc;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<double> Box;
typedef std::pair<Point,uint64IdentProc> boundingPoint;
typedef std::pair<Sphere,uint64IdentProc> boundingSphere;
typedef std::pair<Box,uint64IdentProc> boundingBox;

struct NaluOneDMatrixSystem {
  int globalId_;
  double Aw_, Ap_, Ae_;
  double rhs_;
  double dof_;
  NaluOneDMatrixSystem(unsigned int globalId)
    : globalId_(globalId), Aw_(0.0), Ap_(0.0), Ae_(0.0), rhs_(0.0), dof_(0.0) {}

  void zero() { Aw_ = 0.0; Ap_ = 0.0; Ae_ = 0.0; rhs_ = 0.0; dof_ = 0.0; }
};

//=============================================================================
// Class Definition
//=============================================================================
// PointInfo
//=============================================================================
/**
 * * @par Description:
 * - class to manage pointinformation.
 *
 * @par Design Considerations:
 * -
 */
//=============================================================================
class PointInfo {

 public:

  // constructor and destructor
  PointInfo(
    boundingPoint bPoint,
    const uint64_t localPointId,
    Point &ipCoords,
    Point &pointCoords,
    const int nDim,
    const bool odeActive);
  ~PointInfo();
  
  // if this is an ODE-based approach, we need to initialize the points
  void initialize_ode_points();
  
  // helper solver functions
  void scatter(
    NaluOneDMatrixSystem *mL, NaluOneDMatrixSystem *mR, const double *A, const double *R);

  void dirichlet(
    NaluOneDMatrixSystem *mS, double du);

  void compute_visc(
    const double &tauWall, const double &rhoWall, const double &muWall, 
    const double kappa, std::vector<double> &coordsVec, std::vector<double> &viscVec);

  void assemble( 
    std::vector<NaluOneDMatrixSystem> &AmatrixSystemVec,
    std::vector<std::pair<int,int> > &lrscv, 
    std::vector<double> &coordsVec, 
    std::vector<double> &velocityVec, 
    std::vector<double> &viscVec);
 
  void tdma(
    int mSize, NaluOneDMatrixSystem *ns, double *u);
  
  double solve(
    const double &uExchange, const double &uWall, 
    const double &rhoWall, const double &muWall, 
    const double &kappa, const double &tauWallProvided,
    bool &converged);

  boundingPoint bPoint_;
  // should be able to extract this below from bPoint, right?
  const uint64_t localPointId_;
  const Point ipCoordinates_;
  const Point pointCoordinates_;
  const int nDim_;
  const bool odeActive_;

  stk::mesh::Entity owningElement_;

  double bestX_;
  const double bestXRef_;

  int elemIsGhosted_;

  // master element for background mesh
  MasterElement *meSCS_;

  std::vector<double> isoParCoords_;

  // ODE-based specifications
  double odeFac_;
  const int numPoints_;
  const double bias_;
  double tauWall_;
  int iterations_;
  bool debugOutput_;

  // data structures for ODE-based approach
  std::vector<double> coordsVec_;
  std::vector<double> velocityVec_;
  std::vector<double> viscVec_;
  std::vector<std::pair<int,int> >lrscv_;
  std::vector<NaluOneDMatrixSystem> AmatrixSystemVec_;
};

} // end sierra namespace
} // end Acon namespace

#endif
