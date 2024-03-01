/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <PointInfo.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{


//==========================================================================
// Class Definition
//==========================================================================
// PointInfo - contains point -> donor elements mapping
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
PointInfo::PointInfo(
  boundingPoint bPoint,
  const uint64_t localPointId,
  Point &ipCoords,
  Point &pointCoords,
  const int nDim,
  const bool odeActive)
  : bPoint_(bPoint),
    localPointId_(localPointId),
    ipCoordinates_(ipCoords),
    pointCoordinates_(pointCoords),
    nDim_(nDim),
    odeActive_(odeActive),
    owningElement_(),
    bestX_(1.0e16),
    bestXRef_(1.0e16),
    elemIsGhosted_(0),
    meSCS_(NULL),
    odeFac_(0.0),
    numPoints_(50),
    bias_(200),
    tauWall_(-1.0),
    iterations_(0),
    debugOutput_(false)
{
  // manage ODE set of points
  if ( odeActive_ ) {
    odeFac_ = 1.0;
    initialize_ode_points();
  } 
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
PointInfo::~PointInfo()
{
  // nothing to delete
}

//--------------------------------------------------------------------------
//-------- initialize_ode_points -------------------------------------------
//--------------------------------------------------------------------------
void
PointInfo::initialize_ode_points()
{
  // derived number of segments (number of edges)
  const int n = numPoints_-1;
  const int nm1 = n - 1;

  // derived amplification, r = Li/Li-1 = bias^(1/nm1)
  const double r = std::pow(bias_, 1.0/(double)nm1);  

  // total length based on projected and ip point
  double Ltotal = 0.0;
  for (int i = 0; i < nDim_; ++i ) {
    Ltotal += (pointCoordinates_[i] - ipCoordinates_[i])*(pointCoordinates_[i] - ipCoordinates_[i]);
  }
  Ltotal = std::sqrt(Ltotal);

  // first distance
  const double L1 = Ltotal*(1.0-std::pow(bias_, 1.0/((double)nm1)))/(1.0-std::pow(bias_, (double)n/nm1));

  // create vector of points and deltas
  /*
    0    1     2     3     4 ===> coordsVec
    o----o-----o-----o-----o
      0     1     2     3 ======> dyVec
   */

  // populate a dy
  std::vector<double> dyVec;
  dyVec.push_back(L1);
  for ( int k = 1; k < n; ++k ) {
    const double dyNext = dyVec[k-1]*r;
    dyVec.push_back(dyNext);
  }
  
  // populate coords; assume a zero point location
  coordsVec_.push_back(0.0);
  double ccoords = 0.0;
  for ( size_t k = 0; k < dyVec.size(); ++k ) {
    ccoords += dyVec[k];
    coordsVec_.push_back(ccoords);
  }

  // initialize velocity
  for ( int k = 0; k < numPoints_; ++k ) {
    velocityVec_.push_back(0.0);
  }

  // initialize viscosity
  for ( int k = 0; k < numPoints_; ++k ) {
    viscVec_.push_back(0.0);
  }

  // start left/right mapping for each virtial edge
  for ( int k = 0; k < n; ++k ) {
    lrscv_.push_back(std::make_pair(k,k+1));
  }

  // create a matrix object for each node
  for ( int k = 0; k < numPoints_; ++k ) {
    AmatrixSystemVec_.push_back(NaluOneDMatrixSystem(k));
  }
}

//--------------------------------------------------------------------------
//-------- scatter ---------------------------------------------------------
//--------------------------------------------------------------------------
void 
PointInfo::scatter(
  NaluOneDMatrixSystem *mL, NaluOneDMatrixSystem *mR, const double *A, const double *R) {
  mL->Ap_ += A[0];
  mL->Ae_ += A[1];
  mR->Aw_ += A[2];
  mR->Ap_ += A[3];
  mL->rhs_ += R[0];
  mR->rhs_ += R[1];
}

//--------------------------------------------------------------------------
//-------- dirichlet -------------------------------------------------------
//--------------------------------------------------------------------------
void 
PointInfo::dirichlet(NaluOneDMatrixSystem *mS, double du) {
  mS->Ap_ = 1.0;
  mS->Aw_ = 0.0;
  mS->Ae_ = 0.0;
  mS->rhs_ = du;
}

//--------------------------------------------------------------------------
//-------- compute_visc ----------------------------------------------------
//--------------------------------------------------------------------------
void 
PointInfo::compute_visc(
  const double &tauWall, const double &rhoWall, const double &muWall, 
  const double kappa, std::vector<double> &coordsVec, std::vector<double> &viscVec) 
{
  // no user interface to Aplus
  const double Aplus = 17.0;
  const double yWall = coordsVec[0];
  for ( int k = 0; k < numPoints_; ++k ) {
    const double yk = coordsVec[k] - yWall;
    const double yplus = yk*std::sqrt(rhoWall*tauWall)/muWall;
    const double blend = 1.0 - std::exp(-yplus/Aplus);
    const double muTurb = kappa*rhoWall*std::sqrt(tauWall/rhoWall)*yk*blend*blend;
    viscVec[k] = muWall + muTurb;
  }
}

//--------------------------------------------------------------------------
//-------- assemble --------------------------------------------------------
//--------------------------------------------------------------------------
void 
PointInfo::assemble( 
  std::vector<NaluOneDMatrixSystem> &AmatrixSystemVec,
  std::vector<std::pair<int,int> > &lrscv, 
  std::vector<double> &coordsVec, 
  std::vector<double> &velocityVec, 
  std::vector<double> &viscVec)
{
  // size of traversal (points in the mesh, i.e., a row, and num edges)
  const size_t mSize = AmatrixSystemVec.size();
  const size_t numEdges = mSize-1; 

  // zero all of the matrix entries
  for ( size_t k = 0; k < mSize; ++k ) {
    NaluOneDMatrixSystem *thisEntry = &AmatrixSystemVec[k];
    thisEntry->zero();
  }
    
  // local data
  double A[4] = {};
  double R[2] = {};
  
  for ( size_t k = 0; k < numEdges; ++k ) {
    const int il = lrscv[k].first;
    const int ir = lrscv[k].second;
    
    NaluOneDMatrixSystem *mL = &AmatrixSystemVec[il];
    NaluOneDMatrixSystem *mR = &AmatrixSystemVec[ir];
    
    const double dx = (coordsVec[ir] - coordsVec[il]);
    const double lhsfac = -0.5*(viscVec[ir] + viscVec[il])/dx;
    // Ax = b (not currently in residual-form)
    const double diffFlux = lhsfac*(velocityVec[ir] - velocityVec[il])*0.0;
        
    // first left
    A[0] = -lhsfac;
    A[1] = +lhsfac;
    R[0] = -diffFlux;
    
    // now right
    A[2] = +lhsfac;
    A[3] = -lhsfac;
    R[1] = +diffFlux; 
    
    // scatter
    scatter(mL, mR, &A[0], &R[0]);
  }
}

//--------------------------------------------------------------------------
//-------- tdma ------------------------------------------------------------
//--------------------------------------------------------------------------
void 
PointInfo::tdma(
  int mSize, NaluOneDMatrixSystem *ns, double *u) 
{  
  // forward
  for ( int k = 1; k < mSize; ++k ) {
    const double m = ns[k].Aw_/ns[k-1].Ap_;
    ns[k].Ap_ -= m*ns[k-1].Ae_;
    ns[k].rhs_ -= m*ns[k-1].rhs_;
  }

  // backward
  u[mSize-1] = ns[mSize-1].rhs_/ns[mSize-1].Ap_;
  for ( int k = mSize-2; k >= 0; --k ) {
    u[k] = (ns[k].rhs_ - ns[k].Ae_*u[k+1])/ns[k].Ap_;
  }
}

//--------------------------------------------------------------------------
//-------- solve -----------------------------------------------------------
//--------------------------------------------------------------------------
double
PointInfo::solve(
  const double &uExchange, const double &uWall, 
  const double &rhoWall, const double &muWall, 
  const double &kappa, const double &tauWallProvided,
  bool &converged)
{
  const int MAX_ITER = 1000;
  const double TOLERANCE = 1.0e-6;

  // set some diagnostics
  iterations_ = MAX_ITER;
  
  // initial guess; use last if not the first time here
  double tauWallG = ( tauWall_ > 0.0 ) ? tauWall_ : tauWallProvided;

  for ( int ij = 0; ij < MAX_ITER; ++ij ) {  

    // update viscosity
    compute_visc(tauWallG, rhoWall, muWall, kappa, coordsVec_, viscVec_); 
    
    // assemble
    assemble(AmatrixSystemVec_, lrscv_, coordsVec_, velocityVec_, viscVec_);
    
    // correct for Dirichlet; Ax = b (not currently in residual-form)
    dirichlet(&AmatrixSystemVec_[0], uWall);
    dirichlet(&AmatrixSystemVec_[AmatrixSystemVec_.size()-1], uExchange);
    
    tdma(AmatrixSystemVec_.size(), &AmatrixSystemVec_[0], &velocityVec_[0]);
 
    // simple one-sided
    const double tauWallC = muWall
      *(velocityVec_[1] - velocityVec_[0])/(coordsVec_[1] - coordsVec_[0]);
    
    const double norm = (tauWallC-tauWallG)*(tauWallC-tauWallG);
    if ( std::sqrt(norm) < TOLERANCE ) {
      converged = true;
      iterations_ = ij+1;
      break;
    }

    // set the new tauwall
    tauWallG = tauWallC;
  }

  if ( debugOutput_ ) {
    for ( size_t k = 0; k < coordsVec_.size(); ++k ) {
      NaluEnv::self().naluOutput() << coordsVec_[k] << " " << velocityVec_[k] << std::endl;
    }
    throw std::runtime_error("END");
  }

  // return converged utau
  return std::sqrt(tauWallG/rhoWall);
}

} // namespace nalu
} // namespace sierra
