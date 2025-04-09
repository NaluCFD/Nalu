#include <iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<sstream>
#include <stdexcept>

struct NaluOneDMatrixSystem {
  int globalId_;
  double Aw_, Ap_, Ae_;
  double rhs_;
  double dof_;
  NaluOneDMatrixSystem(unsigned int globalId)
    : globalId_(globalId), Aw_(0.0), Ap_(0.0), Ae_(0.0), rhs_(0.0), dof_(0.0) {}

  void zero() { Aw_ = 0.0; Ap_ = 0.0; Ae_ = 0.0; rhs_ = 0.0; dof_ = 0.0; }
};

void compute_courant_peclet(
  const double dt,
  std::vector<std::pair<int,int> > &lrscv, 
  std::vector<double> &coordsVec, 
  std::vector<double> &rhoVec, 
  std::vector<double> &diffFluxCoeffVec, 
  std::vector<double> &velocityVec)
{
  // size of traversal (points in the mesh, i.e., a row, and num edges)
  const size_t mSize = coordsVec.size();
  const size_t nEdges = mSize-1; 

  double courant = -1.0e6;
  double peclet = -1.0e6;
  for ( size_t k = 0; k < nEdges; ++k ) {
    const int il = lrscv[k].first;
    const int ir = lrscv[k].second;
        
    const double rhoIp = 0.5*(rhoVec[ir] + rhoVec[il]);
    const double muIp = 0.5*(diffFluxCoeffVec[ir] + diffFluxCoeffVec[il]);
    
    const double dx = (coordsVec[ir] - coordsVec[il]);
    const double uxIp = 0.5*(velocityVec[ir] + velocityVec[il]);
    courant = std::max(uxIp*dt/dx, courant);
    peclet = std::max(rhoIp*uxIp*dx/muIp, peclet);
  }
  std::cout << "Courant/Peclet number: " << courant << "/" << peclet << std::endl;
}

void scatter_edge(NaluOneDMatrixSystem *mL, NaluOneDMatrixSystem *mR, const double *A, const double *R) {
  mL->Ap_ += A[0];
  mL->Ae_ += A[1];
  mR->Aw_ += A[2];
  mR->Ap_ += A[3];
  mL->rhs_ += R[0];
  mR->rhs_ += R[1];
}

void scatter_node(NaluOneDMatrixSystem *mK, const double impl, const double src) {
  mK->Ap_ += impl;
  mK->rhs_ += src;
}

void dirichlet(NaluOneDMatrixSystem *mS, double du) {
  mS->Ap_ = 1.0;
  mS->Aw_ = 0.0;
  mS->Ae_ = 0.0;
  mS->rhs_ = du;
}

void zero_matrix(
  std::vector<NaluOneDMatrixSystem> &AmatrixSystemVec)
{
  // size of traversal (points in the mesh, i.e., a row, and num edges)
  const size_t mSize = AmatrixSystemVec.size();
  const size_t nEdges = mSize-1; 
  
  // zero all of the matrix entries
  for ( size_t k = 0; k < mSize; ++k ) {
    NaluOneDMatrixSystem *thisEntry = &AmatrixSystemVec[k];
    thisEntry->zero();
  }
}

void assemble_diffusion(
  std::vector<NaluOneDMatrixSystem> &AmatrixSystemVec,
  std::vector<std::pair<int,int> > &lrscv, 
  std::vector<double> &coordsVec, 
  std::vector<double> &scalarQnp1Vec, 
  std::vector<double> &diffFluxCoeffVec)
{
  // size of traversal (points in the mesh, i.e., a row, and num edges)
  const size_t mSize = AmatrixSystemVec.size();
  const size_t nEdges = mSize-1; 
    
  // local data
  double A[4] = {};
  double R[2] = {};
  
  for ( size_t k = 0; k < nEdges; ++k ) {
    const int il = lrscv[k].first;
    const int ir = lrscv[k].second;
    
    NaluOneDMatrixSystem *mL = &AmatrixSystemVec[il];
    NaluOneDMatrixSystem *mR = &AmatrixSystemVec[ir];
    
    const double dx = (coordsVec[ir] - coordsVec[il]);
    const double lhsfac = -0.5*(diffFluxCoeffVec[ir] + diffFluxCoeffVec[il])/dx;
    // Ax = b (not currently in residual-form and zeroed, M*dx = b - Ax)
    const double diffContribution = lhsfac*(scalarQnp1Vec[ir] - scalarQnp1Vec[il])*0.0;
    
    // first left
    A[0] = -lhsfac;
    A[1] = +lhsfac;
    R[0] = -diffContribution;
    
    // now right
    A[2] = +lhsfac;
    A[3] = -lhsfac;
    R[1] = +diffContribution; 
    
    // scatter
    scatter_edge(mL, mR, &A[0], &R[0]);
  }
}

void assemble_central_advection(
  std::vector<NaluOneDMatrixSystem> &AmatrixSystemVec,
  std::vector<std::pair<int,int> > &lrscv, 
  std::vector<double> &coordsVec, 
  std::vector<double> &scalarQnp1Vec, 
  std::vector<double> &densityVec,
  std::vector<double> &velocityVec)
{
  // size of traversal (points in the mesh, i.e., a row, and num edges)
  const size_t mSize = AmatrixSystemVec.size();
  const size_t nEdges = mSize-1; 
    
  // local data
  double A[4] = {};
  double R[2] = {};
  
  for ( size_t k = 0; k < nEdges; ++k ) {
    const int il = lrscv[k].first;
    const int ir = lrscv[k].second;
    
    NaluOneDMatrixSystem *mL = &AmatrixSystemVec[il];
    NaluOneDMatrixSystem *mR = &AmatrixSystemVec[ir];
    
    const double area = 1.0;
    const double mdot = 0.5*(densityVec[il]*velocityVec[il] + densityVec[ir]*velocityVec[ir])*area;
    const double lhsfac = 0.5*mdot;
    
    // Ax = b (not currently in residual-form, M*dx = b - Ax, and zeroed)
    const double advContribution = lhsfac*(scalarQnp1Vec[ir] + scalarQnp1Vec[il])*0.0;
    
    // first left
    A[0] = +lhsfac;
    A[1] = +lhsfac;
    R[0] = -advContribution;
    
    // now right
    A[2] = -lhsfac;
    A[3] = -lhsfac;
    R[1] = +advContribution; 
    
    // scatter
    scatter_edge(mL, mR, &A[0], &R[0]);
  }
}

void assemble_upwind_advection(
  std::vector<NaluOneDMatrixSystem> &AmatrixSystemVec,
  std::vector<std::pair<int,int> > &lrscv, 
  std::vector<double> &coordsVec, 
  std::vector<double> &scalarQnp1Vec, 
  std::vector<double> &densityVec,
  std::vector<double> &velocityVec)
{
  // size of traversal (points in the mesh, i.e., a row, and num edges)
  const size_t mSize = AmatrixSystemVec.size();
  const size_t nEdges = mSize-1; 
    
  // local data
  double A[4] = {};
  double R[2] = {};
  
  for ( size_t k = 0; k < nEdges; ++k ) {
    const int il = lrscv[k].first;
    const int ir = lrscv[k].second;
    
    NaluOneDMatrixSystem *mL = &AmatrixSystemVec[il];
    NaluOneDMatrixSystem *mR = &AmatrixSystemVec[ir];
    
    const double area = 1.0;
    const double mdot = 0.5*(densityVec[il]*velocityVec[il] + densityVec[ir]*velocityVec[ir])*area;
    const double lhsfacL = 0.5*(mdot+std::abs(mdot));
    const double lhsfacR = 0.5*(mdot-std::abs(mdot));

    // Ax = b (not currently in residual-form, M*dx = b - Ax, and zeroed)
    const double advContribution = (lhsfacL*scalarQnp1Vec[il] + lhsfacR*scalarQnp1Vec[il])*0.0;

    throw std::runtime_error("Implicit LHS contributions using upwind");

    // first left
    A[0] = 0.0; // fill in
    A[2] = 0.0; // fill in
    R[0] = -advContribution;
    
    // now right
    A[1] = 0.0; // fill in
    A[3] = 0.0; // fill in
    R[1] = +advContribution; 
    
    // scatter
    scatter_edge(mL, mR, &A[0], &R[0]);
  }
}

void assemble_time(
  std::vector<NaluOneDMatrixSystem> &AmatrixSystemVec,
  std::vector<double> &coordsVec,
  std::vector<double> &scalarQnp1Vec,
  std::vector<double> &scalarQnVec,
  std::vector<double> &densityVec,
  std::vector<double> &dualVolVec,
  const double dt)
{  
  // size of traversal (points in the mesh, i.e., a row, and num edges)
  const size_t mSize = AmatrixSystemVec.size();
  
  for ( size_t k = 0; k < mSize; ++k ) {
    
    NaluOneDMatrixSystem *mK = &AmatrixSystemVec[k];

    // for residual form
    double time = -(scalarQnp1Vec[k] - scalarQnVec[k])*dualVolVec[k]/dt*0.0;

    const double rhsfac = time; // fill in;
    const double lhsfac = 0.0;  // fill in;
    
    throw std::runtime_error("Implicit RHS and LHS contributions using a first-order Backward Euler scheme: rho*(phi^n+1 - phi^n)/dt*volume");

    // scatter
    scatter_node(mK, lhsfac, rhsfac);
  }
}

void output_state(
  const double time,
  std::vector<double> &coordsVec,
  std::vector<double> &scalarQnp1Vec) {

  std::ostringstream intToString;
  intToString << time;
  std::string fileNumAndType = intToString.str() + ".dat"; 
  std::string outputFileName = "fluids_" + fileNumAndType;

  std::ofstream outFile;
  outFile.open(outputFileName.c_str());

  // file open
  for ( size_t k = 0; k < scalarQnp1Vec.size(); ++k ) {
    outFile << time << " " << coordsVec[k] << " " << scalarQnp1Vec[k] << std::endl;
  }
  outFile.close();
}

void tdma(int mSize, NaluOneDMatrixSystem *ns, double *u) {
  
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

int main(){

  const bool debugOutput_ = false;
  const bool upwind = false;
  const double pi_ = std::acos(-1.0);
  const double omega_ = 1.0;
  
  // data layout
  /*
    0    1     2     3     4 ===> nPoints (coords, dofs, diffFluxCoeff)
    o----o-----o-----o-----o
      0     1     2     3 ======> nEdges (left/right node structure)
  */

  // user specified number of points: nPoints
  const int nPoints = 80;
  
  // derived number of segments: nEdges
  const int nEdges = nPoints - 1;
  
  // total domain size and time
  const double Ltotal = 1.0;  
  const double tFinal = 1.0;
  std::cout << "Ltotal and tTotal: " << Ltotal << " " << tFinal << std::endl;
  
  // user specified origin
  const double L0 = 0.0;
    
  // create vector of points and deltas
  std::vector<double> coordsVec_;
  std::vector<double> dxVec_;
  
  const double uniformS = Ltotal/(double)nEdges;
  dxVec_.push_back(uniformS);
  for ( int k = 1; k < nEdges; ++k ) {
    const double dyNext = uniformS;
    dxVec_.push_back(dyNext);
  }
  
  // populate coords
  coordsVec_.push_back(L0);
  double ccoords = 0.0;
  for ( size_t k = 0; k < dxVec_.size(); ++k ) {
    ccoords += dxVec_[k];
    coordsVec_.push_back(ccoords);
  }
  
  std::cout << "Total length is: " << ccoords << " " << Ltotal << std::endl;
  // populate scalarQ - a Gaussian
  std::vector<double> scalarQnVec_;
  std::vector<double> scalarQnp1Vec_;
  const double height_ = 1.0;
  const double center_ = 0.5;
  const double sigma_ = 0.05;

  for ( size_t k = 0; k < coordsVec_.size(); ++k ) {
    const double bracket = (coordsVec_[k] - center_);
    const double theV = height_*std::exp(-0.5*bracket*bracket/(sigma_*sigma_));
    scalarQnVec_.push_back(theV);
    scalarQnp1Vec_.push_back(theV);
  }
  
  // populate velocity, diffFluxCoeff and density vector
  std::vector<double> velocityVec_;
  std::vector<double> diffFluxCoeffVec_;
  std::vector<double> densityVec_;
  for ( size_t k = 0; k < coordsVec_.size(); ++k ) {
    velocityVec_.push_back(1.0);
    diffFluxCoeffVec_.push_back(0.01);
    densityVec_.push_back(1.0);
  }
  
  // output
  if ( debugOutput_ ) {
    std::cout << "Point/Coords" << std::endl;
    for ( size_t k = 0; k < coordsVec_.size(); ++k ) {
      std::cout << k << " " << coordsVec_[k] << std::endl;
    }
    
    // sanity check
    std::cout << "sanity check on bias: "
              << " " << dxVec_[dxVec_.size() -1] / dxVec_[0] << std::endl;
  }
  
  // start left/right mapping
  std::vector<std::pair<int,int> >lrscv_;
  for ( size_t k = 0; k < nEdges; ++k ) {
    lrscv_.push_back(std::make_pair(k,k+1));
  }
  
  if ( debugOutput_ ) {
    for ( size_t k = 0; k < nEdges; ++k ) {
      const int il = lrscv_[k].first;
      const int ir = lrscv_[k].second;
      std::cout << "Edge " << k << " left/right: " << il << "/" << ir << std::endl;
      std::cout << "  with coordinates: left/right: " << coordsVec_[il] << "/" << coordsVec_[ir] << std::endl;
    }
  }
  
  std::vector<double> dualVolVec_(nPoints, 0.0);
  if ( dualVolVec_.size() != coordsVec_.size() )
    std::cout << "there is a problem with dualVolVec_ size " << std::endl;
  
  for ( size_t k = 0; k < nEdges; ++k ) {
    const int il = lrscv_[k].first;
    const int ir = lrscv_[k].second;
    const double dx = coordsVec_[ir] - coordsVec_[il];
    // assemble dual volume
    dualVolVec_[il] += 0.5*dx;
    dualVolVec_[ir] += 0.5*dx;
  }

  // create a matrix object for each node
  std::vector<NaluOneDMatrixSystem> AmatrixSystemVec_;
  for ( int k = 0; k < nPoints; ++k ) {
    AmatrixSystemVec_.push_back(NaluOneDMatrixSystem(k));
  }
  
  // start the assemble and solve process (no need for iteration)
  const double dt_ = 0.01;
  double time_ = 0.0;
  output_state(time_, coordsVec_, scalarQnp1Vec_);

  while ( time_ <= tFinal ) {
    time_ += dt_;
    
    // set values at the left (for bc)
    const double QwallLeft = 0.0*sin(omega_*pi_*time_);
    scalarQnp1Vec_[0] = QwallLeft;

    compute_courant_peclet(dt_, lrscv_, coordsVec_, densityVec_, diffFluxCoeffVec_, velocityVec_);
    
    // zero matrix
    zero_matrix(AmatrixSystemVec_);
    
    // diffusion
    assemble_diffusion(AmatrixSystemVec_, lrscv_, coordsVec_, scalarQnp1Vec_, diffFluxCoeffVec_);

    // advection
    if ( upwind )
      assemble_upwind_advection(AmatrixSystemVec_, lrscv_, coordsVec_, scalarQnp1Vec_, densityVec_, velocityVec_);
    else
      assemble_central_advection(AmatrixSystemVec_, lrscv_, coordsVec_, scalarQnp1Vec_, densityVec_, velocityVec_);

    // time
    assemble_time(AmatrixSystemVec_, coordsVec_, scalarQnp1Vec_, scalarQnVec_, densityVec_, dualVolVec_, dt_);

    // Ax = b (not currently in residual-form, M*dx = b - Ax)
    dirichlet(&AmatrixSystemVec_[0], QwallLeft);

    // open boundary condirtion: advection out (unity area)
    const size_t mSizeM1 = AmatrixSystemVec_.size() - 1;
    NaluOneDMatrixSystem *mK = &AmatrixSystemVec_[mSizeM1];
    const double lhsAdvection = densityVec_[mSizeM1]*velocityVec_[mSizeM1];
    const double advOpenContribution = lhsAdvection*scalarQnp1Vec_[mSizeM1]*0.0;
    // Ax = b (not currently in residual-form, M*dx = b - Ax)
    mK->Ap_ += lhsAdvection;
    mK->rhs_ -= advOpenContribution;

    // solve
    tdma(AmatrixSystemVec_.size(), &AmatrixSystemVec_[0], &scalarQnp1Vec_[0]);

    // provide output
    output_state(time_, coordsVec_, scalarQnp1Vec_);

    // copy state
    for ( int k = 0; k < nPoints; ++k ) {
      scalarQnVec_[k] = scalarQnp1Vec_[k];
    }
    
  }
  
}
