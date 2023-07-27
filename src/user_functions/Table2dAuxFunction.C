/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/Table2dAuxFunction.h>
#include <MeshMotionInfo.h>
#include <PecletFunction.h>
#include <SolutionOptions.h>

// basic c++
#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

Table2dAuxFunction::Table2dAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  std::vector<double> params,
  std::vector<std::string> stringParams) :
  AuxFunction(beginPos, endPos),
  fieldComp_(0),
  interpX_(0),
  interpY_(0),
  sizeOfX_(0),
  sizeOfY_(0),
  timeBlending_(0.0),
  baseFileName_("na.dat")
{
  // first extract the data
  if ( params.size() != 6 ) {
    throw std::runtime_error("Table2dAuxFunction requires params: interpComp, xComp, yComp, sizeX, sizeY, and time blending parameters");
  }
  else {
    fieldComp_ = params[0];
    interpX_ = params[1];
    interpY_ = params[2];
    sizeOfX_ = params[3];
    sizeOfY_ = params[4];
    timeBlending_ = params[5];
  }

  // second, extract file name
  if ( stringParams.size() != 1 ) {
    throw std::runtime_error("Table2dAuxFunction requires string params: file name");
  }
  else {
    baseFileName_ = stringParams[0];
  }

  /* 
     sample table data: ordered from min to max for both x and y

       x   y   h
     -100 -1   1
     -100  0   2
     -100 +1   3
      -50 -1   4
      -50  0   5
      -50 +1   6
       0  -1   7
       0  -0  20
       0  +1   8
      +50 -1   9
      +50  0   10
      +50 +1   11
     +100 -1   12
     +100  0   13
     +100 +1   14


     transformed data A_[m][n]: row-major access 

     y-axis (m-rows)
     |
     v
     ---------------------
     -1     1   4  7  9 12
      0     2   5 20 10 13 
     +1     3   6  8 11 14
     ---------------------
         -100 -50  0 50 100    <-  x-axis (n-columns)
  */

  // open file and load data
  std::ifstream inFile;
  inFile.open(baseFileName_.c_str());
  if ( !inFile.is_open() ) {
    throw std::runtime_error("Table2dAuxFunction file specified not available/opened-correctly");
  }

  // ignore the first line
  std::string line;
  std::getline(inFile, line);

  // create temp data
  std::vector<std::vector<double> > data;

  // read the data from the file
  double x, y, h;
  while ( !inFile.eof() ) {
    inFile >> x >> y >> h;
    std::vector<double> tmpData;
    tmpData.push_back(x);
    tmpData.push_back(y);
    tmpData.push_back(h);
    data.push_back(tmpData);
  }
  inFile.close();

  // fill in x-coordinates
  for ( size_t k = 0; k < data.size(); k += sizeOfY_ ) {
    std::vector<double> tmpData = data[k];
    const double xT = tmpData[0];
    tableX_.push_back(xT);
  }

  // fill in y coordinates
  for ( size_t k = 0; k < sizeOfY_; ++k ) {
    std::vector<double> tmpData = data[k];
    const double yT = tmpData[1];
    tableY_.push_back(yT);
  }

  // size A_
  A_.resize(sizeOfY_);
  for ( size_t i = 0; i < sizeOfY_; ++i )
    A_[i].resize(sizeOfX_);

  // load A_
  int countX = 0;
  for ( size_t k = 0; k < data.size(); ++k ) {
    std::vector<double> tmpData = data[k];
    const double h = tmpData[2];
    size_t ySlot = k % sizeOfY_;
    size_t xSlot = ySlot == sizeOfY_-1 ? countX++ : countX;
    A_[ySlot][xSlot] = h;
  }
}

Table2dAuxFunction::~Table2dAuxFunction()
{
  // nothing
}

void
Table2dAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  const double timeFac = std::min(time/timeBlending_, 1.0);
  for(unsigned p=0; p < numPoints; ++p) {

    double xp = coords[interpX_];
    double yp = coords[interpY_];

    int x1p = 0;
    int y1p = 0;
    int x2p = 0;
    int y2p = 0;
    
    find_entry(xp, x1p, x2p, tableX_);
    find_entry(yp, y1p, y2p, tableY_);
    
    const double x1 = tableX_[x1p];
    const double x2 = tableX_[x2p];
    const double y1 = tableY_[y1p];
    const double y2 = tableY_[y2p];
    
    const double Q11 = A_[y2p][x1p];
    const double Q12 = A_[y1p][x1p];
    const double Q21 = A_[y2p][x2p];
    const double Q22 = A_[y1p][x2p];
    
    const double R1 = Q11*((x2-xp)/(x2-x1)) + Q21*((xp-x1)/(x2-x1));
    const double R2 = Q12*((x2-xp)/(x2-x1)) + Q22*((xp-x1)/(x2-x1));
    const double interpValue = R1*((yp-y1)/(y2-y1)) + R2*((y2-yp)/(y2-y1));
      
    // first assign to zero
    for ( unsigned i = 0; i < fieldSize; ++i )
      fieldPtr[i] = 0.0;
    
    // assign
    fieldPtr[fieldComp_] = interpValue*timeFac;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

void
Table2dAuxFunction::find_entry(double &x, int &low, int &high, const std::vector<double> &table) const
{
  // check for low and high bounds
  if ( x <= table[0] ) {
    low = 0;
    high = 1;
    x = table[0];
  }
  else if ( x  >= table[table.size()-1] ) {
    high = table.size()-1;
    low = high - 1;
    x = table[table.size()-1];
  }
  else {
    // find the lower bound
    std::vector<double >::const_iterator it = std::lower_bound( table.begin(), table.end(), x );    
    if ( it == table.end() ) {
      std::cout<<"Value not found" << std::endl;
    }
    else {
      high = std::distance(table.begin(), it);
      low = high - 1;
    }
  }
}

} // namespace nalu
} // namespace Sierra
