/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ABLProfileFunction_h
#define ABLProfileFunction_h

namespace sierra{
namespace nalu{

class ABLProfileFunction
{
 public:

  ABLProfileFunction();
  virtual ~ABLProfileFunction();
  virtual double velocity(const double znorm) const = 0;
  virtual double temperature(const double znorm) const = 0;
};

class StableABLProfileFunction : public ABLProfileFunction
{
 public:
  StableABLProfileFunction(double gamma_m, double gamma_h);
  virtual ~StableABLProfileFunction();
  double velocity(const double znorm) const;
  double temperature(const double znorm) const;

 private:
  double gamma_m_;
  double gamma_h_;
};

class UnstableABLProfileFunction : public ABLProfileFunction
{
 public:
  UnstableABLProfileFunction(double beta_m, double beta_h);
  virtual ~UnstableABLProfileFunction();
  double velocity(const double znorm) const;
  double temperature(const double znorm) const;

 private:
  double beta_m_;
  double beta_h_;
  double pi_;
};

class NeutralABLProfileFunction : public ABLProfileFunction
{
 public:
  NeutralABLProfileFunction();
  virtual ~NeutralABLProfileFunction();
  double velocity(const double znorm) const;
  double temperature(const double znorm) const;
};

} // namespace nalu
} // namespace Sierra

#endif
