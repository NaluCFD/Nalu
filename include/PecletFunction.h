/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PecletFunction_h
#define PecletFunction_h

namespace sierra{
namespace nalu{

class PecletFunction
{
public:

  PecletFunction();
  virtual ~PecletFunction();
  virtual double execute(const double pecletNumber) = 0;
  /*virtual void update_values(Realm *realm) = 0;*/
};

class ClassicPecletFunction : public PecletFunction
{
public:
  ClassicPecletFunction(double A, double hf);
  virtual ~ClassicPecletFunction();
  double execute(const double pecletNumber);
  double A_;
  double hf_;
};

class TanhPecletFunction : public PecletFunction
{
public:
  TanhPecletFunction( double c1, double c2 );
  virtual ~TanhPecletFunction();
  double execute(const double pecletNumber);
  double c1_; // peclet number at which transition occurs
  double c2_; // width of the transtion
  double shift_;
  double delta_;
};

} // namespace nalu
} // namespace Sierra

#endif
