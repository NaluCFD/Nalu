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

template<typename T>
class PecletFunction
{
public:

  PecletFunction();
  virtual ~PecletFunction();
  virtual T execute(const T pecletNumber) = 0;
  /*virtual void update_values(Realm *realm) = 0;*/
};

template<typename T>
class ClassicPecletFunction : public PecletFunction<T>
{
public:
  ClassicPecletFunction(T A, T hf);
  virtual ~ClassicPecletFunction();
  T execute(const T pecletNumber);
  T A_;
  T hf_;
};

template<typename T>
class TanhFunction : public PecletFunction<T>
{
public:
  TanhFunction( T c1, T c2 );
  virtual ~TanhFunction();
  T execute(const T indVar);
  T c1_; // peclet number at which transition occurs
  T c2_; // width of the transtion
};

} // namespace nalu
} // namespace Sierra

#endif
