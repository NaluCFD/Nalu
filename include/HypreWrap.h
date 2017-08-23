/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef HYPREWRAP_H
#define HYPREWRAP_H

namespace sierra {
namespace nalu {

namespace hypre {

/** Defines the supported solvers and preconditioners from the Hypre package.
 *
 *  Note that some can be used as both solvers and/or preconditioners, whereas
 *  others are strictly either solver or a preconditioner. The choice is
 *  indicated by the sierra::nalu::hypre::HYPRE_Chooser type. See Hypre Manual
 *  for more details on the solvers and preconditioners.
 *
 *  This section of the code has been adapated from the xSDKTrilinos package.
 */
enum HYPRE_Solver{
  BoomerAMG,
  ParaSails,
  Euclid,
  AMS,
  Hybrid,
  PCG,
  GMRES,
  FlexGMRES,
  LGMRES,
  BiCGSTAB
};

/** Flag differentiating Hypre Solvers and Preconditioners for options etc.
 *
 *  See sierra::nalu::hypre::HYPRE_Solver for currently available
 *  solver/preconditioner options.
 */
enum Hypre_Chooser {
  Solver,
  Prec
};

/** Wrapper class to store HYPRE input options to be passed to solver and
 * preconditioner.
 *
 *  This class is modeled after the FunctionParameter class in
 */
class HypreOptions
{
public:
  //! Single int constructor.
  HypreOptions(
    Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, int), int param1)
    : chooser_(chooser), option_(0), int_func_(funct_name), int_param1_(param1)
  {}

  //! Single double constructor.
  HypreOptions(
    Hypre_Chooser chooser,
    int (*funct_name)(HYPRE_Solver, double),
    double param1)
    : chooser_(chooser),
      option_(1),
      double_func_(funct_name),
      double_param1_(param1)
  {}

  //! Single double, single int constructor.
  HypreOptions(
    Hypre_Chooser chooser,
    int (*funct_name)(HYPRE_Solver, double, int),
    double param1,
    int param2)
    : chooser_(chooser),
      option_(2),
      double_int_func_(funct_name),
      int_param1_(param2),
      double_param1_(param1)
  {}

  //! Two ints constructor.
  HypreOptions(
    Hypre_Chooser chooser,
    int (*funct_name)(HYPRE_Solver, int, int),
    int param1,
    int param2)
    : chooser_(chooser),
      option_(3),
      int_int_func_(funct_name),
      int_param1_(param1),
      int_param2_(param2)
  {}

  //! Int pointer constructor.
  HypreOptions(
    Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, int*), int* param1)
    : chooser_(chooser),
      option_(4),
      int_star_func_(funct_name),
      int_star_param_(param1)
  {}

  //! Double pointer constructor.
  HypreOptions(
    Hypre_Chooser chooser,
    int (*funct_name)(HYPRE_Solver, double*),
    double* param1)
    : chooser_(chooser),
      option_(5),
      double_star_func_(funct_name),
      double_star_param_(param1)
  {}

  //! char pointer constructor
  HypreOptions(
    Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, char*), char* param1)
    : chooser_(chooser),
      option_(6),
      char_star_func_(funct_name),
      char_star_param_(param1)
  {}

private:
  Hypre_Chooser chooser_;
  int option_;

  int (*int_func_)(HYPRE_Solver, int);
  int (*double_func_)(HYPRE_Solver, double);
  int (*double_int_func_)(HYPRE_Solver, double, int);
  int (*int_int_func_)(HYPRE_Solver, int, int);
  int (*int_star_func_)(HYPRE_Solver, int*);
  int (*double_star_func_)(HYPRE_Solver, double*);
  int (*char_star_func_)(HYPRE_Solver, char*);
  int int_param1_;
  int int_param2_;
  double double_param1_;
  int *int_star_param_;
  double *double_star_param_;
  char* char_star_param_;
};

}

}  // nalu
}  // sierra


#endif /* HYPREWRAP_H */
