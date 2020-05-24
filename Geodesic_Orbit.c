#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "Effective_Potential.h"

#define PI 3.14159265
#define G 4.30091252525*pow(10,-3) //grav constant, in (parsec*km^2)/(Ms*sec^2)
#define c  300000 //speed of light, km/sec


void Geodesic_Orbit(int* time_size, double M, double periapsis, double e)
{
    int dimension = 2;		/* number of differential equations */

    double eps_abs = 1.e-8;	/* absolute error requested */
    double eps_rel = 1.e-10;	/* relative error requested */

    /* define the type of routine for making steps: */
    const gsl_odeiv_step_type *type_ptr = gsl_odeiv_step_rkf45;

    /*
       allocate/initialize the stepper, the control function, and the
       evolution function.
    */
    gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dimension);
    gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (eps_abs, eps_rel);
    gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dimension);

    gsl_odeiv_system my_system;	/* structure with the rhs function, etc. */

    double E = E_r(M, periapsis, e);
    double L = L_r(e,periapsis,M);
    double r[*time_size];			/* current solution vector */

    double T, T_next;		/* current and next independent variable */
    double Tmin, Tmax, delta_t;	/* range of t and step size for output */

    double h = 1e-6;		/* starting step size for ode solver */
}
