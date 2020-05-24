#ifndef RESEARCH_PROJECT_EFFECTIVE_POTENTIAL_H
#define RESEARCH_PROJECT_EFFECTIVE_POTENTIAL_H

double reduced_mass(double m1, double m2);

double gravitational_potential_energy(double m1, double m2, double r);

double tangential_speed(double m1,double m2, double r, double semi_major_axs);

void plot_effective_potential(int* time_size, double eccentricity, double semi_major_a, double m1, double m2);

double r(double semi_a, double eccentricity, double theta);

double angular_momentum( double r, double m, double theta, double tangential_speed);

double effective_potential(double angular_momentum, double m, double M, double r);

 double effective_potential_2(double L, double m, double M, double r);

 double orbital_energy(double m, double M, double semi_major_a, double L, double r);

 double Schwarzchild_radius(double mass_of_black_hole);


double E_r(double M, double periaps_r, double e_r);

double L_r(double e_r, double periaps_r, double M);

#endif //RESEARCH_PROJECT_EFFECTIVE_POTENTIAL_H
