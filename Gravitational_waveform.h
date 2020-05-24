#include "Orbit.h"
#ifndef RESEARCH_PROJECT_GRAVITATIONAL_WAVEFORM_H
#define RESEARCH_PROJECT_GRAVITATIONAL_WAVEFORM_H

void Shape_Function(double eccentricity, double period, double eccentric_anomaly[],double true_anomaly[],double m1, double m2 ); //initializes

double H(double m1, double m2, double R, double eccentricity, double semi_major_a); //masses in Solar masses, R in parsecs, semi-major a in AU
double A0(double theta);
double B0(double theta);
double A1(double theta);
double B1(double theta);
double A2(double theta);
double B2(double theta);
double h_cross(double theta, double H, double eccentricity);
double h_plus(double theta, double H, double eccentricity);


void plot_waveform(double eccentricity, double semi_major_a, double mass_of_BH, double mass_of_OO, double R);
//eccentricity between 0 and 1
//semi major a, in AU
//period in seconds
//mass of the Black Hole, in Solar masses
//mass of the orbiting object, in Solar masses
//distance to the origin R, in parsecs

double AU_to_km(double AU);
double sec_to_yrs(double sec);
double yrs_to_sec(double yrs);
double SunM_to_kg(double SunM);
double kg_to_SunM(double kg);
double AU_to_meters(double AU);
double km_to_AU(double km);
double meters_to_AU(double meters);
double parsec_to_meter(double parsec);
double AU_to_parsec(double AU);
double parsec_to_AU(double parsec);
double kmsec_to_parsecsec(double kmsec);
double km2_to_parsec2(double km2);
double parsec2_to_AU2(double parsec2);
double km2_to_AU2(double km2);

//Kepler 3rd law
double period(double semi_major_a, double m1, double m2); //takes semi-major a in AU, returns period in seconds
double semi_major_a(double period, double m1, double m2); //takes period in seconds, returns semi-major a in AU
#endif //RESEARCH_PROJECT_GRAVITATIONAL_WAVEFORM_H
