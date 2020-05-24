#ifndef RESEARCH_PROJECT_LIBRARY_H
#define RESEARCH_PROJECT_LIBRARY_H


//Orbit constructor, computes coordinates and writes them down to the file
void Orbit(double eccentricity, double apoapsis, double period);

void Orbit2(double eccentricity, double orbital_period, double m1, double m2);

double Bisection_Kepler(double eccentric_anomaly[], double mean_anomaly, double eccentricity);

struct Coordinate
{
    double x ;
    double y ;
};

double Kepler_Function(double eccentric_anomaly, double mean_anomaly, double eccentricity);

#endif