#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "Orbit.h"
#include "Gravitational_waveform.h"
#include "Effective_Potential.h"
//#include "Geodesic_Orbit.h"
#define dT 30 // timestep in seconds
#define Total_Time 0.008 // in years

int main() {
    double rs = Schwarzchild_radius(10^6); // in AU
    double eccentricity = 0.7;
    double M = 10^6;
    double m = 10^5;
    double semi_major_axis = 0.03291; // in AU
    double periapsis = semi_major_axis*(1-eccentricity);
    int* time_size = malloc(sizeof(int));
    *time_size = ceil (yrs_to_sec(Total_Time) / dT); //period here has to be in sec! this size is used for calculating the whole period
    if (periapsis <= rs)
    {
        printf("The object appears to be inside the black hole \n");
        return 1;
    }
    plot_waveform(eccentricity,semi_major_axis,M,m, 8000);
    plot_effective_potential(time_size,eccentricity,semi_major_axis,M,m);
    free(time_size);
    return 0;
}
