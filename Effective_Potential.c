#include "Effective_Potential.h"
#include "Gravitational_waveform.h"
#include "Orbit.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>

#define PI 3.14159265
#define G 4.30091252525*pow(10,-3) //grav constant, in (parsec*km^2)/(Ms*sec^2)
#define c  300000 //speed of light, km/sec

double max(double [],int);

double max(double x[],int k)
{
    double t;
    t=x[0];
    for(int i=1;i<k;i++)
    {
        if(x[i]>t)
            t=x[i];
    }
    return(t);
}

//where m1 is the COM (BH) , and m2 is the orbiting body, in solar masses
//semi-major axes is in AU
//masses are in Solar masses
void plot_effective_potential(int* time_size, double eccentricity, double semi_major_a, double m1, double m2)
{
    double period_ = period(semi_major_a,m1,m2); //takes semi-major a in AU, masses in solar masses returns period in seconds
    double eccentric_anomaly[*time_size];
    double true_anomaly[*time_size];
    double r_vector[*time_size]; // in AU
    double effective_potential_N[*time_size];
    double effective_potential_S[*time_size];


    double periapsis = semi_major_a*(1-eccentricity); // in AU
    double Vp = tangential_speed(m1, m2, periapsis, semi_major_a); // km/sec
    //double angular_momentum_p = periapsis*km_to_AU(Vp)*m2;
    double angular_momentum_p = L_r(eccentricity, periapsis,m1);
    Shape_Function(eccentricity, period_,eccentric_anomaly,true_anomaly, SunM_to_kg(m1), SunM_to_kg(m2) );

    FILE* file1;
    FILE* file2;
    FILE* file3;
    FILE* energy;
    FILE* ang_mom;
    FILE* tang_speed;
    file1 = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\Newtonian.csv","w");
    file2 = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\Relativistic.csv","w");
    file3 = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\full.csv","w");
    tang_speed = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\tangential_speed.csv","w");
    energy = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\energy.csv","w");
    ang_mom = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\angular_momentum.csv","w");

    //fprintf(energy, "%e", orbital_energy(m2,m1,semi_major_a,angular_momentum_p,periapsis));
    double E = E_r(m1,periapsis,eccentricity);
    fprintf(ang_mom, "%e", km2_to_AU2(angular_momentum_p)); //output angular momentum is in km^2 / sec -> convert for further computations
    //fprintf(energy, "%e", max(effective_potential_S, *time_size));
    fprintf(energy, "%e", E);
    fprintf(tang_speed, "%f", km_to_AU(Vp)); //original Vp is in km per sec
    fclose(tang_speed);

    for (int i = 0; i <= *time_size; i++)
    {
        r_vector[i] =  r( semi_major_a, eccentricity,  true_anomaly[i]); // in AU
        effective_potential_N[i] = effective_potential(angular_momentum_p, m2, m1,r_vector[i]);
        effective_potential_S[i] = effective_potential_2( angular_momentum_p, m2, m1, r_vector[i]);
        fprintf(file1, "%f %e \n",r_vector[i] , effective_potential_N[i]);
        fprintf(file2, "%f %e \n",r_vector[i] , effective_potential_S[i]);

    }
    fclose(file1);
    fclose(file2);
    fclose(file3);
    fclose(energy);



}

double reduced_mass(double m1, double m2)
{
    return (m1*m2)/(m1+m2);
}

// both masses are in solar masses
//r in parsecs
double gravitational_potential_energy(double m1, double m2, double r)
{
    return (G*m1*m2)/r;
}

double r(double semi_a, double eccentricity, double theta)
{
    return (semi_a*(1-pow(eccentricity,2))) / (1+ eccentricity*cos(theta));
}

// Output is in km/sec
double tangential_speed(double m1, double m2,   double r, double semi_major_axs)
{
    double p1 = G*(m1+m2);
    double p2 = (2/AU_to_parsec(r))- (1/AU_to_parsec(semi_major_axs));
    return sqrt(p1*p2);
}

//r in AU
//m in solar masses
//tangential speed in AU/sec
double angular_momentum( double r, double m, double theta, double tangential_speed)
{
    return r*sin(theta)*m*tangential_speed;
}

double orbital_energy(double m, double M, double semi_major_a, double L, double r)
{
   double U_eff =  effective_potential_2(L,m,M,r);
   double v = tangential_speed(m,M,r,semi_major_a);
   return 0.5*m*km2_to_AU2(pow(v,2)) + U_eff ;
}

double effective_potential(double L, double m , double M, double r)
{
    //Newtonian
    double reduced = reduced_mass(m,M);
    return ((pow(L,2) / (2*reduced*pow(r,2))) - ((km2_to_AU2(G)*M*m) /(AU_to_parsec(r)))) / M;
    //return (pow(L,2) / (2*m*pow(r,2))) - ((km2_to_AU2(G)*M*m)/(AU_to_parsec(r)));

}

//L in AU^2*Mss / sec
//G in (parsec*km^2)/(Ms*sec^2)
// r in AU
// masses in Ms
double effective_potential_2(double L, double m, double M, double r)
{
 //Swarzchild AU^2 / sec^2*M
    double reduced = reduced_mass(m,M);
    double relativistic_correction = (parsec_to_AU(G)*(M+m)*pow(L,2))/(pow(c,2)*pow(r,3)*reduced);
    double V =  (pow(L,2) / (2*reduced*pow(r,2))) - ((km2_to_AU2(G)*M*m) /(AU_to_parsec(r))) - relativistic_correction;
    return V;


 //Swarzchild  (dimensionless)
 // https://hepweb.ucsd.edu/ph110b/110b_notes/node80.html
  // double rs = Schwarzchild_radius(M); // in AU
   //double dimensionless = (-rs/r) + (pow(L,2) / (pow(m,2)*pow(km_to_AU(c),2)))*((1/pow(r,2)) - (rs / pow(r,3)));
  //return dimensionless;

}

// G is in (parsec*km^2)/(Ms*sec^2)
//Mass is in Ms
//speed of light in km/sec
//the function converts Schwarzchild radius from parsec to AU
double Schwarzchild_radius(double mass_of_black_hole)
{
    return parsec_to_AU(2*G*mass_of_black_hole / pow(c,2));
}

double E_r(double M, double periaps_r, double e_r)
{
    double GM = parsec_to_AU(G) * M;
    double p1 = (1+e_r)*pow(c,2)*periaps_r;
    double num = GM*(1-e_r)*(p1 - 4*GM);
    double denum = pow(c,2)*periaps_r* (p1 - GM*(3 + pow(e_r,2)));
    double s = sqrt(1- (num/denum));
    return s*pow(c,2);
}

//input periapsis is in AU
// output angular momentum is in km^2 / sec
double L_r(double e_r, double periaps_r, double M)
{
    double p1 = 1 + e_r;
    double p2 = G*M;
    return km_to_AU((p1*c*sqrt(p2)*AU_to_km(periaps_r)) / sqrt (p1 * pow(c,2) * AU_to_parsec(periaps_r) - p2*(3 + pow(e_r,2))));
}






