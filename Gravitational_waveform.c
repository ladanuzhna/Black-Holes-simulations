#include "Gravitational_waveform.h"
#include "Orbit.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>

#define PI 3.14159265
#define ERROR 0.001
#define G 6.67430151515*pow(10,-11) //grav constant, in m^3/(kg*s^2)
#define c 299792.458 //speed of light
#define Theta_n 0 //longitude of the ascending node
#define Theta_p 0 // argument of periapsis
#define Phi 0 //polarization angle
#define i 0 //inclination of the plane
#define dT 30 // timestep in seconds
#define Total_Time 0.008 // in years

void plot_waveform(double eccentricity, double semi_major_a, double m1, double m2, double R) // semi-major in AU, masses in Solar masses, distance R in parsecs
{
    double period_ = period(semi_major_a,m1,m2);
   size_t time_size = ceil (yrs_to_sec(Total_Time) / dT); //period here has to be in sec! this size is used for calculating the whole period
    double eccentric_anomaly[time_size];
    double true_anomaly[time_size];
    double h_cross_[time_size];
    double h_plus_[time_size];
    double H_ = H(m1, m2, R,eccentricity, semi_major_a);

    Shape_Function(eccentricity, period_,eccentric_anomaly,true_anomaly, m1, m2 );

    //do we calculate it for true anomaly???
    FILE* file;
    file = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\gravitational_waveform.csv","w");
    for (int a = 0; a <= time_size; a++)
    {
        h_cross_[a]= h_cross(true_anomaly[a],H_, eccentricity);
        h_plus_[a] = h_plus(true_anomaly[a],H_, eccentricity);
        fprintf(file, "%e %e \n", h_cross_[a], h_plus_[a]);
    }
    fclose(file);

}

void Shape_Function (double eccentricity, double period,double eccentric_anomaly[],double true_anomaly[], double m1, double m2 )
{
    size_t time_size = ceil (yrs_to_sec(Total_Time) / dT);
    size_t size_of_ea_ordered = ceil(2*PI / ERROR);
    double ea_ordered[size_of_ea_ordered];
    double time[time_size];
    double mean_anomaly [time_size];
    double shape [time_size];
    struct Coordinate c1[time_size];

    double Ra1 = semi_major_a(period, m1, m2) * (1-eccentricity);

    for (int j = 0; j <= size_of_ea_ordered; j++)
    {
        ea_ordered[j]  = ERROR * j;
    }

    FILE* file;
    file = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\time_vector.csv","w"); // time written in seconds
    for (int k = 0; k <= time_size; k++)
    {
        time[k] = dT * k;
        fprintf(file,"%f \n", time[k]);
        mean_anomaly[k] = (2*PI / period) * time[k];
    }
    fclose(file);


    //Calls bisection, computes position of the object, writes it to the file
    for (int k = 0; k <= time_size; k++)
    {
        eccentric_anomaly[k] = Bisection_Kepler(ea_ordered, mean_anomaly[k], eccentricity);
    }

    //Finds true anomaly using values of found eccentric anomalies
    for (int k = 0; k <= time_size; k++)
    {
        double right_side = sqrt((1+eccentricity) / (1-eccentricity)) * tan(eccentric_anomaly[k]/2);
        true_anomaly[k] = 2*atan(right_side);
    }

    FILE* file1;
    file1 = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\pos1.csv","w");
    //computes shape equation of object 1 and 2
    for (int k = 0; k <= time_size; k++)
    {
        shape[k] = Ra1*(1-eccentricity) / (1 + eccentricity*cos(true_anomaly[k]));
    }
    //computes positions for orbit 1 and orbit 2, writes them down to the file
    for (int p = 0; p <= time_size; p++)
    {
        c1[p].x = shape[p] * cos(true_anomaly[p]);
        c1[p].y = shape[p] * sin(true_anomaly[p]);
        fprintf(file1, "%f,%f \n", c1[p].x, c1[p].y );
    }

    fclose (file1);
}

double H (double m1, double m2, double R, double eccentricity, double semi_major_a) //masses in Solar masses, R in parsecs, semi-major a in AU
{
    double HH = (4*pow(G,2)*SunM_to_kg(m1)*SunM_to_kg(m2))/(pow(c,4)*AU_to_meters(semi_major_a)*(1-pow(eccentricity,2))*parsec_to_meter(R));
    return HH;
}

double A0(double theta)
{
    return -0.5*(1+pow(cos(i),2))*cos(2*(theta-Theta_n));
}

double B0(double theta)
{
    return -cos(i)*sin(2*(theta-Theta_n));
}

double A1(double theta)
{
   double part_1 =  0.25*pow(sin(i),2)*cos(theta-Theta_p);
   double part_2 = ((1+ pow(cos(i),2))*(5*cos(theta-2*Theta_n+Theta_p)))/8;
   double part_3 = cos(3*theta-2*Theta_n-Theta_p);
   return part_1-part_2+part_3;
}

double B1(double theta)
{
    double brackets =  5*sin(theta-2*Theta_n+Theta_p)+sin(3*theta-2*Theta_n-Theta_p);
    return -0.25*cos(i)*brackets;
}

double A2(double theta)
{
    return 0.25*pow(sin(i),2)-0.25*(1+pow(cos(i),2))*cos(2*(Theta_n-Theta_p));
}

double B2(double theta)
{
    return 0.5*cos(i)*sin(2*(Theta_n-Theta_p));
}

double h_cross(double theta,double H, double eccentricity)
{
    return H*(sin(2*Phi)*(A0(theta)+eccentricity*A1(theta)+pow(eccentricity,2)*A2(theta))+cos(2*Phi)*(B0(theta)+eccentricity*B1(theta)+pow(eccentricity,2)*B2(theta)));
}

double h_plus(double theta, double H, double eccentricity)
{
    return H*(cos(2*Phi)*(A0(theta)+eccentricity*A1(theta)+pow(eccentricity,2)*A2(theta))-sin(2*Phi)*(B0(theta)+eccentricity*B1(theta)+pow(eccentricity,2)*B2(theta)));
}

double sec_to_yrs(double sec)
{
    return sec*3.17*pow(10,-8);
}

double yrs_to_sec(double yrs)
{
    return yrs*3.154*pow(10,7);
}

double period(double semi_major_a, double m1, double m2) //takes semi-major a in AU, masses in solar masses returns period in seconds
{
    double period =   pow(pow(AU_to_meters(semi_major_a),3)*4*pow(PI,2)/(G*(SunM_to_kg(m1+m2))),0.5);
    return period;
}

double semi_major_a(double period, double m1, double m2) //takes period in seconds, masses in solar masses, returns semi-major a in AU
{
     return meters_to_AU(pow(pow(period,2)*G*(SunM_to_kg(m1+m2))/(4*pow(PI,2)),0.3333333));
}

double SunM_to_kg(double SunM)
{
    return SunM*1.98847*pow(10,30);
}

double kg_to_SunM(double kg)
{
    return kg*5.02785431*pow(10,-31);
}

double AU_to_meters(double AU)
{
    return AU*149597870691;
}

double AU_to_km(double AU)
{
    return AU*1.496*pow(10,8);
}

double meters_to_AU(double meters)
{
    return meters*6.68458712267059*pow(10,-12);
}

double parsec_to_meter(double parsec)
{
    return parsec*3.086*pow(10,16);
}

double AU_to_parsec(double AU)
{
 return AU*0.0000048481367817234;
}

double parsec_to_AU(double parsec)
{
    return parsec*206264.80624538;
}


double km2_to_parsec2(double km2)
{
    return km2*1.0503*pow(10,-27);
}

double parsec2_to_AU2(double parsec2)
{
    return parsec2 * 4.25451703 * pow(10, 10);
}

double km2_to_AU2(double km2)
{
 return km2*4.46837*pow(10,-17);
}

double km_to_AU(double km)
{
    return km*6.684587124056*pow(10,-9);
}