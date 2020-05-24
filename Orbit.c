#include "Orbit.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include "Gravitational_waveform.h"
#include "Effective_potential.h"


#define PI 3.14159265
#define ERROR 0.001
#define dT 30 // timestep in seconds
#define Total_Time 0.008 // in years



void Orbit(double eccentricity, double apoapsis, double period)
{
    size_t time_size = ceil (yrs_to_sec(Total_Time) / dT); //time vector is big enough to make a full period --> possibly make this vector bigger in the future
    double time[time_size];
    size_t size_of_ea_ordered = ceil(2*PI / ERROR);
    double ea_ordered[size_of_ea_ordered];
    double mean_anomaly [time_size];
    double eccentric_anomaly [time_size];
    double true_anomaly [time_size];
    struct Coordinate c[time_size];
    double periapsis = (apoapsis * (1-eccentricity)) / (1+ eccentricity);


    for (int j = 0; j <= size_of_ea_ordered; j++)
    {
        ea_ordered[j]  = ERROR * j;
    }


    for (int i = 0; i <= time_size; i++)
    {
        time[i] = dT * i;
        mean_anomaly[i] = (2*PI / period) * time[i];
    }


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

    FILE* file;
    file = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\position_data.txt","w");
    //computes position of the object
    for (int i = 0; i <= time_size; i++)
    {
        c[i].x = apoapsis*(cos(eccentric_anomaly[i]) - eccentricity);
        c[i].y = apoapsis*sqrt( 1 - pow(eccentricity,2.0)) * sin(eccentric_anomaly[i]);
        fprintf(file, "%f %f \n", c[i].x, c[i].y );
    }
    fclose (file);

}

void Orbit2(double eccentricity, double a1, double m1, double m2)
{
    size_t time_size = ceil (yrs_to_sec(Total_Time) / dT); //time vector is big enough to make a full period --> possibly make this vector bigger in the future
    double time[time_size];
    size_t size_of_ea_ordered = ceil(2*PI / ERROR);
    double ea_ordered[size_of_ea_ordered];
    double mean_anomaly [time_size];
    double eccentric_anomaly_one [time_size];
    double true_anomaly_one [time_size];
    double shape1 [time_size];
    struct Coordinate c1[time_size];

    double Ra1 = a1 * (1-eccentricity);


    for (int j = 0; j <= size_of_ea_ordered; j++)
    {
        ea_ordered[j]  = ERROR * j;
    }


    for (int i = 0; i <= time_size; i++)
    {
        time[i] = dT * i;
        mean_anomaly[i] = (2*PI / pow(a1,1.5)) * time[i];
    }


    //Calls bisection, computes position of the object, writes it to the file
    for (int k = 0; k <= time_size; k++)
    {
        eccentric_anomaly_one[k] = Bisection_Kepler(ea_ordered, mean_anomaly[k], eccentricity);
    }

    //Finds true anomaly using values of found eccentric anomalies
    for (int k = 0; k <= time_size; k++)
    {
        double right_side = sqrt((1+eccentricity) / (1-eccentricity)) * tan(eccentric_anomaly_one[k]/2);
        true_anomaly_one[k] = 2*atan(right_side);
    }

    FILE* file1;
    //FILE* file2;
    file1 = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\pos1.csv","w");
    //file2 = fopen("C:\\Users\\ladan\\Documents\\Research_Project\\pos2.csv","w");
    //computes shape equation of object 1 and 2
    for (int k = 0; k <= time_size; k++)
    {
        shape1[k] = Ra1*(1-eccentricity) / (1 + eccentricity*cos(true_anomaly_one[k]));
}
    //computes positions for orbit 1 and orbit 2, writes them down to the file
    for (int p = 0; p <= time_size; p++)
    {
        c1[p].x = shape1[p] * cos(true_anomaly_one[p]);
        c1[p].y = shape1[p] * sin(true_anomaly_one[p]);
        fprintf(file1, "%f %f \n", c1[p].x, c1[p].y );
    }

   /* shape2 = m1*shape1[1]/m2;
    c2.x = shape2 * cos(true_anomaly_one[1]+PI);
    c2.y = shape2 * sin(true_anomaly_one[1]+PI);
    fprintf(file2, "%f %f \n", c2.x, c2.y ); */

    fclose (file1);
    //fclose (file2);
}


double Bisection_Kepler(double eccentric_anomaly[], double mean_anomaly, double eccentricity)
{
    int start = 0;
    int end = ceil(2*PI / ERROR);
    int mid = (start + end) / 2;
    double start_value = Kepler_Function(eccentric_anomaly[start], mean_anomaly, eccentricity);
    double end_value = Kepler_Function(eccentric_anomaly[end], mean_anomaly, eccentricity);
    double mid_value = Kepler_Function(eccentric_anomaly[mid], mean_anomaly, eccentricity);

    if (start_value * end_value  > 0)
    {
        return -1;
    }

    while ( fabs(start_value - end_value) > ERROR)
    {
        if (fabs(start_value)  <= ERROR )
        {
            return eccentric_anomaly[start];
        }
        else if (fabs(end_value)  <= ERROR)
        {
            return eccentric_anomaly[end];
        }
        else if (fabs(mid_value)  <= ERROR)
        {
            return eccentric_anomaly[mid];
        }

        //RECOMPUTE START AND END AND MID INDEX + START AND END AND MID VALUES
        if (start_value * mid_value < 0)
        {
            end = mid;
            mid = (start + end) / 2;
            end_value = Kepler_Function(eccentric_anomaly[end], mean_anomaly, eccentricity);
            mid_value = Kepler_Function(eccentric_anomaly[mid], mean_anomaly, eccentricity);
        }
        else if (end_value * mid_value < 0)
        {
            start = mid;
            mid = (start + end )/ 2;
            start_value = Kepler_Function(eccentric_anomaly[start], mean_anomaly, eccentricity);
            mid_value = Kepler_Function(eccentric_anomaly[mid], mean_anomaly, eccentricity);
        }
    }

    return eccentric_anomaly[start];

};


double Kepler_Function(double eccentric_anomaly, double mean_anomaly, double eccentricity)
{
    double val = -eccentric_anomaly + (eccentricity*sin(eccentric_anomaly)) + fmod(mean_anomaly,2*PI);
    return val;
};