#include <stdio.h>
#include <iostream>
#include <math.h>

#include "fileutils.h"
#include "stringutils.h"
#include "gnuplot.h"

#define g0 9.8067 //default gravity
#define G = 6.673e-11 // Gravitational Constant

#define eps 1e-8

#define mass_M87 6.5e9 // solar mass of Messier 87
#define mass_sun 1
#define mass_earth 3e-6
#define mass_moon 3.69e-8
#define mass_jupiter 9.5e-4
#define mass_spaceship 5e-28 // mass of a 1000 kg spaceship

//LJ parameters
#define epsilon 1.0
#define sigma 1.0


typedef struct {
	double x,y;
} point2D;

typedef struct {
    double mass;
    double rx,ry,rz;
    double vx,vy,vz;
} body;


//---------------------------------------------------------------------------

void init(point2D *r,point2D *v){
        // Define initial positions and velocities

    //Earth
    r[0].x = 147.09e9;//r1_x
    r[0].y = 0.0;    //r1_y      
    v[0].x = 0.0;    //v1_x     
    v[0].y = 30290;  //v1_y   
    
    //Moon
    r[1].x = 147.1e9; //r2_x      
    r[1].y = 0.3844e9;//r2_y     
    v[1].x = 0.0;     //v2_x      
    v[1].y = 1022.0;  //v2_y      

    //Sun
    r[2].x = 0.0; //r3_x            
    r[2].y = 0.0; //r3_y     
    v[2].x = 0.0; //v3_x
    v[2].y = 0.0; //v3_y           


    // Print the initial conditions
    printf("Initial positions and velocities:\n");
    printf("Earth: r:(%f,%f) v:(%f,%f)\n", r[0].x,r[0].y,v[0].x,v[0].y);
    printf("Moon: r:(%f,%f) v:(%f,%f)\n",r[1].x,r[1].y,v[1].x,v[1].y);
    printf("Sun: r:(%f,%f) v:(%f,%f)\n", r[2].x,r[2].y,v[2].x,v[2].y);
}


//---------------------------------------------------------------------------
double F1(double f1){
    return f1;
}
double F2(double f2){
    return f2;
}
double F3(double f3){
    return f3;
}


//---------------------------------------------------------------------------


//---------------------------------------------------------------------------


int main() {
    //initial variables
	int N=100;
    int Nt=5000;
    double m=1.0;
    double ht=0.001;
    fHandle fT;
    fHandle fD;


    gnuplot *gp=new gnuplot();

    int n,i,h,j,fn,nc;
    double t;
    double x,y;
    int num_bodies = 3;//dimesions

    point2D *r,*v,*f;
    
    //object 1 (earth)
    r=new point2D[num_bodies];
    v=new point2D[num_bodies];
    f=new point2D[num_bodies]; 
    
    init(r,v);

	for(n=0;n<Nt;n++) {
        t=n*ht;
        h = 2;
        //time step
        //TODO: 1. Calculate Forces
        //TODO: 2. Update Momentum
        //TODO: 3. Update Position
        //TODO: 4. Write out time, pos 1, pos 2, pos 3

	}

    //graph the Temp. vs Temp Graph
    // gp->plotfile("temperature.txt","u 1:2 w l t 'Temp'");
    // gp->show();
    // gp->addcommand("reset");
    //graph the Temp. vs Temp Graph
    // gp->plotfile("density.txt","u 1:2 w l t 'Density'");
    // gp->show();
    // gp->addcommand("reset");


    // delete gp;
    //create movie from png images created by gnuplot
    // system("ffmpeg -y -i OUT/MD_%06d.png OUT/MD.m4v");
	
	delete[] r;
	delete[] v;
    delete[] f;
	delete gp;
	return 0;
}
