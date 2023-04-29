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
#define num_bodies 3

using namespace std;

//---------------------------------------------------------------------------

class RungeKutta {
public:
    static void calculate(double h, vector<double>& u, vector<double> (*derivative)()) {
        double a[4] = {h/2, h/2, h, 0};
        double b[4] = {h/6, h/3, h/3, h/6};
        vector<double> u0(u.size()), ut(u.size());
        int dimension = u.size();

        for (int i = 0; i < dimension; i++) {
            u0[i] = u[i];
            ut[i] = 0;
        }

        for (int j = 0; j < 4; j++) {
            vector<double> du = derivative();

            for (int i = 0; i < dimension; i++) {
                u[i] = u0[i] + a[j]*du[i];
                ut[i] = ut[i] + b[j]*du[i];
            }
        }

        for (int i = 0; i < dimension; i++) {
            u[i] = u0[i] + ut[i];
        }
    }
};

//---------------------------------------------------------------------------

void init(double *r, double *v, double *m, double *u,int num_bodies){
    // Define initial positions and velocities

    //Earth
    r[0].x = 147.09e9;//r1_x
    r[0].y = 0.0;    //r1_y      
    v[0].x = 0.0;    //v1_x     
    v[0].y = 30290;  //v1_y   
    m[0] = mass_earth;
    
    //Moon
    r[1].x = 147.1e9; //r2_x      
    r[1].y = 0.3844e9;//r2_y     
    v[1].x = 0.0;     //v2_x      
    v[1].y = 1022.0;  //v2_y      
    m[1] = mass_moon;

    //Sun
    r[2].x = 0.0; //r3_x            
    r[2].y = 0.0; //r3_y     
    v[2].x = 0.0; //v3_x
    v[2].y = 0.0; //v3_y           
    m[2] = mass_sun;

    //instead of using two arrays for position and velocity, define one array for both.
    u = [r[0].x,r[0].y,v[0].x,v[0].y,r[1].x,r[1].y,v[1].x,v[1].y,r[2].x,r[2].y,v[2].x,v[2].y]



    // Print the initial conditions
    printf("Initial positions and velocities:\n");
    printf("Earth: r:(%f,%f) v:(%f,%f)\n", r[0].x,r[0].y,v[0].x,v[0].y);
    printf("Moon: r:(%f,%f) v:(%f,%f)\n",r[1].x,r[1].y,v[1].x,v[1].y);
    printf("Sun: r:(%f,%f) v:(%f,%f)\n", r[2].x,r[2].y,v[2].x,v[2].y);
}


//---------------------------------------------------------------------------
//calculate the derivative 

double* derivative() {
    double du;
    du = new double[num_bodies * 4];

    // Loop through the bodies
    for (int iBody = 0; iBody < initialConditions.bodies; iBody++) {
        // Starting index for current body in the u array
        int bodyStart = iBody * 4;

        du[bodyStart + 0] = state.u[bodyStart + 0 + 2]; // Velocity x
        du[bodyStart + 1] = state.u[bodyStart + 0 + 3]; // Velocity y
        du[bodyStart + 2] = acceleration(iBody, 0); // Acceleration x
        du[bodyStart + 3] = acceleration(iBody, 1); // Acceleration y
    }

    return du;
}

//---------------------------------------------------------------------------
//calculate the acceleration


// Calculates the acceleration of the body 'iFromBody'
// due to gravity from other bodies,
// using Newton's law of gravitation.
//   iFromBody: the index of body. 0 is first body, 1 is second body.
//   coordinate: 0 for x coordinate, 1 for y coordinate

double acceleration(int from, int coord,double *u, double *m){
  double result = 0.0;
  int iFrom = from * 4; //time 4 because the index of the pos/velocity array holds 4 values per body
  double distX, distY, dist, overDist3;

  //loop through the bodies
  for(int iTo = 0; iTo < num_bodies; iTo++){ //iterates through all of the bodies
    if(iFrom == iTo) { continue; } // iterates the next body
    int iTo = i * 4;
    // Distance between the two bodies
    distX = u[iTo + 0] - u[iFrom + 0]; //separation of each object in X
    distY = u[iTo + 1] - u[iFrom + 1]; //separation of each object in Y
    dist = sqrt(distX*distX + distY*distY); //calculate the total distance between the two objects
    overDist3 = 1/(dist*dist*dist); // 1/distance^3
    result += G*m[iTo](u[iTo + coord] - u[iFrom - coord])*overDist3; //the net force
  }
  return result;
}



//---------------------------------------------------------------------------

void updatePosition(timestep){
  RungeKutta::calculate(timestep, u,derivative);
}

void calculateNewPosition() {
  //loop through the bodies
  for(int i = 0; i < num_bodies; i++){
    int start = i * 4;

  }
}


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

    point2D *r,*v,*f;
    
    r=new double[num_bodies];
    v=new double[num_bodies];
    m=new double[num_bodies];
    u=new double[4*num_bodies];


    init(r,v,m,u,num_bodies);

	for(n=0;n<Nt;n++) {
        t=n*ht;
        h = 2;
        //time step
        //TODO: 1. Calculate Forces
        acceleration();
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
