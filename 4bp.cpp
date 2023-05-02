#include <stdio.h>
#include <iostream>
#include <math.h>

#include "fileutils.h"
#include "stringutils.h"
#include "gnuplot.h"
// #define G = 6.673e-11 // Gravitational Constant

#define eps 1e-8

#define mass_M87  4.7734080027e+39// solar mass of Messier 87
#define mass_sun 1.98855e30
#define mass_earth 5.972e24
#define mass_moon 7.347e22 
#define mass_jupiter 1.898e28
#define mass_spaceship 1e4 // mass of a 1000 kg spaceship

#define num_bodies 4

using namespace std;

// don't really need this one since we'll eventually have the state structure
typedef struct {
	double x,y; //position x,y  
} point2D;

typedef struct {
	double x,y; //position x,y 
    double mass; //mass 
} state;


//---------------------------------------------------------------------------

void init(point2D *r, point2D *v, double *m, double *u){
    // Define initial positions and velocities

    //Earth
    r[0].x = 147.09e8;//r1_x
    r[0].y = 0.0;    //r1_y      
    v[0].x = 0.0;    //v1_x     
    v[0].y = 30290;  //v1_y   
    m[0] = mass_earth;
    
    //Moon
    r[1].x = 147.1e8; //r2_x      
    r[1].y = 0.3844e8;//r2_y     
    v[1].x = 0.0;     //v2_x      
    v[1].y = 1022.0;  //v2_y      
    m[1] = mass_moon;

    //Sun
    r[2].x = 0.0; //r3_x            
    r[2].y = 0.0; //r3_y     
    v[2].x = 1000.0; //v3_x
    v[2].y = 0.0; //v3_y           
    m[2] = mass_sun;

    //Jupiter
    r[3].x = 5.0605e+10; //r3_x            
    r[3].y = 0.0; //r3_y     
    v[3].x = 0.0; //v3_x
    v[3].y = 20.0; //v3_y           
    m[3] = mass_M87;

   
    //instead of using two arrays for position and velocity, define one array for both.
    // u = {r[0].x,r[0].y,v[0].x,v[0].y,r[1].x,r[1].y,v[1].x,v[1].y,r[2].x,r[2].y,v[2].x,v[2].y};
    for(int i = 0; i < num_bodies; i++){
        u[4*i] =  r[i].x;
        u[(4*i)+1] = r[i].y;
        u[(4*i)+2] = v[i].x;
        u[(4*i)+3] = v[i].y;
    }

}

//---------------------------------------------------------------------------
//calculate the acceleration

//credit for this design goes to Evgenii 

//acceleration for each component, is called a total: 6 times for 3 bodies
double acceleration(int fromBody, int coord, double *u, double *m){
    double result = 0.0;
    int fromBodyStart = fromBody * 4; //time 4 because the index of the pos/velocity array holds 4 values per body
    double distX, distY, dist, overDist3;

    double G = 6.673e-11;
    // double G = 1.0;

    //loop through the bodies
    for(int iToBody = 0; iToBody < num_bodies; iToBody++){ //iterates through all of the bodies
        if(fromBody == iToBody) { continue; } // is the same, iterate to the next body
        int iToBodyStart = iToBody * 4;

        // Distance between the two bodies
        distX = u[iToBodyStart + 0] - u[fromBodyStart + 0]; //separation of each object in X
        distY = u[iToBodyStart + 1] - u[fromBodyStart + 1]; //separation of each object in Y
        
        dist = sqrt(distX*distX + distY*distY); //calculate the total distance between the two objects
        
        // printf("distance = %f\n",dist);

        //okay this works
        overDist3 = 1/(dist*dist*dist); // 1/distance^3
        // printf("overDist3 = %.50f\n",overDist3);

        result += G*m[iToBody]*(u[iToBodyStart + coord] - u[fromBodyStart + coord])*overDist3; //the net force
    }

    // printf("acceleration = %f\n",result);
    return result;
}



//---------------------------------------------------------------------------
//calculate the derivative 

double * derivative(double *u, double *m) {
    double *du;
    du = new double[num_bodies * 4];

    // Loop through the bodies
    for (int iBody = 0; iBody < num_bodies; iBody++) {
        // Starting index for current body in the u array
        int bodyStart = iBody * 4;
        du[bodyStart + 0] = u[bodyStart + 2]; // dr_(bodyStart)_x = v_(bodyStart)_x 
        du[bodyStart + 1] = u[bodyStart + 3]; // dr_(bodyStart)_y = v_(bodyStart)_y
        
        //PROBLEM !!!!
        du[bodyStart + 2] = acceleration(iBody, 0, u, m); // Acceleration x
        du[bodyStart + 3] = acceleration(iBody, 1, u, m); // Acceleration y
        // printf("du[bodyStart + 2] here is %f\n",du[bodyStart + 3]);
    
    }
    return du;
}

//---------------------------------------------------------------------------

//RK4 to update position and velocity all at once
static void calculate(double h, double *u, double *m) {
    double a[4] = {h/2, h/2, h, 0}; // for RK4
    double b[4] = {h/6, h/3, h/3, h/6}; // for RK4
    double *u0, *ut;
    int uSize = sizeof(u);

    u0 = new double[uSize];
    ut = new double[uSize];

    int dimensionOfArray = uSize;

    for (int i = 0; i < dimensionOfArray; i++) {
        u0[i] = u[i]; // keep our initial value the same
        // if(i==0) printf("u[i] here is %f\n",u[i]); // FINE HERE
        ut[i] = 0.0; // we want to reset and previously existing value here
    }

    for (int j = 0; j < 4; j++) { // over the number of steps for RK4

        //PROBLEM
        double *du = derivative(u,m); // need the derivatives for both position and velocity
    
    
        for (int i = 0; i < dimensionOfArray; i++) { // goes to 4*numBodies = (12 for 3 bodies, 16 for 4 bodies)
            // printf("du[i] here is %f\n",du[i]);
            u[i] = u0[i] + a[j]*du[i]; // initial value
            ut[i] = ut[i] + b[j]*du[i]; // time stepped value
        }
    }

    for (int i = 0; i < dimensionOfArray; i++) {
        u[i] = u0[i] + ut[i];
        
    }
}

//---------------------------------------------------------------------------


int main() {
    //initial variables
    int N = num_bodies;
	int Nt=10000000;
    double ht=0.5;


    fHandle f;
    f = FileCreate("4bp.txt");

    gnuplot *gp=new gnuplot();

    int n,fn;
    double t;
    double x,y;
    string s;

    point2D *r,*v;
    double *u,*m;
    
    r=new point2D[num_bodies];
    v=new point2D[num_bodies];
    u=new double[4*N];
    m=new double[num_bodies];

    init(r,v,m,u); //initialize it first here

	for(n=0;n<Nt;n++) {
        t=n*ht;
        init(r,v,m,u); //reinitilize it after every single time step
        calculate(t, u, m);
        //TODO: 4. Write out time, pos 1, pos 2, pos 3
        s = FloatToStr(t)+"\t"+ FloatToStr(u[0]) + "\t"+FloatToStr(u[1]) + "\t"
                    + FloatToStr(u[4]) + "\t"+FloatToStr(u[5]) + "\t"
                    + FloatToStr(u[8]) + "\t"+FloatToStr(u[9]) + "\t"
                    + FloatToStr(u[12]) + "\t"+FloatToStr(u[13]) + "\n";
        FileWrite(f,s.c_str(),s.length());
	}
    FileClose(f);

    //graph the posX vs posY Graph
    gp->plotfile("4bp.txt","u 2:3 w l t 'Earth'");
    gp->replotfile("4bp.txt","u 4:5 w l t 'Moon'");
    gp->replotfile("4bp.txt","u 6:7 w l t 'Sun'");
    gp->replotfile("4bp.txt","u 8:9 w l t 'M87'");
    gp->addcommand("reset");
    gp->show();


    delete gp;
    //create movie from png images created by gnuplot
    // system("ffmpeg -y -i OUT/3bp_%06d.png OUT/3bp.m4v");
	
	delete[] r;
	delete[] v;
	delete[] u;
    delete[] m;
	delete gp;
	return 0;
}