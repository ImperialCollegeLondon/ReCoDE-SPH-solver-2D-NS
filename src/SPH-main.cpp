//This is a main programme to run the SPH simulation

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include "class.h"
#include <mpi.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;


/**defining functions for the validation cases 
 *(explantations are provided at their implementation bellow the main programme)
 **/

void ic_one_particle(int n, SPH & sph);

void ic_two_particles(int n, SPH & sph);

void ic_three_particles(int n, SPH & sph);

void ic_four_particles(int n, SPH & sph);

void ic_dam_break(int n, SPH & sph);

void ic_block_drop(int n, int n1, int n2, SPH & sph);

void ic_droplet(int n, SPH & sph);

int dropletn(int n);

//Start of the main programme
int main(int argc, char* argv[]){
    
    //Process to obtain the directions provided by the user through the command line
    po::options_description desc("Allowed options");
     desc.add_options()
        ("ic-one-particle","take an initial condition for the case of one particle")
        ("ic-two-particles","take an initial condition for the case of two particles")
        ("ic-three-particles","take an initial condition for the case of three particles")
        ("ic-four-particles","take an initial condition for the case of four particles")
        ("ic-dam-break","take an initial condition for the case of dam break")
        ("ic-block-drop","take an initial condition for the case of block drop")
        ("ic-droplet","take an initial condition for the case of droplet")
        ("T", po::value<double>(), "take integration time")
        ("dt", po::value<double>(), "take time-step")
        ("h", po::value<double>(), "take radius of influence")
    ;

    po::variables_map vm;    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    int n,n1,n2;                            //n is defined as the number of particles, while n1 and n2 are only used for the droplet case and block drop cases

    double T1 = vm["T"].as<double>();       //Total integration time
    double dt = vm["dt"].as<double>();      //Time step dt
    double h = vm["h"].as<double>();        //Radius of influence

    int T = int(T1/dt) + 1;                 //Transforming the time in seconds to iterations
 
    if (vm.count("ic-one-particle")) { n = 1; }
    
    if (vm.count("ic-two-particles")) { n = 2; }

    if (vm.count("ic-three-particles")) { n = 3; }

    if (vm.count("ic-four-particles")) { n = 4; }

    if (vm.count("ic-dam-break")) { n = 400; }

    if (vm.count("ic-block-drop")) { 

        //n1 and n2 are such that dx and dy are equal, so that there can be a unifrom distribution
        n1=17;
        n2=25;

        n=n1*n2; // Total number of particles
        
    
     }
    
    if (vm.count("ic-droplet")) { 
        
        n = 400;
    
        /**Firstly the exact nuber of particles needs to be determined.
         * This is because the initial condition is formed by firstly 
         * creating a square grid and then tranformig it into a circle grid.
         * Therefore, the number of particles that remain in the final grid is
         * determined by the function dropletn
         **/

        //The initial number of particles is saved to be used later
        n1 = n;

        //Replacing n with the actual particles that will be passed inside the class
        n = dropletn(n);
        
    }

    //Defining the solver object called sph for the first time.
    //In its definition, the number of particles is required
    SPH sph(n);

    /**After the number of partciles was introduced inside the class and
     * therefore the appropriate matrices were initialized, the particles 
     * are ordered in the correct positions
     **/

    if (vm.count("ic-one-particle")) {

        ic_one_particle(n,sph);
        sph >> dt;
        sph < h;

    }

    if (vm.count("ic-two-particles")) {

        ic_two_particles(n,sph);
        sph >> dt;
        sph < h;
    }

    if (vm.count("ic-three-particles")) {


        ic_three_particles(n,sph);
        sph >> dt;
        sph < h;
    }

    if (vm.count("ic-four-particles")) {

        ic_four_particles(n,sph);
        sph >> dt;
        sph < h;
    }

    if (vm.count("ic-dam-break")) {

        ic_dam_break(n,sph);
        sph >> dt;
        sph < h;
    }

    if (vm.count("ic-block-drop")) {

        ic_block_drop(n,n1,n2,sph);
        sph >> dt;
        sph < h;
    }

    if (vm.count("ic-droplet")) {
    
        ic_droplet(n1,sph);
        sph >> dt;
        sph < h;
        
    }
    
    /**This is a step that helps to split the original matrix through which the values of the 
     * initial coordinates and the initial velocities were introduced inside the class.
     **/
    sph.x0();
    sph.y0();
    sph.vx0();
    sph.vy0();

    //The function to calculate the mass of the partciles is called
    sph.mass();

    ofstream vOut("Positions-x-y.txt", ios::out | ios::trunc);
    ofstream vOut2("Energy-File.txt", ios::out | ios::trunc);
    

    vOut.precision(5);
    vOut << "x" << "          " << "y"  << "\n" ;
    vOut2.precision(5);
    vOut2 << "t" << "      " << "Ek" << "       " << "Ep" << "     " << "Etotal"  << "\n" ;



    //Time integration loop
    for (int t = 0; t < T; t++){

        //passing the specific "time" of the loop inside the class
        sph > t;
        
        //In each iteration the disatnces between the particles are recalculated, as well as their densities
        sph.rVec(); 
        sph.den();   
        sph.spatial();
        //sph.getdata();
     
  

        //Writing the energies on the Energy-File
        vOut2 << t*dt << "  " << sph.Ek() << "  " << sph.Ep() << "  " << sph.Ep() + sph.Ek()  << "\n" ; 

            //Getting the posistions after the integration is completed
            if(t==T-1){


                for (int l = 0; l < n; l++){

                    vOut << sph.retx(l) << " " << sph.rety(l)  << "\n" ; 

                }
            }
              
    }



    
    return 0;
    
}

//Initial condition with one particle
void ic_one_particle(int n, SPH & sph){


    sph(0,0) = 0.5;
    sph(1,0) = 0.5;
    sph(2,0) = 0.0;
    sph(3,0) = 0.0;

}

//Initial condition with two particles
void ic_two_particles(int n, SPH & sph){


    sph(0,0) = 0.5;
    sph(0,1) = 0.5;

    sph(1,0) = 0.5;
    sph(1,1)= 0.01;

    sph(2,0) = 0.0;
    sph(2,1) = 0.0;    

    sph(3,0) = 0.0;
    sph(3,1) = 0.0;

}

//Initial condition with three particles
void ic_three_particles(int n, SPH & sph){

    sph(0,0) = 0.5;
    sph(0,1) = 0.495;
    sph(0,2) = 0.505;

    sph(1,0) = 0.5;
    sph(1,1)= 0.01;
    sph(1,2)= 0.01;

    sph(2,0) = 0.0;
    sph(2,1) = 0.0; 
    sph(2,2) = 0.0;   

    sph(3,0) = 0.0;
    sph(3,1) = 0.0;
    sph(3,2) = 0.0;

}

//Initial condition with four particles
void ic_four_particles(int n, SPH & sph){

    sph(0,0) = 0.505;
    sph(0,1) = 0.515;
    sph(0,2) = 0.51;
    sph(0,3) = 0.5;

    sph(1,0) = 0.5;
    sph(1,1) = 0.5;
    sph(1,2) = 0.45;
    sph(1,3) = 0.45;

    sph(2,0) = 0.0;
    sph(2,1) = 0.0; 
    sph(2,2) = 0.0; 
    sph(2,3) = 0.0;  

    sph(3,0) = 0.0;
    sph(3,1) = 0.0;
    sph(3,2) = 0.0;
    sph(3,3) = 0.0;

}

//Dam break 
void ic_dam_break(int n, SPH & sph){

    int el = pow(n,0.5);

    //Defining the initial distance between the particles in both directions
    double step = 0.19 / (el-1);

    //Starting position in x
    double posx = 0.01;
    double posy;

    //Assing tha values in x for all particles
    for (int i = 0; i < el; i++){

        for (int j = 0; j < el; j++){

            sph(0,i*el + j) = posx + double(rand())/ RAND_MAX/100000;
            sph(2,i*el + j) = 0.0;
        }

        posx += step;
                
    }

    //For uniform distribution the step in y has to be equal to the step in x
    step = 0.19 / (el-1);

    //Assing tha values in y for all particles
    for (int i = 0; i < el; i++){

        posy = 0.01;
        for (int j = 0; j < el; j++){

            sph(1,i*el + j) = posy + double(rand())/ RAND_MAX/100000;
            sph(3,i*el + j) = 0.0 ;
            posy += step;
            
        }      

    }

}

//Block Drop
void ic_block_drop(int n, int n1, int n2,SPH & sph){

    //Defining the distance between neighbouring particles in x and y
    //0.2 is the total distance in x and 0.3 in y
    double dx = 0.2 / double((n1-1));
    double dy = 0.3 / double((n2-1));

    //Starting position in x
    double posx = 0.1;
    double posy;
    int kx,ky;

    //Assing tha values in x for all particles
    for (int i = 0; i < n1; i++){

        for (int j = 0; j < n2; j++){

            kx = i*n2 + j;
            sph(0,kx) = posx + double(rand())/ RAND_MAX/100000;
            sph(2,kx) = 0.0;

        }

        posx += dx;
                
    }

   
    //Assing tha values in y for all particles
    for (int i = 0; i < n1; i++){

        posy = 0.3;
        for (int j = 0; j < n2; j++){

            ky = i*n2 + j;

            sph(1,ky) = posy + double(rand())/ RAND_MAX/100000;
            sph(3,ky) = 0.0 ;
            posy += dy;
            
        }      

    }

}

//Droplet
void ic_droplet(int n, SPH & sph){
    
    double * xi = new double[n];
    double * yi = new double[n];
    int el = pow(n,0.5);
    int kx;

    //For uniform distribution the step in y has to be equal to the step in x
    double step = 0.2 / (el-1);

    double posx = 0.4;    //Starting position in x
    double posy;          //Starting position in y

    for (int i = 0; i < el; i++){

        for (int j = 0; j < el; j++){

            xi[i*el + j] = posx ;
            
        }

        posx += step;
                
    }

    step = 0.2 / (el-1);

    for (int i = 0; i < el; i++){

        posy = 0.6;
        for (int j = 0; j < el; j++){

            yi[i*el + j] = posy ;
            
            posy += step;
            
        }      

    }
    kx = 0;
    for (int i = 0; i < el; i++){

        for (int j = 0; j < el; j++){

          if(sqrt(pow((yi[i*el + j]-0.7),2)+pow((xi[i*el + j]-0.5),2)) <= 0.1){

                
                sph(0,kx) = xi[i*el + j]+ double(rand())/ RAND_MAX/100000; 
                sph(1,kx) = yi[i*el + j]+ double(rand())/ RAND_MAX/100000;
                sph(2,kx) = 0; 
                sph(3,kx) = 0;
                kx++;
             }

            else{

                }
        }      

    }

    delete[] xi;
    delete[] yi;
}

//This function has the aim of defining the number of particles that will be in the circular region
int dropletn(int n){

    //The bellow porcess is similar to the dam break and creates an initial square
    double * xi = new double[n];
    double * yi = new double[n];
    int el = pow(n,0.5);
    double step = 0.2 / (el-1);
    double posx = 0.4;
    double posy;

    for (int i = 0; i < el; i++){

        for (int j = 0; j < el; j++){

            xi[i*el + j] = posx;
            
        }

        posx += step;
                
    }

    step = 0.2 / (el-1);

    for (int i = 0; i < el; i++){

        posy = 0.6;
        for (int j = 0; j < el; j++){

            yi[i*el + j] = posy;
            
            posy += step;
            
        }      

    }

    //After the initial square was created, the number of particles that are in that square and from a distance 
    //from the centre less or equal to the radius of the circle is calculated
    int count = 0;
    for (int i = 0; i < el; i++){

        for (int j = 0; j < el; j++){

           if(sqrt(pow((yi[i*el + j]-0.7),2)+pow((xi[i*el + j]-0.5),2)) <= 0.1){

               count ++;
           }

           else{

               count += 0;
           }
        }      

    }
    
    delete[] xi;
    delete[] yi;

    return count;
}