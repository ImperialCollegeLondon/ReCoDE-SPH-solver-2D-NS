// This is a main program to run the SPH simulation
#include "class.h"
#include "initial_conditions.h"
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>

namespace po = boost::program_options;

SPH initialise(int &n, int &T, double &h, double &dt);
void init_output_files(std::ofstream &vOut, std::ofstream &vOut2);
void time_integration(SPH &sph,int &n, int &T, double &h, double &dt,std::ofstream &vOut, std::ofstream &vOut2);
// Start of the main programme
int main(int argc, char *argv[]) {

  // Declare the parameters of the problem
  int n;
  int T;
  double h;
  double dt;

  // Read input files, initialise the sph class and the parameters of the problem
  SPH sph = initialise(n,T,h,dt);

  // Declare and initialise the output files
  std::ofstream vOut("Positions-x-y.txt", std::ios::out | std::ios::trunc);
  std::ofstream vOut2("Energy-File.txt", std::ios::out | std::ios::trunc);
  init_output_files(vOut,vOut2);

  // Time integration loop
  time_integration(sph,n,T,h,dt,vOut,vOut2);

  return 0;
}

SPH initialise(int &n, int &T, double &h, double &dt){

  // Process to obtain the directions provided by the user 
  po::options_description desc("Allowed options");
  desc.add_options()("init_condition", po::value<std::string>(),
                     "take an initial condition")("T", po::value<double>(),
                                                  "take integration time")(
      "dt", po::value<double>(), "take time-step")("h", po::value<double>(),
                                                   "take radius of influence");

  po::variables_map vm;
  std::ifstream inputFile;
  inputFile.open("../exec/inputs/case.txt");

  if (inputFile.is_open()) {
    po::store(po::parse_config_file(inputFile, desc), vm);
    inputFile.close();
  } else {
    std::cerr << "Error opening file: inputs.txt" << std::endl;
  }

  po::notify(vm);

  int n1 = 17;  // required for ic-block-drop
  int n2 = 25;  // required for ic-block-drop
  int n3 = 400; // required for ic-dam-break and ic-droplet

  double T1 = vm["T"].as<double>();  // Total integration time
  dt = vm["dt"].as<double>(); // Time step dt
  h = vm["h"].as<double>();   // Radius of influence

  T = int(T1 / dt) + 1; // Transform time in seconds to iterations

  std::map<std::string, int> initConditionToParticlesMap = {
      {"ic-one-particle", 1},      {"ic-two-particles", 2},
      {"ic-three-particles", 3},   {"ic-four-particles", 4},
      {"ic-dam-break", n3},        {"ic-block-drop", n1 * n2},
      {"ic-droplet", dropletn(n3)}};

  n = initConditionToParticlesMap[vm["init_condition"].as<std::string>()];


  // Define the solver object (called sph)
  // In its definition, the number of particles is required
  SPH sph(n);

  /**After the number of partciles is introduced inside the class and
   * therefore the appropriate matrices are initialized, the particles
   * are ordered in the correct positions
   **/

  // Create map to associate function names with function pointers
  std::map<std::string, std::function<void(int, SPH &)>> functionMap = {
      {"ic-one-particle", ic_one_particle},
      {"ic-two-particles", ic_two_particles},
      {"ic-three-particles", ic_three_particles},
      {"ic-four-particles", ic_four_particles},
      {"ic-dam-break", ic_dam_break},
      {"ic-droplet", ic_droplet}};

  // Get the function pointer from the map
  auto initFunc = functionMap.find(vm["init_condition"].as<std::string>());
  if (initFunc != functionMap.end()) {
    int n_particles = n;

    // The ic-droplet case requires a different n argument.
    if (vm["init_condition"].as<std::string>() == "ic-droplet") {
      n_particles = n3;
    }
    initFunc->second(n_particles, sph);
    sph >> dt;
    sph < h;
  } else {
    /**The ic-block-drop case is not in the map because it has two
     * additional parameters, so it requires a different case.
     **/
    if (vm["init_condition"].as<std::string>() == "ic-block-drop") {
      ic_block_drop(n, n1, n2, sph);
      sph >> dt;
      sph < h;
    } else {
      std::cerr << "Error: Function not found!" << std::endl;
    }
  }

    /** Split the original matrix through which the values of the
   * initial coordinates and the initial velocities were introduced
   * inside the class.
   **/
    sph.x0();
    sph.y0();
    sph.vx0();
    sph.vy0();

    // Calculate the mass of the particles
    sph.mass();
  
  return sph;
}

void init_output_files(std::ofstream &vOut, std::ofstream &vOut2){

  vOut.precision(5);
  vOut << "x"
       << "          "
       << "y"
       << "\n";
  vOut2.precision(5);
  vOut2 << "t"
        << "      "
        << "Ek"
        << "       "
        << "Ep"
        << "     "
        << "Etotal"
        << "\n";

}

void time_integration(SPH &sph,int &n, int &T, double &h, double &dt,std::ofstream &vOut, std::ofstream &vOut2){

  for (int t = 0; t < T; t++) {

    // Pass the specific "time" of the loop inside the class
    sph > t;

    // In each iteration the disatnces between the particles are recalculated,
    // as well as their densities
    sph.rVec();
    sph.den();
    sph.spatial();
    // sph.getdata();

    // Write energies on the Energy-File
    vOut2 << t * dt << "  " << sph.Ek() << "  " << sph.Ep() << "  "
          << sph.Ep() + sph.Ek() << "\n";

    // Get the posistions after integration is completed
    if (t == T - 1) {

      for (int l = 0; l < n; l++) {

        vOut << sph.retx(l) << " " << sph.rety(l) << "\n";
      }
    }
  }

}