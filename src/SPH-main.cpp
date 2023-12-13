// This is a main program to run the SPH simulation

#include "class.h"
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>

namespace po = boost::program_options;

/**Functions for the validation cases
 *(explantations are provided at their implementation bellow the main programme)
 **/

void ic_one_particle(int n, SPH &sph);

void ic_two_particles(int n, SPH &sph);

void ic_three_particles(int n, SPH &sph);

void ic_four_particles(int n, SPH &sph);

void ic_dam_break(int n, SPH &sph);

void ic_block_drop(int n, int n1, int n2, SPH &sph);

void ic_droplet(int n, SPH &sph);

int dropletn(int n);

// Start of the main programme
int main(int argc, char *argv[]) {

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
    return 1;
  }

  po::notify(vm);

  int n;        // number of particles
  int n1 = 17;  // required for ic-block-drop
  int n2 = 25;  // required for ic-block-drop
  int n3 = 400; // required for ic-dam-break and ic-droplet

  double T1 = vm["T"].as<double>();  // Total integration time
  double dt = vm["dt"].as<double>(); // Time step dt
  double h = vm["h"].as<double>();   // Radius of influence

  int T = int(T1 / dt) + 1; // Transform time in seconds to iterations

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
      return 1;
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

  // Calculate the mass of the partciles
  sph.mass();

  std::ofstream vOut("Positions-x-y.txt", std::ios::out | std::ios::trunc);
  std::ofstream vOut2("Energy-File.txt", std::ios::out | std::ios::trunc);

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

  // Time integration loop
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
  return 0;
}

// ========== Initial Conditions ==========

void ic_one_particle(int n, SPH &sph) {

  sph(0, 0) = 0.5;
  sph(1, 0) = 0.5;
  sph(2, 0) = 0.0;
  sph(3, 0) = 0.0;
}

void ic_two_particles(int n, SPH &sph) {

  sph(0, 0) = 0.5;
  sph(0, 1) = 0.5;

  sph(1, 0) = 0.5;
  sph(1, 1) = 0.01;

  sph(2, 0) = 0.0;
  sph(2, 1) = 0.0;

  sph(3, 0) = 0.0;
  sph(3, 1) = 0.0;
}

void ic_three_particles(int n, SPH &sph) {

  sph(0, 0) = 0.5;
  sph(0, 1) = 0.495;
  sph(0, 2) = 0.505;

  sph(1, 0) = 0.5;
  sph(1, 1) = 0.01;
  sph(1, 2) = 0.01;

  sph(2, 0) = 0.0;
  sph(2, 1) = 0.0;
  sph(2, 2) = 0.0;

  sph(3, 0) = 0.0;
  sph(3, 1) = 0.0;
  sph(3, 2) = 0.0;
}

void ic_four_particles(int n, SPH &sph) {

  sph(0, 0) = 0.505;
  sph(0, 1) = 0.515;
  sph(0, 2) = 0.51;
  sph(0, 3) = 0.5;

  sph(1, 0) = 0.5;
  sph(1, 1) = 0.5;
  sph(1, 2) = 0.45;
  sph(1, 3) = 0.45;

  sph(2, 0) = 0.0;
  sph(2, 1) = 0.0;
  sph(2, 2) = 0.0;
  sph(2, 3) = 0.0;

  sph(3, 0) = 0.0;
  sph(3, 1) = 0.0;
  sph(3, 2) = 0.0;
  sph(3, 3) = 0.0;
}

void ic_dam_break(int n, SPH &sph) {

  int el = pow(n, 0.5);
  // Initial distance between the particles in both directions
  double step = 0.19 / (el - 1);
  // Starting position in x
  double posx = 0.01;
  double posy;
  // Assing the values in x for all particles
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      sph(0, i * el + j) = posx + double(rand()) / RAND_MAX / 100000;
      sph(2, i * el + j) = 0.0;
    }
    posx += step;
  }

  // For uniform distribution the step in y has to be equal to the step in x
  step = 0.19 / (el - 1);

  // Assing values in y for all particles
  for (int i = 0; i < el; i++) {
    posy = 0.01;
    for (int j = 0; j < el; j++) {
      sph(1, i * el + j) = posy + double(rand()) / RAND_MAX / 100000;
      sph(3, i * el + j) = 0.0;
      posy += step;
    }
  }
}

void ic_block_drop(int n, int n1, int n2, SPH &sph) {

  // Distance between neighbouring particles in x and y
  // 0.2 is the total distance in x and 0.3 in y
  double dx = 0.2 / double((n1 - 1));
  double dy = 0.3 / double((n2 - 1));

  // Starting position in x
  double posx = 0.1;
  double posy;
  int kx, ky;

  // Assing the values in x for all particles
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      kx = i * n2 + j;
      sph(0, kx) = posx + double(rand()) / RAND_MAX / 100000;
      sph(2, kx) = 0.0;
    }
    posx += dx;
  }

  // Assing the values in y for all particles
  for (int i = 0; i < n1; i++) {
    posy = 0.3;
    for (int j = 0; j < n2; j++) {
      ky = i * n2 + j;
      sph(1, ky) = posy + double(rand()) / RAND_MAX / 100000;
      sph(3, ky) = 0.0;
      posy += dy;
    }
  }
}

// Droplet
void ic_droplet(int n, SPH &sph) {

  double *xi = new double[n];
  double *yi = new double[n];
  int el = pow(n, 0.5);
  int kx;

  // For uniform distribution the step in y has to be equal to the step in x
  double step = 0.2 / (el - 1);
  double posx = 0.4; // Starting position in x
  double posy;       // Starting position in y

  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      xi[i * el + j] = posx;
    }
    posx += step;
  }

  step = 0.2 / (el - 1);

  for (int i = 0; i < el; i++) {
    posy = 0.6;
    for (int j = 0; j < el; j++) {
      yi[i * el + j] = posy;
      posy += step;
    }
  }
  kx = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (sqrt(pow((yi[i * el + j] - 0.7), 2) +
               pow((xi[i * el + j] - 0.5), 2)) <= 0.1) {
        sph(0, kx) = xi[i * el + j] + double(rand()) / RAND_MAX / 100000;
        sph(1, kx) = yi[i * el + j] + double(rand()) / RAND_MAX / 100000;
        sph(2, kx) = 0;
        sph(3, kx) = 0;
        kx++;
      }
    }
  }
  delete[] xi;
  delete[] yi;
}

// Defines the number of particles that will be in the circular region
int dropletn(int n) {

  // Process similar to dam break. Creates an initial square
  double *xi = new double[n];
  double *yi = new double[n];
  int el = pow(n, 0.5);
  double step = 0.2 / (el - 1);
  double posx = 0.4;
  double posy;

  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      xi[i * el + j] = posx;
    }
    posx += step;
  }

  step = 0.2 / (el - 1);

  for (int i = 0; i < el; i++) {
    posy = 0.6;
    for (int j = 0; j < el; j++) {
      yi[i * el + j] = posy;
      posy += step;
    }
  }

  // After the initial square is created, the number of particles that are in
  // that square and from a distance from the centre less or equal to the radius
  // of the circle is calculated
  int count = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (sqrt(pow((yi[i * el + j] - 0.7), 2) +
               pow((xi[i * el + j] - 0.5), 2)) <= 0.1) {
        count++;
      } else {
        count += 0;
      }
    }
  }
  delete[] xi;
  delete[] yi;

  return count;
}