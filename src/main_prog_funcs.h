// Definition of the functions called by the main program
#ifndef MAIN_H
#define MAIN_H
#include <boost/program_options.hpp>
#include <memory>

#include "fluid.h"
#include "sph_solver.h"

namespace po = boost::program_options;

void initialise(std::unique_ptr<Fluid>& fluidPtr, SphSolver& solver);

std::tuple<std::ofstream, std::ofstream, std::ofstream> initOutputFiles(
    const std::string& folderPath);

void retrieveInputsFromFile(const std::string& fileName,
                            const std::string& icCase,
                            const po::options_description& desc,
                            po::variables_map& vm);

void handleInputErrors(const po::variables_map& caseVm,
                       const po::variables_map& domainVm,
                       const po::variables_map& constantsVm,
                       const po::variables_map& icVm);

void setInitialConditions(const std::string& icCase,
                          std::unique_ptr<Fluid>& fluidPtr,
                          const po::variables_map& icVm,
                          const po::variables_map& domainVm);

void createDirectory(std::string folderPath);

void storeToFile(Fluid& fluid, std::string type, std::ofstream& targetFile,
                 double dt = 0.0, double currentTime = 0.0);
#endif