// Standard library:
#include <iostream>
#include <cstdlib>

// Bayeux:
#include <datatools/clhep_units.h>
#include <genbb_help/primary_event.h>
#include <mygsl/rng.h>
#include <mygsl/histogram.h>

// This project:
#include <exaegir/RH_spectra_generator.h>

void run1();

int main(void)
{
  int error_code = EXIT_SUCCESS;
  run1();
  exit(EXIT_SUCCESS);
}

void run1()
{
  int error_code = EXIT_SUCCESS;

  // Random number generator:
  int32_t prngSeed(314159); 
  mygsl::rng prng(prngSeed);
  
  // Event generator configuration:
  datatools::properties EvGenConfig;
  
  // Read the event generator configuration file:
  std::string EvGenConfigPath = "@exaegir:config/config.conf";
  datatools::fetch_path_with_env(EvGenConfigPath);
  datatools::properties::read_config(EvGenConfigPath, EvGenConfig);
  {
    std::cerr << "[debug] Event generator configuration:" << std::endl;
    boost::property_tree::ptree pOptions;
    pOptions.put("indent", "[debug] ");
    EvGenConfig.print_tree(std::cerr, pOptions);
  }
  
  // The event generator:
  exaegir::RH_spectra_generator EvGen;
  EvGen.set_external_random(prng);

  // Initialize the event generator:
  EvGen.initialize_standalone(EvGenConfig);

  // Use the event generator:
  genbb::primary_event event;
  bool computeClassification = true;
  unsigned int nbEvents = 5;
  for (unsigned int iEvent = 0; iEvent < nbEvents; iEvent++) {
    std::clog << "Event #" << iEvent << std::endl;
    EvGen.load_next(event, computeClassification);
    event.print_tree(std::cerr);
  }

  // Terminate the event generator:
  EvGen.reset();
  
  exit(EXIT_SUCCESS);
}
