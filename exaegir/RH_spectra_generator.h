// -*- mode: c++ ; -*-
/** \file exaegir/RH_spectra_generator.h
 * Author(s) : Francois Mauger <mauger@lpccaen.in2p3.fr>
 * Creation date: 2023-09-21
 * Last modified: 2023-09-21
 */

#ifndef EXAEGIR_RH_SPECTRA_GENERATOR_H
#define EXAEGIR_RH_SPECTRA_GENERATOR_H 1

// - Bayeux:
#include <genbb_help/i_genbb.h>
#include <datatools/utils.h>
#include <datatools/clhep_units.h>

// - ROOT:
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TRandom3.h> 
#include <TMath.h> 
#include <TVector3.h>

namespace exaegir {

  struct momentum_vector
  {   
    double pMag;
    double px;
    double py;
    double pz;
    double phi;
    double theta;
  };

  class RH_spectra_generator
    : public ::genbb::i_genbb
  {

  public:

    /// Constructor
    RH_spectra_generator();

    /// Destructor
    virtual ~RH_spectra_generator();

    /// Check if an external random engine can be plugged in
    bool can_external_random() const override;

    /// Check initialization status
    bool is_initialized() const override;

    /// Initialize the generator from configuration properties
    void initialize(const datatools::properties & config,
		    datatools::service_manager & srv_mgr_,
		    genbb::detail::pg_dict_type & pg_dict_) override;

    /// Reset the object
    void reset() override;

    /// Check is a next event is available
    bool has_next() override;

    double get_bin_width(TTree *t, double &e2);
    double get_pMag(double T, double m);
    double sample_phi();
    double sample_theta();
    double sample_theta_diff(double k);
    double get_px(double pMag, double phi, double theta);
    double get_py(double pMag, double phi, double theta);
    double get_pz(double pMag, double theta);
    momentum_vector get_first_vector(double T, double m);
    momentum_vector get_pPrime(double T, double m, double k,momentum_vector p1);
    momentum_vector get_second_vector(double T, double m, double k, momentum_vector p1);

 
  protected:

    /// Load/generate a new event
    void _load_next(::genbb::primary_event & event_,
		    bool compute_classification_ = true) override;

  private:
  
    void _at_init_();

    void _at_reset_();

    double get_bin_width(TTree *t);

  private:
      
    bool _initialized_ = false; //!< Initialization flag

    // Configuration parameters:
    const double ELECTRON_MASS = 0.5109989461;
    std::string _spectrum_input_path;
    double _angular_correlation;


    // Working data:
    TH2D *h2d_spectrum;

    GENBB_PG_REGISTRATION_INTERFACE(RH_spectra_generator)
      
  };
      
} // end of namespace exaegir

#endif // EXAEGIR_RH_SPECTRA_GENERATOR_H
