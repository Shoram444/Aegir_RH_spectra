// exaegir/RH_spectra_generator.cc

// Ourselves
#include <exaegir/RH_spectra_generator.h>

// Standard library:
#include <cmath>

// Third party:
// - Bayeux:
#include <datatools/exception.h>
#include <datatools/logger.h>
#include <datatools/properties.h>
#include <datatools/units.h>
#include <genbb_help/primary_event.h>
#include <genbb_help/pdg_particle_tools.h>
#include <genbb_help/kinematics.h>
#include <mygsl/rng.h>


// - ROOT:
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TRandom3.h> 
#include <TMath.h> 
#include <TVector3.h>

namespace exaegir {

  GENBB_PG_REGISTRATION_IMPLEMENT(RH_spectra_generator,
				  "exaegir::RH_spectra_generator")
    
  RH_spectra_generator::RH_spectra_generator()
  {
    _initialized_ = false;
    _at_reset_();
    return;
  }

  RH_spectra_generator::~RH_spectra_generator()
  {
    if (_initialized_) {
      reset();
    }
    return;
  }

  bool RH_spectra_generator::can_external_random() const
  {
    return true;
  }

  bool RH_spectra_generator::is_initialized() const
  {
    return _initialized_;
  }

  void RH_spectra_generator::initialize(const datatools::properties & config_,
					  datatools::service_manager & srv_mgr_,
					  genbb::detail::pg_dict_type & pg_dict_)
  {
    DT_THROW_IF(
      _initialized_, std::logic_error,
		  "Operation prohibited! Object is already initialized!"
    );

    _initialize_base(config_);

    _spectrum_input_path = config_.fetch_string("spectrum_input_path");

    DT_LOG_DEBUG(
      get_logging_priority(),
		  "Tabulated kinetic energy spectrum file path: '" << _spectrum_input_path << "'"
    );

    _angular_correlation = config_.fetch_real("angular_correlation");
   
    _at_init_();
    _initialized_ = true;
    return;
  }

  void RH_spectra_generator::reset()
  {
    if (!_initialized_) {
      return;
    }
    _initialized_ = false;
    _at_reset_();
    return;
  }

  bool RH_spectra_generator::has_next()
  {
    return true;
  }

  void RH_spectra_generator::_load_next(
    ::genbb::primary_event & event_,
		bool compute_classification_)
  {
    DT_LOG_TRACE_ENTERING(get_logging_priority());

    double T1, T2; // kinetic energy of the two electrons
    h2d_spectrum->GetRandom2(T1, T2);

    TVector3 v1, v2;
    momentum_vector p1, p2;

    p1 = get_first_vector(T1, ELECTRON_MASS);
    p2 = get_second_vector(T2, ELECTRON_MASS, _angular_correlation, p1);  

    event_.reset();
    event_.set_time(0.0 * CLHEP::second);
    event_.set_label("RH");

    // first electron
    ::genbb::primary_particle & electron1 = event_.add_particle();
    electron1.set_type(::genbb::primary_particle::particle_type::ELECTRON);
    electron1.set_pdg_code(::genbb::pdg::particle::ELECTRON);
    electron1.set_time(0.0 * CLHEP::second);

    electron1.grab_momentum().setX(p1.px * CLHEP::MeV);
    electron1.grab_momentum().setY(p1.py * CLHEP::MeV);
    electron1.grab_momentum().setZ(p1.pz * CLHEP::MeV);

    // // second electron
    ::genbb::primary_particle & electron2 = event_.add_particle();
    electron2.set_type(::genbb::primary_particle::particle_type::ELECTRON);
    electron2.set_pdg_code(::genbb::pdg::particle::ELECTRON);
    electron2.set_time(0.0 * CLHEP::second);

    electron2.grab_momentum().setX(p2.px * CLHEP::MeV);
    electron2.grab_momentum().setY(p2.py * CLHEP::MeV);
    electron2.grab_momentum().setZ(p2.pz * CLHEP::MeV);

    if (compute_classification_) {
      event_.compute_classification();
    }
    
    DT_LOG_TRACE_EXITING(get_logging_priority());
    return;
  }
 
  void RH_spectra_generator::_at_init_()
  {
    // Read input spectrum into a ttree, the input file must be delimited with a single white-space " "
    TTree *t = new TTree("t", "tree");
    std::cout << "_spectrum_input_path = " << _spectrum_input_path << std::endl;

    t->ReadFile(_spectrum_input_path.c_str(), "e1/D:e2/D:rate/D");

    double e1, e2, rate;
    int n_entries = t->GetEntries();

    t->SetBranchAddress("e1", &e1);
    t->SetBranchAddress("e2", &e2);
    t->SetBranchAddress("rate", &rate);

    // convert the input spectrum into a 2D histogram (for this we need: xmin, xmax, ymin, ymax, nbinsx, nbinsy)
    double e1_min = t->GetMinimum("e1");  
    double e1_max = t->GetMaximum("e1"); 
    double e2_min = t->GetMinimum("e2"); 
    double e2_max = t->GetMaximum("e2"); 
    double bin_width = get_bin_width(t, e2);

    std::cout << "bin_width = " << bin_width << std::endl;

    double n_bins_e1 = (e1_max - e1_min) / bin_width +1 ;
    double n_bins_e2 = (e2_max - e2_min) / bin_width +1 ;

    double half_bin = bin_width / 2.0; 
    std::cout << "n_bins_e1 = " << n_bins_e1 << std::endl;
    std::cout << "n_bins_e2 = " << n_bins_e2 << std::endl;

    h2d_spectrum = new TH2D(
        "h2", 
        "h2", 
        n_bins_e1, 
        e1_min - half_bin, 
        e1_max + half_bin, 
        n_bins_e2, 
        e2_min - half_bin, 
        e2_max + half_bin
    );

    for( int i = 0; i < n_entries; i++)
    { 
        t->GetEntry(i); 
        h2d_spectrum->Fill(e1, e2, rate);
    }

    return;
  }
  
  void RH_spectra_generator::_at_reset_()
  {
    return;
  }


  double RH_spectra_generator::get_bin_width(TTree *t, double &e2)
  {
    // Computes the spacing used in the input spectra. This is done simply by subtracting the second row energy from first.
    // A uniform spacing is assumed!!!

      double first, bin_width;

      t->GetEntry(0);
      first = e2;
      t->GetEntry(1);

      bin_width = e2 - first;
      return bin_width;
  }

  double RH_spectra_generator::get_pMag(double T, double m)
  {
      //  Computes magnitude of momentum as:
      //     |p| = √((m + T)² - m²)
      //  
      //  Input arguments are:
      //  * double T  --> kinetic energy of the electron in [MeV]
      //  * double m  --> rest mass in [MeV/c²]

      return std::sqrt( std::pow(m+T ,2) - std::pow(m, 2));
  }

  double RH_spectra_generator::sample_phi()
  {
      //  Samples an azimuthal angle uniformly distributed on a sphere in [rad].

      TRandom3 r(0);
      return 2.0 * TMath::Pi() * r.Uniform();
  }

  double RH_spectra_generator::sample_theta()
  {
      //  Samples a polar angle uniformly distributed on a sphere in [rad].

      TRandom3 r(0);
      return TMath::ACos(1 - 2 * r.Uniform());
  }

  double RH_spectra_generator::sample_theta_diff(double k)
  {
      //  Samples theta_diff angle between the two electrons in [rad].
      //
      //  The pdf for the sampling is given as:
      //  pdf(x) = 0.5 * (1 - k*x) ; where x is subtited for cos(theta_diff)
      //
      //  The argument `k` is the angular correlation factor (i.e. for standard 2nubb: k = -0.88, for RH Se82 theoretical k = 0.37)
      // 
      //  Sampling itself is done using the Inverse CDF method 
      //  (Inverting leads to a quadratic formula, so the proper solution must be chosen, hence the if condition at the end.)

      TRandom3 r(0);

      double a = k/4.0;                   // coefficient near x^2
      double b = 0.5;                     // coefficient near x^1
      double c = b - a - r.Uniform();     // coefficient near x^0

      double d = std::sqrt(b*b - 4*a*c);  // discriminant
      double x1 = (-b - d) / (2 * a);
      double x2 = (-b + d) / (2 * a);

      if(std::abs(x1) <= 1.0) // check whether first solution is within bounds of acos
      {
          return TMath::ACos(x1);
      }
      else if(std::abs(x2) <= 1.0) // check whether second solution is within bounds of acos
      {
          return TMath::ACos(x2);
      }
      else
      {
          return -10000.0; // something fishy here
      }
  }

  double RH_spectra_generator::get_px( double pMag, double phi, double theta)
  {
      //  Computes x-component of the momentum as:
      //     pₓ = |p|*cos(ϕ)*sin(θ) 
      //
      //  Input arguments are:
      //  * pMag : momentum magnitude in [MeV]
      //  * phi  : azimuth angle in [rad]
      //  * theta: polar angle in [rad]

      return pMag * TMath::Cos(phi) * TMath::Sin(theta);    
  }

  double RH_spectra_generator::get_py( double pMag, double phi, double theta)
  {
      //  Computes y-component of the momentum as:
      //     p𝔂 = |p|*sin(ϕ)*sin(θ) 
      //
      //  Input arguments are:
      //  * pMag : momentum magnitude in [MeV]
      //  * phi  : azimuth angle in [rad]
      //  * theta: polar angle in [rad]

      return pMag * TMath::Sin(phi) * TMath::Sin(theta);    
  }

  double RH_spectra_generator::get_pz( double pMag, double theta)
  {
      //  Computes z-component of the momentum as:
      //     p𝐳 = |p|*cos(θ) 
      //
      //  Input arguments are:
      //  * pMag : momentum magnitude in [MeV]
      //  * theta: polar angle in [rad]

      return pMag * TMath::Cos(theta);    
  }

  momentum_vector RH_spectra_generator::get_first_vector(double T, double m)
  {
      //  Returns the momentum vector of the first electron. 
      //  Momentum magnitude is calculated from the electron kinetic energy (T). 
      //  The direction is generated uniformly on a sphere.
      //  
      //  The momentum components are (in order):
      //  * pMag - Magnitude of the momentum vector (given by the electron Energy)
      //  * px   - x-component
      //  * py   - y-component
      //  * pz   - z-component
      //  * θ    - polar angle  
      //  * ϕ    - azimuthal angle 


      double pMag = get_pMag(T, m);
      double phi = sample_phi();
      double theta = sample_theta();
      double px = get_px(pMag, phi, theta);
      double py = get_py(pMag, phi, theta);
      double pz = get_pz(pMag, theta);

      return momentum_vector({pMag, px, py, pz, phi, theta});
  }

  momentum_vector RH_spectra_generator::get_pPrime(double T, double m, double k,momentum_vector p1)
  {
      //  Returns the momentum (prime) vector of the second electron. The vector direction is generated as a random point from a circle projected onto the x-y plane. 
      //  The circle radius is given by r = sin(θdif)|p1|, where θdif is the generated angle between the two electrons and |p1| 
      //  is the magnitude of the momentum vector of the first electron. The process assumes p1 to be directed along the z-axis (hence projection onto x-y plane).

      double pMag = get_pMag(T, m);
      double phi = sample_phi();
      double theta = sample_theta_diff(k);
      double px = get_px(pMag, phi, theta);
      double py = get_py(pMag, phi, theta);
      double pz = get_pz(pMag, theta);

      return momentum_vector({pMag, px, py, pz, phi, theta});
  }

  momentum_vector RH_spectra_generator::get_second_vector(double T, double m, double k, momentum_vector p1)
  {
      //  Returns the momentum vector of the second electron. Momentum magnitude is calculated from the electron kinetic energy (T). 
      //  The direction is generated from sampling θdif angle (see `sample_theta_dif`). In order to calculate the direction: first a p2Prime vector is generated
      //  by sampling uniformly from a circle projected onto x-y plane, where the circle is determined by 
      //  the p1 vector directed along z-axis. Afterwards, p2Prime is rotated first along the y-axis, second along the z-axis to obtain p2 in the direction 
      //  of a cone around p1, where the cone is given by θdif. The equation for p2 is:
      //
      //  p2 = Rz*Ry*pPrime
      momentum_vector pPrime = get_pPrime(T, m, k, p1);
      
      double pMag = pPrime.pMag;
      double px = TMath::Cos(p1.phi) * (pPrime.px * TMath::Cos(p1.theta) + pPrime.pz * TMath::Sin(p1.theta)) - pPrime.py * TMath::Sin(p1.phi);
      double py = TMath::Sin(p1.phi) * (pPrime.px * TMath::Cos(p1.theta) + pPrime.pz * TMath::Sin(p1.theta)) + pPrime.py * TMath::Cos(p1.phi);
      double pz = pPrime.pz * TMath::Cos(p1.theta) - pPrime.px * TMath::Sin(p1.theta);

      double phi = TMath::ATan2(py, px);
      double theta = TMath::ACos(pz / pMag);

      return momentum_vector({pMag, px, py, pz, phi, theta});
  }



} // end of namespace exaegir




