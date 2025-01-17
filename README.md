# Aegir RH Spectra event generator

# INTRODUCTION

This plug-in can be used to simulate the two-neutrino double beta decay ($\nu_R\nu_L\beta\beta$) with **one** of the neutrinos being Right-Handed ($\nu_R$) for Se-82 isotope. See [here](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.171801) for reference.

The input spectrum (theoretical) was provided by L. Graf. This spectrum can be found in `resources/data/RH_input_data.data` file. The file is a standard `CSV` file with three columns separated by a **single white space**. The three columns represent the kinetic energy of the first electron $T_1$ and the second electron $T_2$, respectively. The last colum represents the differential phase-space factor $\frac{dG}{dT_1dT_2}$. 

The provided spectrum uses a 10 keV spacing between the points. This should be enough for SuperNEMO usecase, but may be improved upon request. 

# SAMPLING THE ENERGY AND ANGLE

In order to provide the required input for Falaise simulation - the initial momentum vector of the simulated particles - two components are required. First, the energy is of the electron (pair of electrons) is sampled using the provided input spectrum. This is done via converting the spectrum to a `TH2D` histogram which has a method `GetRandom2` which draws a sample of 2 random numbers (the kinetic energies of the two electrons). Second, the angular distribution of the decay angle **between the two electrons** is sampled using a three step approach:
1. A random direction of the **first** electron is sampled uniformly on a sphere (isotropic direction).
2. A $\theta_{diff}$ angle is sampled using the inverted CDF method of the angular distribution pdf provided in [here](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.171801). The shape of the pdf is: $pdf(\theta) \approx 1 - k\cos{\theta}$; **where $k$ is the angular correlation factor -> The free parameter of this process which has to be fitted from measured data.** A default value of $k=0.37$ can be used as per the supplemental material of the publication. 
3. From the cone created by drawing a dircle around the first electron's direction, where the size of the cone is given by $\theta_{diff}$ a random direction is chosen as the direction of the **second** electron. 

# INSTALLATION AND USAGE
(There's a lot of paths that need to be specified, this can probably be done simpler and better but I haven't bothered :( )

0. (Optional, if you haven't done this already :-) ) Load Falaise v5.x: `source /sps/nemo/sw/snswmgr/snswmgr.conf && snswmgr_load_setup falaise@5.1.2`
1. Build the library using: `./build.bash`
2. In the `simu.conf` file that you use for `flsimualte` application, you must specify the path to `lib64`, an example is shown in the `Example_simu.conf` example file in this repository.
3. In the `generators.conf` file you **must specify the path to the input data and desired k factor!!** An example is shown in the `Example_generators.conf` file.
4. In the `variant.profile` file you must specify the path to `generators.conf` file and the name of the generator "RH". Again an example is provided in `Example_variant.profile` file.  
5. Use `flsimulate` as usual with the config of `simu.conf`. That's it! :) 

(Optionally the `variant.profile` can be created using the `flsimulate-configure` application. Here, in the **Primary events** tab you must select the `aegir` "value". In the "generators_file" row you need to specify the path to the `generators.conf` file which includes the path to the input data and k factor. In the "selected" row you need to specify the name of the generator "RH".) 