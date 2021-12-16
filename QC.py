#!/usr/bin/python3
import numpy as np
import argparse
import linecache

class Units:
    """ Bag of data containing units and parameters for the simulation system.

    This class uses the Borg pattern to ensure that it can be instantiated multiple times with each
    instant containing the same values for all the constants."""
    # Shared state for simulation parameters
    _params = {}

    # These are constants of the universe and can never change
    kb_eV = 8.617333262e-5 # Boltzmann's constant
    e_charge = 1.602176634e-19 # Elementary charge
    mF_per_cm2_constant = 6.241418e40 # Constant to convert default units to mF/cm^-2
    F_per_gram_constant = 3.75867e42 # Constant to convert default units to F/g
    def __init__(self):
        self.__dict__ = self._params

class DOS:
    """Top-level class for reading input and manipulating density of states data. Can be overridden to 
    parse different file formats."""
    def __init__(self, infile, debug = False):
        # Read the data file into a set of three numpy arrays. We need to pass "unpack = True" so that
        # the data comes out in the right order for the comma operator unpacking
        self.read_input(infile)

        # Now sum the absolute values of the DOS terms to get the total density of states
        self.TDOS = abs(self.up) + abs(self.down)

        # Keep track of whether to print debugging output
        self.debug = debug
    
    def read_input(self, infile):
        self.E, self.up, self.down = np.loadtxt(infile, unpack=True)

    def print_QC(self, phi_min, phi_max, num_points = 100):
        """ Print the Quantum capacitance for external potentials in the range [phi_min, phi_max].

        This function builds an array of linearly spaced points between phi_min and phi_max. If no value
        is provided for num_points then this function defaults to calculating for 100 different
        potentials.
        """

        units = Units()
        # Print calculation parameters and header
        print(f"#T = {units.T:.3f} K, area = {units.S:.3f} Angstrom^2, mass = {units.M:.3f} amu")
        print("#Phi(eV)    QC(C^2 eV^-1 cell^-1)    QC(mF cm^-2)     QC(F g^-1)")

        # Make the array of potentials
        phispace = np.linspace(phi_min, phi_max, num_points)

        for phi in phispace:
            self.integrate(phi)

    def integrate(self, phi):
        """ Calculates the Fermi-Dirac distribution function this class's energy distribution. 
        
        This function creates a new array by evaluating the Fermi-Dirac function for all elements in 
        self.E"""

        units = Units()
        etherm = (self.E - phi)/(2*units.kbT)
        # Now mask off any values too small to be useful. Values of etherm < -330 or > 330 cause
        # FD ~= sech^2(etherm) to underflow. Fortunately, their contribution is so small that 
        # we can remove them from the calculation without loss of accuracy.
        # NOTE: It's feasible that there are some situations in which this will result in a loss of
        # accuracy. Consider the case where we have a large positive or negative etherm, resulting in a
        # very small value for the Fermi-Dirac distribution, but a very large density of states.
        # Analytically, these two values would produce a reasonably sized integrand when multiplied
        # together, but this masking procedure would give a value of zero. It's extremely unlikely that
        # this would occur, as it would require an extraordinarily large density of states (which is
        # likely to be unphysical), but it's probably worth keeping in mind.
        #
        # TODO: should probably refine these bounds later.
        etherm_masked = np.ma.masked_outside(etherm, -330, 330)

        # Get the mask and use it on the other arrays
        mask = etherm_masked.mask
        # Print the number of masked-off values if debugging is enabled
        if self.debug:
            print("#Number of underflowing elements (Fermi-Dirac function too small to represent):")
            print(f"#{mask.size}")

        TDOS_masked = np.ma.masked_array(self.TDOS, mask=mask)
        E_masked = np.ma.masked_array(self.E, mask=mask)

        # Now we can do the calculation (using a hyperbolic trig identity, since numpy doesn't have a
        # sech function)
        FD = 2/(np.cosh(2*etherm_masked) + 1)

        # Calculate the integrals using the trapezoidal method
        integral = np.trapz(FD*TDOS_masked, x = E_masked)
        if self.debug:
            print(f"#Unnormalised integral = {integral}")
        
        # Finally, print the output in multiple different units.
        # First, get the "default" units by multiplying by e^2/kbT
        integral *= units.e_charge**2/(4*units.kbT)

        # Now get the value of the integral in mF cm^-2
        integral_mF_per_cm2 = integral * units.mF_per_cm2_constant/units.S
        # And in F/g
        integral_F_per_gram = integral * units.F_per_gram_constant/units.M

        print(f"{phi}     {integral}     {integral_mF_per_cm2}     {integral_F_per_gram}")

class DOSCAR(DOS):
    def read_input(self, infile):
        # Read straight from the VASP output format (DOSCAR)
        # Get info about the number of states and Fermi level
        energy_params = linecache.getline(infile,6).split()
        num_states = int(energy_params[2])
        E_fermi = float(energy_params[3])

        self.E = np.zeros(num_states)
        self.up = np.zeros(num_states)
        self.down = np.zeros(num_states)

        # Loop through the NEDOS lines describing the ion
        # DOSCAR contains 7 lines before the total DOS is plotted
        for line in range(7,7+num_states):
            # E_index starts counting at 0 wherever we are in the file
            E_index = line - 7
            # Break the line up into spin-orbital components
            data = linecache.getline(infile,line).split()

            # Set zero of energy scale at the Fermi level by subtracting
            # Ef from each point
            self.E[E_index] = (float(data[0])) - E_fermi

            # Get total DOS at each energy
            DOSup = float(data[1])        
            self.up[E_index] = DOSup

            DOSdown = float(data[2])        
            self.down[E_index] = DOSdown
######################################################################################################
if __name__ == "__main__":

    # Initialise the command-line argument parser
    parser = argparse.ArgumentParser(description="A program to calculate quantum capacitance from \
density of states. \
The program reads input from a file: the spin-up and spin-down DOS in units of states eV-1 supercell-1 \
as a function of energy in eV, and outputs the external potential in eV and the \
quantum capacitance in various units")
    parser.add_argument("datafile", help="Name of file containing density of states data.")
    parser.add_argument("-T", "--temperature", help="Temperature of the capacitor (in Kelvin).",
                        type=float, default=300.0)
    parser.add_argument("-N", "--num-points", help="Number of points to calculate capacitance for.",
                        type=int, default=100)
    parser.add_argument("--phi", help="End points of the range of potentials to use (in eV). The \
    capacitance will be calculated for potentials in the range [-phi,phi].",
                        type=float, default=5.0)
    parser.add_argument("-S", "--surface-area", help="Cross-sectional area of the electrode (in Angstrom^2).",
                        type=float, default=5.241)
    parser.add_argument("-M", "--mass", help="Mass of the electrode (in amu).",
                        type=float, default=24.0)

    # This argument controls whether or not to print intermediate results (unnormalised integrals,
    # number of masked-out values)
    parser.add_argument("-d", "--debug", help="Print intermediate results of the calculation.",
                        action="store_true")

    # Parse and store arguments
    args = parser.parse_args()
    units = Units()

    # Initialise parameters and units
    units.T = args.temperature
    units.kbT = units.T*units.kb_eV
    units.S = args.surface_area
    units.M = args.mass

    # Define the range of potentials we want to run over
    phi_min = -1*args.phi
    phi_max = args.phi

    # First, read in the data file as numpy arrays
    data = DOSCAR(args.datafile, debug=args.debug)
    data.print_QC(phi_min, phi_max, args.num_points)
    
