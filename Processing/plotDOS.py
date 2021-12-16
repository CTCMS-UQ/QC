#!/usr/bin/python3
"""
Extract total DOS from spin-polarized DOSCAR
"""

import linecache

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Set a large font size for DOS plots for figures/papers

font = {'family' : 'arial',
    'weight' : 'bold',
    'size'   : 16}

plt.rc('font', **font)

#import sys

#file_location = 'C:/Users/uqhhoh/OneDrive - The University of Queensland/Desktop'
file_location = "./"
doscar = file_location + '/DOSCAR'
output = file_location + '/output.txt'
#ion_number = sys.argv[1]

#first_parameters = linecache.getline(doscar,1)

#n_ions,n_real_ions,PDOS_bool,ncdij = first_parameters.split()

name = linecache.getline(doscar,5)

energy_params = linecache.getline(doscar,6).split()

NEDOS = int(energy_params[2])
E_f = float(energy_params[3])


"""
Create lists of zeros to store data
"""
energies = [0] * NEDOS

total_DOSup= [0] * NEDOS
total_DOSdown= [0] * NEDOS


def get_ion_PDOS(file,points,energies,Ef):
# Loop through the NEDOS lines describing the ion
# DOSCAR contains 7 lines before the total DOS is plotted
    for E_point in range(7,7+NEDOS):
# E_index starts counting at 0 wherever we are in the file
        E_index = E_point - 7
# Break the line up into spin-orbital components
        all_PDOS = linecache.getline(file,E_point).split()

# Set zero of energy scale at the Fermi level by subtracting
# Ef from each point
        energies[E_index] = (float(all_PDOS[0])) - Ef

# Get total DOS at each energy
        DOSup = float(all_PDOS[1])        
        total_DOSup[E_index] += DOSup

        DOSdown = float(all_PDOS[2])        
        total_DOSdown[E_index] += DOSdown

    return np.asarray((energies,total_DOSup,total_DOSdown))

dos_array = get_ion_PDOS(doscar,NEDOS,energies,E_f)

def plot_total_DOSup(data_array,fermi_energy,E_lower,E_upper):
    # Set range of energy axis
    # Set energy data as a variable for neatness
    E = data_array[0]
    # Plot the DOS
    plt.plot(E,data_array[1])
    
    plt.xlim(E_lower,E_upper)
    ymin,ymax = plt.ylim(-20,20)
    
def plot_total_DOSdown(data_array,fermi_energy,E_lower,E_upper):
    # Set range of energy axis
    # Set energy data as a variable for neatness
    E = data_array[0]
    # Plot the DOS
    plt.plot(E,(data_array[2])*-1)
    
    plt.xlim(E_lower,E_upper)
    ymin,ymax = plt.ylim(-20,20)
    
# Add a dotted line across the plot to mark zero DOS    
#    ymin,ymax = plt.ylim()
    plt.hlines(0,ymin,ymax,linestyle='-',color='black')
    
# Add a line across the plot to mark the Fermi energy    
#    ymin,ymax = plt.ylim()
    plt.vlines(0,ymin,ymax,linestyle='--',color='black')

    # Labels, annotations etc    
    plt.xlabel('Energy / eV')
    plt.ylabel('DOS')
#    plt.annotate("$E_F$",xy=(fermi_energy,ymin),xytext=((fermi_energy-0.3),(ymin-0.15))  
    plt.show()

def print_to_file(data_array, filename):
    # Write the elements of a data array to a plain text file
    with open(filename, "w") as ofp:
        # First, write a header
        ofp.write("#Energy           TDOS-UP         TDOS-DOWN\n")
        # Now write the data
        for ii in range(data_array.shape[1]):
            # Change the "5" in braces to however many decimal places you want in the output file
            ofp.write("{:5f}        {:5f}        {:5f}\n".format(data_array[0][ii], data_array[1][ii], data_array[2][ii]))


print_to_file(dos_array, "output.txt")
plot_total_DOSup(dos_array,E_f,-2,2)
plot_total_DOSdown(dos_array,E_f,-2,2)
