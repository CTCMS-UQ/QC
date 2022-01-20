# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 08:55:41 2022

@author: uqhhoh
"""

import pandas as pd
import matplotlib.pyplot as plt


#Units output by program is quantum capacitance per supercell (C^2/eV·supercell)					
#1 eV =	1.60E-19 J		
#1 F  = 1 C2/J		
#Lattice constants of supercell
#a=2.467724085 #convert angstroms to cm
#b=2.467724085
#c=20

"""
#For orthorrhombic cells
#1 supercell = a*b	Å^2 = ab*10e-16 cm^2 = cross-section area
# supercell=a*b
#print(supercell)
"""				

#1 C^2/eV·supercell= 6.2415e+40 divide by supercell µF/cm^2
constant=6.2415e+40
supercell=20.9775518 #cross-section area, Volume/c

qc_data = "C:/Users/uqhhoh/OneDrive - The University of Queensland/Desktop/QH_QC_5V.dat"

df = pd.read_csv(qc_data,sep='\s+',skiprows=(1))


x = df.loc[:99,'#Phi']
raw_y = df.loc[:99,"QC"]
y = (raw_y)*constant/(supercell)


# Set a large font size for DOS plots for figures/papers
font = {'family' : 'arial',
    'weight' : 'bold',
    'size'   : 16}

plt.rc('font', **font)


plt.plot(x,y)
plt.xlabel('Potential (V)')
plt.ylabel('Quantum Capacitance $\mathregular{(µF/cm^2)}$')  
# need to change units $\mathregular{µF/cm^2}$
#plt.title('Simple Line Plot')
plt.xlim(-2,2)
ymin,ymax = plt.ylim(0,9000)
plt.ylim(ymin,ymax)

# Add a vertical line to mark the zero potential    
plt.vlines(0,ymin,ymax,linestyle='--',color='black')