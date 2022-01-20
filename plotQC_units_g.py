# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 19:37:10 2022

@author: uqhhoh
"""


import pandas as pd
import matplotlib.pyplot as plt


#Units output by program is quantum capacitance per supercell (C^2/eV·supercell)					
#1 eV =	1.60E-19 J		
#1 F  = 1 C2/J		
#1 C^2/eV·supercell= 3.75867E+42 divide by M  F/g
constant=3.75867E+42
# 1 supercell = M g/mol = Mr/6.022*10^23 g
M = 24 #2 C atoms
#M= 169  #1H , 14C atoms = 169 amu

qc_data = "C:/Users/uqhhoh/OneDrive - The University of Queensland/Desktop/QH_QC_5V.dat"
#qc-data = "C:/Users/uqhhoh/OneDrive - The University of Queensland/Documents/UQ/Magnus data/test_g1x1/kpts21/g1x1_QC.dat"

df = pd.read_csv(qc_data,sep='\s+',skiprows=(1))


x = df.loc[:99,'#Phi']
raw_y = df.loc[:99,"QC"]
y = (raw_y)*constant/M


# Set a large font size for DOS plots for figures/papers
font = {'family' : 'arial',
    'weight' : 'bold',
    'size'   : 16}

plt.rc('font', **font)


plt.plot(x,y)
plt.xlabel('Potential (V)')
plt.ylabel('Quantum Capacitance $\mathregular{(F/g)}$')  
# need to change units $\mathregular{µF/cm^2}$
#plt.title('Simple Line Plot')
plt.xlim(-2,2)
ymin,ymax = plt.ylim(0,)
plt.ylim(ymin,ymax)

# Add a vertical line to mark the zero potential    
plt.vlines(0,ymin,ymax,linestyle='--',color='black')