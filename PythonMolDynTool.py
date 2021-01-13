###################### IMPORTING MODULES - MAKE SURE NUMPY, PANDAS, PYTRAJ AND MDANALYSES ARE DOWNLOADED ##############
import numpy as np
import pandas as pd
from numpy import ndarray
import MDAnalysis as mda
import MDAnalysis.analysis.rms
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis import align
from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)
import matplotlib.pyplot as plt

####################### DEFINING TRAJECTORIES, ATOMGROUPS, REFERENCE STRUCTURES ETC ###################################

BlaC = mda.Universe("BlaC_Structure.pdb","BlaC.dcd") #This step took a while in deciding which file types to use - MDAnalysis supports multiple file types but different file types are better suited to different analyses (e.g. .dcd file types can't seem to be read by MDAnalysis for charge types and thus can't be used for Hydrogen Bond analysis). Initially I tried to use a multi-structure PDB format for the trajectory, but running the code seemed to bring up a lot of errors - potentially this could have been due to two PDB file formats used. An NC file type for the trajectory also didn't seem to function properly in the code, so I settled for fa DCD file format.

protein= BlaC.select_atoms("protein")

reference_coordinates = BlaC.trajectory.timeseries(asel=protein).mean(axis=1)
reference = mda.Merge(protein).load_new(reference_coordinates[:, None, :], order="afc") #here a reference structure is built using the trajectory. This took a lot of time to correctly identify the arguments that were required for this - specifically understanding what to put for the .load_new(reference_coordinates... in order to create an average structure of the trajectory for RMSF analysis.

Calpha = BlaC.select_atoms('name CA') #defines the Calpha atoms as those atoms in the protein structure file and the trajectory file as CA
FULL = 'backbone and resid 1-265'  #In these two lines you are able to define part of the protein to be used in analyses later on

Nterm = BlaC.select_atoms('resid 1 and name N') #defines the N terminus of the protein as the N atom on Residue 1
Cterm = BlaC.select_atoms('resid 265 and name C') #defines the C terminus of the protein as the C atom on Residue 265 - the last residue in the protein



############################# RUNNING ANALYSES ########################################################################

                                        #######RMSD ANALYSIS#######

R = MDAnalysis.analysis.rms.RMSD(Calpha, Calpha, select=FULL, ref_frame=0)  #this RMSD calculation - this step required a lot of trouble shooting as many of my first attempts saw just the first and last structures' RMSD plotted on a graph with a straight line drawn between the two. A lot of time was spent working out how to get the script to run through the entire trajectory and calculate RMSD of each frame and plot that. Time was spent trying to run the script without ref_frame=0 argument - this here specifies that you are calculating the deviation of the structure through the trajectory compared to the structure in the first frame. Initially I use 'reference' (see above) the place of this argument - but this is the average structure of the trajectory and subsequently not useful as an RMSD reference.

R.run() #command to run the RMSD analysis

df = pd.DataFrame(R.rmsd, columns=['Frame', 'Time(ns)', 'FULL']) #creates a dataframe of the RMSD analysis - trouble shooting needed to make sure the number and titles of columns matched the values put in the columns and the arguments called in plotting section (see below)

                                        #######RMSF ANALYSIS#######

rmsfer = RMSF(Calpha).run() #calculates the RMSF of all the CA atoms in all the frames in the trajectory
aligner = align.AlignTraj(BlaC, reference, select="protein and name CA", in_memory=True) #This compares the average structure across the trajectory (reference - defined above) to each frame in the trajectory. For a long time I did not realise the align.AlignTraj module had to be used for RMSF analysis and could not work out why the script wouldn't work. 

                                        #######RADIUS OF GYRATION ANALYSIS#######

def radgyr(atomgroup, masses, total_mass=None): #This was another section of the code that took a lot of time - there is no predefined module that can calculate the Radius of Gyration, so I had to define a module in the code myself.
    coordinates = atomgroup.positions #This defines the coordinates of the protein as the position of every atom
    center_of_mass = atomgroup.center_of_mass() #the centre of mass is the centre of mass of the protein

    sq_dist = (coordinates-center_of_mass)**2 
    sq = np.sum(sq_dist, axis=1)
    sq_yz = np.sum(sq_dist[:,[1,2]], axis=1) 
    sq_xz = np.sum(sq_dist[:,[0,2]], axis=1) 
    sq_xy = np.sum(sq_dist[:,[0,1]], axis=1)

    sq_rs = np.array([sq, sq_yz, sq_xz, sq_xy])

    RoG_sq = np.sum(masses*sq_rs, axis=1)/total_mass
    return np.sqrt(RoG_sq)

RoG = MDAnalysis.analysis.base.AnalysisFromFunction(radgyr, BlaC.trajectory, protein, protein.masses, total_mass=np.sum(protein.masses)) #Calculates the Radius of Gyration, using the predifined radyr module, the trajectory, the atomgroup 'protein', the masses of the atoms int that atomgroup and then the total mass of the atom group

RoG.run() #cycles through calculating the radius of gyration throughout the trajectory

                                        #######END-TO-END DISTANCE ANALYSIS#######

dist = []
for frame in BlaC.trajectory:
    dist.append(np.linalg.norm(Nterm.positions - Cterm.positions)) #uses the numpy linalg to go through every frame in the trajecotry and calculate the position of the N terminus - position of the C terminus and save it in dist list. Much troubleshooting time was spent making sure my list creation was correct and that values were being added into the list. Further time was spent research how to get the Nterm.positions - Cterm.positions to be added to the list without errors.

dist = np.array(dist) #creates a numpy array from the calculated distance differences

########################################## PLOTTING AND SAVING GRAPHS #################################################

#This section plots the 4 graphs of the 4 analyses the code has completed - this took a lot of trouble shooting when trying to use plotting tools - eventually I settled on Matplotlib. Time was spent also trying to get a dataframe to plot. A lot of time was also spent on trying to get each plot to plot on a separate graph that could be seen, rather than just a single graph with all plots, are multiple overlapping graphs where you could only see the top graph. Further time was spent trying to get it to save into a png file. I realise now it is pretty simple but I have no experience with matplotlib and so took a lot of time to work it out.

fig = df.plot(x='Frame', y=['FULL'], kind='line').get_figure()
fig.savefig('RMSD.png')

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(Calpha.resnums, rmsfer.rmsf)
ax1.set_xlabel("Residue Number")
ax1.set_ylabel(r"RMSF ($\AA$)")
fig1.savefig('RMSF.png')

fig2 = plt.figure()
labels = ['Full Protein']
for col, label in zip(RoG.results.T, labels):
    plt.plot(col, label=label)
plt.legend()
plt.ylabel('Radius of gyration (Å)')
plt.xlabel('Frame');
fig2.savefig('RadGyr.png')

fig3 = plt.figure()
plt.plot( dist, '-k')
plt.xlabel('Frame')
plt.ylabel('end-to-end distance, A')
fig3.savefig('E2E.png')
