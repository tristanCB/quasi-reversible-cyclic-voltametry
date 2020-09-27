#!/usr/bin/env python
import matplotlib.pyplot as plt
import math
import numpy
import time
import argparse

# This function is used for processing simulation data
def find_nearest(array, value):
    array = numpy.asarray(array)
    idx = (numpy.abs(array - value)).argmin()
    return array[idx]

# This is the main function answering question B.5. Bard p. 807
def CV(psi):
	""" 
	Appropriate time windows for different electrochemical techniques:
	Cyclic Voltammetry 
	
	Time parameter: 			RT/F*v
	Usual range of parameter:	v = 0.02-10**6 (V/s)
	Time window					10**-7 - 1 (s)

	Bard p.480

	Potential sweet method:
	more information can be gained in a single experiment by sweeping the potential with time and recording the i-E curve directly. "traversing the three dimentinal i-t-E realm" bard p.226
	Usually potential is varied linearly with time with sweep rates v ranging from 10 mV/s to 1000 V/s with conventinal electrodes and up to 10**6 V/s with UMEs
	In this experiment it is customary to record the current as a function of potential, which is obviously the same as recording current vs. time. 

	wee want to solve numerically 6.2.11 

	The physical interpretation of k0 is straightforward. It simply is a measure of the kinetic facility of a redox couple. Bard p. 96
	
	"""
	######################## PHYSICAL CONSTANTS ########################
	F       = 9.64853*10**4   	# C/mol 	# Faraday's constant
	R       = 8.31447  			# J/mol-K	# Ideal gas constant
	T      	= 298.15 			# K			# Temperature. Default = 298.15
	f       = F/(R*T) 			# V**-1		# Normalized Faraday's constant at room temperature

	######################## Independed Model Params ########################
	alpha 		= 0.5			# NONE					# Transfer Coefficient. The way in which kf and kb depend on potential. 	
	Dm 			= 0.45			# NONE 					# Model diffusion Coefficient 
	l 			= 50			# Number of iterations
	sensitivity = 1				# This is used to 	
	Ei 			= +0.2			# V 		# Initial overpotential
	Ef 			= -0.2			# V  		# Final Overpotential
	DA 			= 1*10**-5		# cm^2/s 	# Diffusion Coefficientof species A 
	k0			= 0.2			# cm/s		# Electrochemical rate constant
	CA			= 1E-06			# mol/cm3	# Bulk concentration of species A
	n 			= 1 			# NONE		# Number of electrons
	A			= 4				# cm^2 		# AREA of the electrode
	# For Ei set to +0.6, iterations should be increased to 500 and sensitivity to 10 inorder from program to produce resonable results
	# Initially following two variables we're used in simulation but errors we're arising
	# They we're removed and equations written in term of overpotential only.
	# Ei0			= 0.015			# E0'		# 
	# Einorm		= Ei0*f 		# NONE		# E0' normalized
	# v 			= 0.05			# V/s 		# sweep rate. Default = 1E-3


	
	######################## Run time parameters ########################
	plot 		= 1	# loop indexing used in plotting
	psiCheck 	= 0	# Will serve to exit a function once a simulation has been run at values of reversal potentials = Ep - 112.5/n mV.
	PrevrevPot 	= 0	# Serves to plot many figures in the same window
	iteration 	= 0	# Logs how many times a simulation was run for a particular psi

	while psiCheck != len(psi):

		# Usually scan rate is set for a specific experiment.
		# This vlaue will be adjusted to run the simulation at specific values of spi which are supplied to the function CV([psi,psi,psi]).
		v = (k0/psi[psiCheck])**2/(math.pi*DA*f)
		
		# Initializing some variable for wave generation
		eta1 = []
		eta2 = []
		eta  = []

		######################## INITIAL CONDTS ########################
		# Electro reactant A is uniformely distributed initially.
		# Potential step is is applied at t = 0 forcing surface conc. of A to zero by converting it faradaically to species B
		# Start by setting up arrays to represent fractional conc. of A and B in each box
		# arrays initiallized to reflect uniform conc. of A and abscence of B
		# Old and new arrays for each species.... 

		######################## MODEL DESCRIPTORS ########################
		psiCHECKCHECK = k0/(math.pi*DA*f*v)**0.5	# Just making sure that psi is really psi
		Lambda = psi[psiCheck]/(math.pi**(-1/2))	# rate of charge transfer to mass transfer
		tk = 2*(Ei-Ef)/v			# Length of the experiment (s)
		delt = tk/l					# Carachteristic time (s)
		delx = (DA*delt/Dm)**0.5	# Characteristic length (cm)


		######################## WAVE GENERATOR ###########################
		# Concatenates two ramps to create sawtooth like potential regime
		# This was inspired from https://petermattia.com/cyclic_voltammetry_simulation/code.html?fbclid=IwAR1c8VkM0oiEy3fHWLuFiwAQn00tfVOv5IXibbHZqeDyERKPJeln2q7Gz0A
		k = range(0,l+1)                						# time index vector
		t = [i*delt for i in k]         						# time vector
		eta1 = [Ei - v*i for i in t]    						# negative scan
		eta2 = [Ei - 2*v*t[int(len(t)/2)] + v*i for i in t]    	# positive scan
		eta.extend(eta1[:int(len(eta1)/2)])						# Combining both scans part 1
		eta.extend(eta2[int(len(eta1)/2):])						# Combining both scans part 2
		Enorm = [i*f for i in eta]    							# normalized overpotential
		####################################################################

		########################### Setting up and boundary conditions ###########################
		# Our initial boundary condition of concentrations of A is uniformely everywhere in solution
		FA 	= [[CA] * l]*(l)# It was chosen to
		# No B in solution is our second boundary condition
		FB 	= [[0] * l]*(l)
		# I assume that I am starting my sweep at a potential where reduction of A is not occuring therefore there will be no flux
		JO  = [0] * (l)	# Flux at electrode surface
		Z 	= [0] * (l)	# Current
		##################################################################################

		####################### Using BV to compute rate constants #######################
		# Rate constants through time / potential regime
		kf = [k0*math.exp(-alpha*n*(Enorm)) for Enorm in Enorm] 	#B.4.11 p799 Bard
		kb = [k0*math.exp((1-alpha)*n*(Enorm)) for Enorm in Enorm]	#B.4.12 p799 Bard

		####################### This begins the numerical simulation part of the code #######################
		for k in range(0, l-1):
			jmax 		= int(math.ceil(4.2*(k**0.5)))	# Bard p. 793
			for j in range(1,jmax-1):
				# Ficks combined first and second law
				FA[k+1][j] = FA[k][j]+Dm*(FA[k][j-1]-2*FA[k][j]+FA[k][j+1])
				FB[k+1][j] = FB[k][j]+Dm*(FB[k][j-1]-2*FB[k][j]+FB[k][j+1]) 
			
			# This flux boundary condition comes from the excel spread sheet provided by:
			# Brown, J. H. (2015). Development and Use of a Cyclic Voltammetry Simulator to Introduce Undergraduate Students to Electrochemical Simulations.
			# Journal of Chemical Education, 92(9), 1490-1496. 
			JO[k+1] = -(kf[k+1]*FA[k+1][1]-kb[k+1]*FB[k+1][1])/(1+(kf[k+1]*delx/DA)+(kb[k+1]*delx/DA))
			# Update surface concentrations also taken from excel spread sheet mentionned above
			FA[k+1][0] = FA[k+1][1] + JO[k+1]*delx/DA
			FB[k+1][0] = FB[k+1][1] - JO[k+1]*delx/DA
			
			# # First Box  from eq. B.4.13 & B.4.14 p.800
			# In appendix B, bard uses these equations. I tried to run a simulation but could not make anything work using this
			# dkf = [tk**0.5*k0/DA**0.5*math.exp(-alpha*n*(k-Einorm)) for k in Enorm] 	#B.4.11 p799 Bard
			# dkb = [tk**0.5*k0/DA**0.5*math.exp((1-alpha)*n*(k-Einorm)) for k in Enorm]	#B.4.12 p799 Bard
			# Z[k+1] = dkf*FA[0][k] - dkb*FA[0][k]
			# FA[0][k+1] = FA[0][k]+Dm*(FA[1][k]-FA[0][k]) - Z[k]*(Dm/l)**(1/2)
			# FB[0][k+1] = FB[0][k]+Dm*(FB[1][k]-FB[0][k]) + Z[k]*(Dm/l)**(1/2)
		
		# Current is computed from flux at electrode
		current = [-n*F*JO*A*0.1 for JO in JO]	# *0.1 to convert from A/m^2 to mA/cm^2

		####################### Post processing of data #######################
		## Extracting peak currents:
		try:
			ipa = min(current)
			ipc = max(current)
		except Exception as e:
			print(f"Error{e}: Try increasing the number of iterations or reconsider your values for psi")
	
		# Calculate peak current as described by Randles-Sevick equation (mA)
		ipeakRS = (2.69*10**5)*n**(3/2)*A*DA**0.5*CA*v**0.5*1000 

		## Extracting Half wave potential
		derivative = numpy.diff(current)
		derivative = derivative.tolist()
		EpHALF = eta[derivative.index(max(derivative))]
		IpHALF = current[derivative.index(max(derivative))]
		tpHALF = t[derivative.index(max(derivative))]

		## Extracting peak potentials from currents
		Epa = eta[current.index(ipa)]	#Peak anodic potential (V)
		Epc = eta[current.index(ipc)]	#Peak cathodic potential (V)

		## Extracting times corresponding to peak currents
		tpa = t[current.index(ipa)]	# (s)
		tpc	= t[current.index(ipc)] # (s)
		
		## Measuring the difference in between peak potentials and converting to mV for comparison to table 6.5.2 p. 243
		Edif = (Epa - Epc)*1000 # Potential difference in (mV)

		## Measuring reversal potential i.e. Elambda
		revPot = Ef*1000	# reversal potential in (mV)
		## Computing EpONEHALF from EpHALF according to equation 6.2.21 p. 231 Bard
		EpONEHALF = EpHALF*1000 - 28/n
		## The following code uses function stated at the begining to extract current corresponding to EpONEHALF
		neareta = find_nearest(eta[:int(len(eta)/2)], EpONEHALF/1000)
		IpONEHALF = current[eta.index(neareta)]
	
		## The difference is computed to be able to rerun simulations at varying potentials to be able to have
		## a simulation which is Elambda = Ep - 112.5/n
		difference = (revPot - Epc*1000) 

		## In the original text bard referencs for table 6.5.2 p.243, Nicholson specifies that all experiments are run at Reversal potential - E1/2.
		## This is the difference we will use to adjust the simulation potentials
		difference2 = revPot - EpONEHALF
		# print("difference2",abs((difference2 + 141/n)))

		## Displays something to the console letting the user know a simulation has been sucsfully completed
		print("##### Iteration:",iteration,". For psi:", psiCHECKCHECK," ///@", "Ei = ",int(Ei*1000)," mV.", " Ef = ", int(Ef*1000)," mV. #####")
		
		# Keeping a variable that has the initial reversal overpotential so every CV generated for psis start consistently
		if iteration == 0:
			Efinitial = Ef
		iteration = iteration+1
		### This compared the difference and adjusts reversal potential to get values close to the 
		### ones which we're used in the derivation for the date of table 6.5.2
		if abs((difference2 + 141/n)) > sensitivity:
			if difference2 < 141/n:
				Ef = Ef+0.001
				# Ei = Ei-0.01
			if difference > 141/n: 
				Ef = Ef-0.001
				# Ei = Ei+0.01

		## Once a simulation comes within resonable range we print out the 
		## simulation variables and output a graph
		else:
			
			print("############################ CASTING PSI:",psiCHECKCHECK," ############################")
			
			# Check for reversibility according to p. 239 Bard
			if (Lambda >= 15) & (k0>=0.3*v**0.5):
				print("Reversible (nernstian)")
			else:
				if (Lambda <= 15) & (Lambda >= 10**(-2*(1+alpha))):
					print("Lambda test check")
					if (k0 <= (0.3*v**0.5)) & (k0 >= 2*10**-5*v**0.5):
						print("Quasireversible")
					else:
						print("NOT Quasireversible")
				else:
					if (Lambda <= 10**(-2*(1+alpha))) & (k0 <= 2*10**-5*v**0.5):
						print("Totally irreversible")
					else:
						print("Out of zone boundaries for LSV")

			## The following just print stuff out during simulation for debugging and analysis ##
			print("----- MODEL AND EXPERIMENTAL PARAMS ----")
			print("Scan rate (V/s):",v)
			print("Lambda is:",Lambda)
			print("tk is (s):",tk)
			print("k0 is (cm/s):",k0)
			print("delx (cm):",delx)
			print("delt (s):",delt)
			print("Dm from other coeffs",DA*delt/delx**2)
			print("----- POTENTIAL AND CURRENT PARAMS ----")
			print("Potential at reversal (mV): ",revPot)
			print("ipa (mA): ",ipa,"ipc (mA): ",ipc)
			print("ipeakRS (mA): ",ipeakRS)
			print("Epa (mV): ",int(Epa*1000),"Epc (mV): ",int(Epc*1000))
			print("Edif (mV)",Edif)
			print("revPot (mV)",revPot)
			print("Ei (mV):",int(Ei*1000),"Ef (mV)",int(Ef*1000))

			###################################################
			################ Plotting results #################
			###################################################
			## Plotting array which are not the same length lead to errors: the pop function 
			##deletes the last item making sure both objects we want to plot have the same length
			while len(eta) > len(current):
				eta.pop()
			while len(eta) < len(current):
				current.pop()
			
			# The following produces the classical graph for CV
			plt.subplot(3, 3, plot)
			plt.plot(eta,current)
			plt.gca().invert_xaxis()
			plotTitle = ("CV with psi of: "+str(psi[psiCheck])+" ---> ")
			plt.title(plotTitle)
			plt.xlabel('Potential (V)')
			plt.ylabel('Current (mA)')
			# This is to specify to plot in another subplot and not in same plot

			# This plots the peak anodic current
			plt.scatter(Epa,ipa,s=80, facecolors='none', edgecolors='r')
			text = "Epa ="+str(int(Epa*1000))+" mV"
			plt.text(Epa, ipa, text, fontsize=8)
			
			# This plots the peak cathodic current
			plt.scatter(Epc,ipc,s=80, facecolors='none', edgecolors='r')
			text = "Epc ="+str(int(Epc*1000))+" mV"
			plt.text(Epc, ipc, text, fontsize=8)

			# This plots the peak half wave current
			plt.scatter(EpHALF,IpHALF,s=80, facecolors='none', edgecolors='b')
			
			plt.scatter(EpONEHALF/1000,IpONEHALF,s=80, facecolors='none', edgecolors='g')
			text = "Epc1/2 ="+ str(int(EpONEHALF))+" mV"
			plt.text(EpONEHALF/1000, IpONEHALF, text, fontsize=8)
			plot = plot + 1

			## Plotting array which are not the same length lead to errors: the pop function 
			##deletes the last item making sure both objects we want to plot have the same length
			while len(t) > len(current):
				t.pop()
			while len(t) < len(current):
				current.pop()
				
			# The first plot is of the same type as 6.5.1 in Bard p. 241
			plt.subplot(3, 3, plot)
			plt.plot(t,current)
			# This plots the peak anodic current
			plt.scatter(tpa,ipa,s=80, facecolors='none', edgecolors='r')
			# This plots the peak cathodic current
			plt.scatter(tpc,ipc,s=80, facecolors='none', edgecolors='r')
			# This plots the peak half wave current
			plt.scatter(tpHALF,IpHALF,s=80, facecolors='none', edgecolors='b')
			# Plot some descriptors on graph
			text = "Epc/2 ="+ str(int(EpHALF*1000))+" mV"
			plt.text(tpHALF, IpHALF, text, fontsize=8)
			
			# Subplot Formatting
			plotTitle = ("Peak Splitting: Epa-Epc = "+str(int(Edif))+" (mV)")
			plt.title(plotTitle)
			plt.xlabel('Time (s)')
			plt.ylabel('Current (mA)')

			## Plotting array which are not the same length lead to errors: the pop function 
			##deletes the last item making sure both objects we want to plot have the same length
			while len(eta) > len(t):
				eta.pop()
			while len(eta) < len(t):
				t.pop()

			# This is to specify to plot in another subplot and not in same plot
			plot = plot + 1

			# This plots potential which the potentiostat imposed to the system
			plt.subplot(3, 3, plot)
			plt.plot(t, eta)

			plotTitle = ("Potential regime. Elambda - Ep1/2 = "+str(int(difference2))+" (mV)")
			plt.title(plotTitle)
			plt.xlabel('Time (s)')
			plt.ylabel('Applied Potential (V)')

			# This is to specify to plot in another subplot and not in same plot
			plot = plot + 1

			# Because we have run a simulation within aggreable conditions we pass 
			#onto running a simulation for the next value of psi specified in function
			###################################################
			########### End of Plotting results ###############
			###################################################
			Ef = Efinitial 				# Reinstating reversal overpotential for the next CV for next psi
			psiCheck = psiCheck + 1		# This changes the value for psi for the next set of simulations
			iteration = 0				# This is just to keep track of things
			print("#####################################################################################")
	
	plt.subplots_adjust(wspace=0.45, hspace=0.6)
	plt.show()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--psi', type=float, nargs=3, help='Three values for psi, the dimensionless intrinsic rate parameter', default=[0.1,1,20])
	args = parser.parse_args()

	# Perform simulations and plots them using matplotlib
	CV(args.psi)


