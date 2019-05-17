# module binary


import sys
sys.path.append('..')
from astropy import units as u
from astropy import constants as const
import astropy
import numpy
import numpy as np
import math
import timeit
import zams
import argparse
import System
from scipy.integrate import trapz


def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Nsys', type = int, required = True)
    parser.add_argument('--binaryfile', type = str, required = True)
    parser.add_argument('--Mexp', type = float, default = 2.35)

    args = parser.parse_args()

    return args

def sampleIMF(args, Mlow,nbin):

    if args.Mexp == 0.0:
        return np.random.uniform(Mlow,100,nbin)
    alpha = args.Mexp
    def IMF(m):
        return m**-alpha
    mm = np.linspace(Mlow,100,1000)
    A=1./trapz(IMF(mm),x=mm)
    def invIMF(I):
        return (I*(1-alpha)/A + Mlow**(1-alpha))**(1./(1-alpha))
    xx=np.random.uniform(0,1,nbin)
    return invIMF(xx)




def roche_lobe(m1, m2):
	if not isinstance(m1, u.Quantity):
		raise ValueError("m1 must be a Quantity with mass units")
	if not isinstance(m2, u.Quantity):
		raise ValueError("m2 must be a Quantity with mass units")

	q_mass1 = m1 / m2
	q_mass3 = q_mass1 ** (1.0 / 3.0)
	q_mass2 = q_mass1 ** (2.0 / 3.0)
	lobe = (0.49 * q_mass2) / (0.6 * q_mass2 + numpy.log(1.0 + q_mass3))
	return lobe  # Dimensionless. In units of orbital separation


def separation_to_period(separation, m1, m2):
	if not isinstance(m1, u.Quantity):
		raise ValueError("m1 must be a Quantity with mass units")
	if not isinstance(m2, u.Quantity):
		raise ValueError("m2 must be a Quantity with mass units")
	if not isinstance(separation, u.Quantity):
		raise ValueError("sep must be a Quantity with length units")

	mbin = m1 + m2
	period = numpy.sqrt(4.0 * math.pi ** 2.0 * separation ** 3.0 / (const.G * mbin))
	return period.to('day')






def random_separation(amin, amax, nsample=1):
	"""Returns a random number for the radius of separation
	of a binary. min cannot be zero
	"""
	if not isinstance(amin, u.Quantity):
		raise ValueError("amin must be a Quantity with mass units")
	if not isinstance(amax, u.Quantity):
		raise ValueError("amax must be a Quantity with mass units")

	amin = amin.to('Rsun')
	amax = amax.to('Rsun')
	if numpy.min(amin.value) <= 0:
		raise ValueError("min cannot be 0 or negative")
	if amax.value.any() < amin.value.any():
		raise ValueError("min must be less than or equal to max")

	y =  numpy.random.random(nsample)*(numpy.log(amax.value) - numpy.log(amin.value))
	return numpy.exp(y)*amin


def random_eccentricity(nsample=1):
	# Generate a random number from the thermal eccentricity distribution f(e)=2e (CDF: F(e)=e^2)
	e = numpy.random.random(nsample)
	return numpy.sqrt(e)








##def f(x, a1 = 1.3, a2 = 2.3, a3 = 2.35, m_min=0.08, m_1 = 0.5, m_2 = 1.0, m_max=150.0):
##	fo_inv = 1./(1.-a1) * (m_1**(1.-a1)-m_min**(1.-a1))/(m_min**(-a1)) \
##		+ 1./(1.-a2) * (m_1/m_min)**(-a1) * (m_2**(1.-a2)-m_1**(1.-a2))/(m_1**(-a2)) \
##		+ 1./(1.-a3) * (m_1/m_min)**(-a1) * (m_2/m_1)**(-a2) * (m_max**(1.-a3)-m_2**(1.-a3))/(m_2**(-a3))
##
##	fo = 1. / fo_inv
##
##	if type(x) is numpy.ndarray:
##		value = numpy.zeros(len(x))
##		idx1 = numpy.where((x >= m_min) * (x < m_1))
##		if len(idx1[0]) > 0:
##			value[idx1] = fo * numpy.power(x[idx1]/m_min,-a1)
##		idx2 = numpy.where((x >= m_1)  * (x < m_2))
##		if len(idx2[0]) > 0:
##			value[idx2] = fo * numpy.power(m_1/m_min,-a1) * numpy.power(x[idx2]/m_1,-a2)
##		idx3 = numpy.where((x >= m_2) * (x < m_max))
##		if len(idx3[0]) > 0:
##			value[idx3] = fo * numpy.power(m_1/m_min,-a1) * numpy.power(m_2/m_1,-a2) * numpy.power(x[idx3]/m_2,-a3)
##		idx4 = numpy.where(x > m_max)
##		if len(idx4[0]) > 0:
##			value[idx4] = 0.0
##		idx5 = numpy.where(x < m_min)
##		if len(idx5[0]) > 0:
##			value[idx5] = 0.0
##
##	elif type(x) is float or type(x) is numpy.float64:
##		value=0.0
##		if (x >= m_min) and (x < m_1):
##			value = fo * numpy.power(x/m_min,-a1)
##		elif (x >= m_1)  and (x < m_2):
##			value = fo * numpy.power(m_1/m_min,-a1) * numpy.power(x/m_1,-a2)
##		elif (x >= m_2) and (x < m_max):
##			value = fo * numpy.power(m_1/m_min,-a1) * numpy.power(m_2/m_1,-a2) * numpy.power(x/m_2,-a3)
##		else:
##			value=0.0
##
##	return value
####
##
##
##def g(x, a1 = 1.3, m_min=0.08, m_max=150.0):
##	fo_inv = 1./(1.-a1) * (m_max**(1.-a1)-m_min**(1.-a1))/(m_min**(-a1)) 
## 
##	fo = 1. / fo_inv
##
##	if type(x) is numpy.ndarray:
##		value = numpy.zeros(len(x))
##		idx1 = numpy.where((x >= m_min) * (x < m_max))
##		if len(idx1[0]) > 0:
##			value[idx1] = fo * numpy.power(x[idx1]/m_min,-a1)
##	elif type(x) is float or type(x) is numpy.float64:
##		value=0.0
##		if (x >= m_min) and (x < m_max):
##			value = fo * numpy.power(x/m_min,-a1)
##
##	return value
##
##
##




def generate_bse_binaries_array(args):
	# Input
	# fileout :: str, name of output file.
	# nbin    :: Number of binary pairs to make.
	fileout = args.binaryfile
	nbin = args.Nsys
	time_to_evolve 	= 100. * u.Myr  # float, Myr. Max evolution time
	z 				= 0.02 # Metallicity
	zams.prepare_coefficients(z)
	sep_max 		= 1.e5  # float.
	sep_max 		= numpy.full((nbin), sep_max*u.Rsun)
	M2_min 			= 0.08

	# Define master array
	data = numpy.zeros([nbin,6])  # numpy.array, floats.

	#Draw M1 and M2
	Mlow 		= 2.0
	Mup 		= 100.0
	Qmin 		= 0.0
	Qmax 		= 1.0
	x 			= numpy.random.rand(nbin)
	data[:,0] 	= sampleIMF(args, Mlow, nbin)#inv_imf_2part(x, Mlow, Mup)		    # Draw M1 stars
	Q 			= (Qmax-Qmin)*numpy.random.rand(nbin)+Qmin
	data[:,1] 	= data[:,0]*Q	# Draw M2 stars from M1 and q-ratio.

	#Initiate re sampling for those M2 < M2_min
	dum 		= numpy.zeros([numpy.size(numpy.where(data[:,1] < M2_min)),3])
	count_limit1 = 1000
	if numpy.size(dum[:,0]) != 0:
		print (numpy.where(data[:,1] < M2_min))
		dum_index = numpy.where(data[:,1] < M2_min)
		dum[:,2] 	= dum_index[0] 			# Array indexes.
		dum[:,0] 	= data[dum[:,2].astype(int),0] 	# M1
		dum[:,1] 	= data[dum[:,2].astype(int),1] 	# M2
		count = 0
		
		while numpy.size(dum) != 0 and count < count_limit1:
			count = count + 1
			#re-sample M2 from M1.
			dum[:,1] 		= dum[:,0]*numpy.random.rand(numpy.size(dum[:,2]))
			index 			= numpy.where((dum[:,1] > M2_min))
			data[dum[index[0],2].astype(int),1] 	= dum[index[0],1]
			dum = dum[numpy.where(dum[:,1] < M2_min),:]
			#re do exercise
			if numpy.size(numpy.where(data[:,1] < M2_min)) != 0:
				where_fail 	= numpy.where(data[:,1] < M2_min)
				dum 		= numpy.zeros([numpy.size(where_fail),3])
				dum[:,2] 	= where_fail[0]
				dum[:,0] 	= data[dum[:,2].astype(int),0] # M1
				dum[:,1] 	= data[dum[:,2].astype(int),1] # M2
		if count > count_limit1:
			print ('while loop 1 exceeded count limit of ', count_limit1)
			
	#Draw eccentricity
	data[:,3] 	= random_eccentricity(nbin)
	ecc 		= data[:,3]
	m1 			= data[:,0]*u.Msun
	m2 			= data[:,1]*u.Msun
	rl1 		= roche_lobe(m1, m2)
	rzams1 		= zams.rzams(m1)
	sep_min 	= rzams1/(rl1*(1.-ecc))
	
	#make check that no sep_min > sep_max
	print ('1st: sep_min > sep_max', numpy.size(numpy.where(sep_min.value > sep_max)))

##
	#1 find index in violation
	index1 = numpy.where(sep_min.value > sep_max)
	count_limit2 = 1000
	if numpy.size(index1[0]) != 0:
				
		#index2 is the one we loop over henceforth
		index2 = index1
		
		#2 select values for resampling
		#make dummy array in which 
		dum 		= numpy.zeros((numpy.size(index2[0]),3))
		dum[:,0] 	= ecc[index2[0]] 							#Ecc
		dum[:,1] 	= sep_min[index2[0]] 						#sep_min
		dum[:,2] 	= index2[0] 								#index
		count = 0
		while numpy.size(index2) != 0 and count < count_limit2:
			count = 1 + count
			dum[:,0] 	= random_eccentricity(numpy.size(index2)) 			#ecc
			dum[:,1] 	= rzams1[index2[0]]/(rl1[index2[0]]*(1.-dum[:,0]))  #sep_min
			re_insert  	= numpy.where(dum[:,1] < sep_max[0])
			data[index2[0][re_insert],3] 	= dum[re_insert,0]

			#update data set
			ecc2 							= data[index2[0][re_insert],3]
			m12 							= data[index2[0][re_insert],0]*u.Msun
			m22								= data[index2[0][re_insert],1]*u.Msun
			rl12	 						= roche_lobe(m12, m22)
			rzams12 						= zams.rzams(m12)
			sep_min[index2[0][re_insert]] 	= rzams12/(rl12*(1.-ecc2))
			#re do exercise
			index2 = numpy.where(sep_min.value > sep_max[0])
			if numpy.size(index2) != 0:
				#re define dum array to new number of index.
				dum 		= numpy.zeros((numpy.size(index2[0]),3))
				dum[:,0] 	= ecc[index2[0]] 							#Ecc
				dum[:,1] 	= sep_min[index2[0]] 						#sep_min
				dum[:,2] 	= index2[0] 								#index
			
			
		if count > count_limit2:
			print ('while loop 2 exceeded count limit of ', count_limit2)
	print ('2nd sep_min > sep_max after while loop', numpy.size(numpy.where(sep_min.value > sep_max[0])))
		
	# Update data set.
	ecc 		= data[:,3]
	m1 			= data[:,0]*u.Msun
	m2 			= data[:,1]*u.Msun
	rl1 		= roche_lobe(m1, m2)
	rzams1 		= zams.rzams(m1)
	sep_min 	= rzams1/(rl1*(1.-ecc))
	
	sep 		= random_separation(sep_min, sep_max * u.Rsun, nbin)
	tb 			= separation_to_period(sep, m1, m2)

	data[:,2] 	= tb
	data[:,4] 	= z
	data[:,5] 	= time_to_evolve


	#write data to file
	numpy.savetxt(fileout,data,fmt='%12.10f %12.10f %12.10f %12.10f %4.3f %4.1f', header=str(int(nbin)), comments='' )


if __name__ == "__main__":
    args=parse_commandline()

    generate_bse_binaries_array(args)

