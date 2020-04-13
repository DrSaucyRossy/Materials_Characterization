# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:22:25 2020

@author: Ross Anthony

Calculating the contrast for graphene
"""
from cmath import exp, pi

import numpy as np
import pandas as pd


#refractive indices
def Intensity(n1,wavelength, dielectric_thickness):
    n0 = 1 # air
    n2 = 1.77 #Al203
    n3 = 5.6-0.4j #Silicon

    d1 = 0.34 #single layer graphene thickness in nm
    d2 = dielectric_thickness
    # phase shift due to changes in the optical path
    phi1= 2*pi*n1*d1/wavelength
    phi2= 2*pi*n2*d2/wavelength
    
    e1 = exp(1j*(phi1+phi2))
    e2 = exp(1j*(phi1-phi2))
    e3 = exp(-1j*(phi1+phi2))
    e4 = exp(-1j*(phi1-phi2))

    #relative indices of refraction
    r1 = (n0-n1)/(n0+n1)
    r2 = (n1-n2)/(n1+n2)
    r3 = (n2-n3)/(n2+n3)
    
    I = (abs((r1*e1+r2*e2+r3*e3+r1*r2*r3*e4)/(e1+r1*r2*e2+r1*r3*e3+r2*r3*e4)))**2
    return(I)
    

def Contrast(n1,wavelength,dielectric_thickness):
    Contrast = (Intensity(1,wavelength,dielectric_thickness)-Intensity(n1, wavelength,dielectric_thickness))/Intensity(1,wavelength,dielectric_thickness)
    return Contrast

 #grapite refractive index varies by 5% from 300 to 590
n1graphite = 2.6 - 1.3j #Graphite
n1air = 1


''' Wavelength range'''
wavelength_data_points = 110
wavelength_start = 300
wavelength_finish = 630
wavelength_step = (wavelength_finish-wavelength_start)/wavelength_data_points
wavelength = []

for i in range(wavelength_data_points):
    wavelength.append(wavelength_start+i*wavelength_step)
wavelength.append(wavelength_finish)

'''
This could also be done using np.arrange, see below
wavelengths = np.arange(wavelength_start,wavelength_finish+wavelength_step,wavelength_step) 
'''

'''Dielectric_Thickness Range'''

dielectric_thickness_data_points= 100
dielectric_thickness_start= 50
dielectric_thickness_finish= 550
dielectric_step = (dielectric_thickness_finish-dielectric_thickness_start)/(dielectric_thickness_data_points)
dielectric_thickness = []

for i in range(dielectric_thickness_data_points):
    dielectric_thickness.append(dielectric_thickness_start+i*dielectric_step)
dielectric_thickness.append(dielectric_thickness_finish)

t = Intensity(n1air,400,550)
u = Intensity(n1graphite,400,550)

Final = np.array([])
for i in range(len(dielectric_thickness)):
    run = np.array([])
    for j in range(len(wavelength)):
        C = Contrast(n1graphite,wavelength[j],dielectric_thickness[i])
        run = np.append(run,C)
    Final = np.append(Final, run)