#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 15:48:34 2016

@author: undarmaabayarsaikhan
"""

# Initial setup
# Basic calculations
import math
import numpy
import random
import matplotlib.pyplot as plt

#x(t)~N[sigma^2/(2k^2)*(1-exp(-kt))^2),sigma^2/(2k)*(1-exp(-2kt))]
# reference http://www-ljk.imag.fr/AMMSI/Manifestations/Toulouse2014/Grall%20AMMSI%2020.01.14.pdf
sigmar=0.02
sigmaB=0.005
kB=0.1
sigmaC=0.01
kC=0.1
lambdaB_0T = 0.01
lambdaC_0T = 0.03
capital_T = 120
M = 50000

def calc_xt_mean(sigma, k, t): #calculates mean of one variable x(t)
#    print('calc_xt_mean')
    return sigma**2/(2*k**2)*(1-math.exp(-k*t))**2

def calc_xt_var (sigma, k, t): #calculates variance of one variable x(t)
#    print('x_var')
    return sigma**2/(2*k)*(1-math.exp(-2*k*t))
    
def gene_xt(sigma, k, t): #generate x(t)
    xt_mean = calc_xt_mean(sigma,k,t)
    xt_std = math.sqrt(calc_xt_var(sigma,k,t))
 #   print('gene_x')
    return random.gauss(xt_mean,xt_std)

def calc_y(sigma, k, t): #calculate y(t)
    return sigma**2/(2*k)*(1-math.exp(-2*k*t))
    
def calc_Gt(k, t, capital_t): #calculate G(t)
    return math.exp(-k*(capital_t - t))/k 

def cal_lambda(sigma, k, t, capital_T, lambda_0T): #calculate lambda
    xt = gene_xt(sigma,k,t)
    yt = calc_y(sigma,k,t)
    Gt = calc_Gt(k,t,capital_T)
    return lambda_0T + math.exp(-k*(capital_T - t)) * (xt + yt * Gt) 
    
    
#lambda_matrix[i] is path of lambda 0 to 120
lambda_matrix = []
for i in numpy.arange(0, M):
    lambda_vector =[]
    for t in numpy.arange(0, capital_T+1):
        lambda_vector.append(cal_lambda(sigmaB, kB, t, capital_T, lambdaB_0T))
    lambda_matrix.append(lambda_vector)
        
