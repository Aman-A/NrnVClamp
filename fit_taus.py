#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 12:48:40 2018

@author: amanaberra
"""
import numpy as np
from scipy.optimize import curve_fit

## regression function
def double_exp(x, a, b, c, d,e,f):
    return a * np.exp(-b * (x - e))  + c * np.exp(-d * (x - f))
def single_exp(x, a, b, c):
    f = a * (1 - np.exp(-(x - c)/b))
    f[x<c] = 0
    return f

def fit_double_exp(t_vec,y_vec):
    fits,covs = curve_fit(double_exp,t_vec,y_vec,p0=(-20.0,2,-1.0,5,0,0))
    tau1 = 1/fits[1]
    tau2 = 1/fits[3]    
    return tau1, tau2

def fit_single_exp(t_vec,y_vec,p0=None,bounds=None):        
    if y_vec[np.abs(y_vec[0]).argmax()] < 0: # inward current, prob Na
        if p0 is None:
            p0 = (y_vec.min(),0.1,t_vec[1])
        if bounds is None:
            bounds = ([1.5*y_vec.min(),0.005,0],[0,0.8,0.1])        
    else: # outward current, prob K
        if p0 is None:
            p0 = (y_vec.max(),1,t_vec[1])
        if bounds is None:
            bounds = ([0,0.005,0],[1.5*y_vec.max(),5,0.1])
    fits,covs = curve_fit(single_exp,t_vec,y_vec,p0=p0,bounds=bounds,method='dogbox')
    tau1 = fits[1]       
    return tau1,fits
    
    