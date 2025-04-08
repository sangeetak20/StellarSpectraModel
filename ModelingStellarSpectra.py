#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import time
import os
import io

import pandas as pd
import requests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from astropy.table import Table
import numpy as np 
from numpy import load

import astropy.units as u
from astropy.coordinates import SkyCoord

import corner
from astropy.io import fits

import pprint
import glob

from pastamarkers import markers

###needs to be run in terminal###
wget -nv -r -nH --cut-dirs=7
-i /Users/iisan1/Downloads/M15.txt
-B https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/stars/apo25m/M15/

wget -nv -r -nH --cut-dirs=7
-i /Users/iisan1/Downloads/N6791.txt
-B https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/stars/apo25m/N6791/

wget -nv -r -nH --cut-dirs=7
-i /Users/iisan1/Downloads/K2_C4_168-21.txt
-B https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/stars/apo25m/K2_C4_168-21/

wget -nv -r -nH --cut-dirs=7
-i /Users/iisan1/Downloads/060+00.txt
-B https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/stars/apo25m/060+00/
#####

M15 = r'/Users/iisan1/Lab2/apo25m/M15/'
N6791 = r'/Users/iisan1/Lab2/apo25m/N6791/'
K2_C4_168 = r'/Users/iisan1/Lab2/apo25m/K2_C4_168-21/'
field060_00 = r'/Users/iisan1/Lab2/apo25m/060+00/'

field060_00_file = np.loadtxt('060+00.txt', dtype = str)
K2_C4_168_file = np.array(np.loadtxt('K2_C4_168-21.txt', dtype = str))
M15_file = np.loadtxt('M15.txt', dtype = str)
N6791_file = np.loadtxt('N6791.txt', dtype = str)

#note: the length of newset is different than the length of the actual apogee star list
newset = np.isin(data['APOGEE_ID'], apogeestar)
matchedID = np.where(newset == True)
#print(matchedID)
names = data['APOGEE_ID'][matchedID]
logg = data['LOGG'][matchedID]
feh = data['FE_H'][matchedID] 
mgfe = data['MG_FE'][matchedID] 
sife = data['SI_FE'][matchedID]  
teff = data['TEFF'][matchedID]
snr = data['SNR'][matchedID]  

mask = (logg != -9999) & (feh != -9999) & (mgfe != -9999) & (sife != -9999)  \
& (teff != -9999) & (snr > 50)  & (teff < 5700) & (logg < 4) & (feh > -1)

cuts = data[matchedID][mask]

path = r'/Users/iisan1/Lab2/apo25m/'
apogee_dir = os.listdir(path)

M15_matched = []
N6791_matched = []
K2_C4_168_21_matched = []
cluster060_00_matched = []

M15_index = []
N6791_index = []
K2_C4_168_21_index = []
cluster060_00_index = []
for i in range(1, len(apogee_dir)): 
    
    #getting the path for the cluster 
    cluster_path = path + '/' + apogee_dir[i]
    cluster_dir = os.listdir(cluster_path)
    #print(apogee_dir)
    
    #getting the apogee name so that it corresponds to the format
    apogee_name_list = []
    for j in range(len(cluster_dir)): 
        apogee_name = cluster_dir[j].split('-')[2].split('.fits')[0]
        apogee_name_list.append(apogee_name)
    
    #seeing where the apogee and apstar list match up 
    apogee_matched = np.isin(apogee_name_list, names)
    #print(apogee_matched)
    for k in range(len(apogee_matched)): 
        
        #getting the index where the apogee and apstar list are the same and adding it to a list
        if apogee_matched[k] == True: 
            
            if apogee_dir[i] == 'M15':
                M15_matched.append(cluster_dir[k])
                M15_index.append(k)
                
            elif apogee_dir[i] == 'N6791': 
                N6791_matched.append(cluster_dir[k])
                N6791_index.append(k)
                
            elif apogee_dir[i] == 'K2_C4_168-21': 
                K2_C4_168_21_matched.append(cluster_dir[k])
                K2_C4_168_21_index.append(k)
                
            elif apogee_dir[i] == '060+00': 
                cluster060_00_matched.append(cluster_dir[k])
                cluster060_00_index.append(k)
        else: 
            continue 
    
for i in range(1, len(apogee_dir)): 
    #getting the path for the cluster 
    cluster_path = path + '/' + apogee_dir[i]
    cluster_dir = os.listdir(cluster_path)
    
    if apogee_dir[i] == 'M15': 
        matched_index = M15_index
        
    elif apogee_dir[i] == 'N6791': 
        matched_index = N6791_index
        
    elif apogee_dir[i] == 'K2_C4_168-21': 
        matched_index = K2_C4_168_21_index
        
    elif apogee_dir[i] == '060+00': 
        matched_index = cluster060_00_index
        
    #looping through the matched index array to find the star from the cuts we did in Q2
    for j in matched_index: 
        current_star = cluster_dir[j]
        star_path = cluster_path + '/' + current_star
        #print(current_star)
        hdulist = fits.open(star_path)
        
        #defining bitmask and error from fits file
        if len(np.shape(hdulist[1].data)) > 1: 
            spec, error, bitmask = hdulist[1].data[0], hdulist[2].data[0], hdulist[3].data[0] 
            all_fits.append([spec,error,bitmask,wavelength])
        else: 
            spec, error, bitmask = hdulist[1].data, hdulist[2].data, hdulist[3].data
            all_fits.append([spec,error,bitmask,wavelength])
            
        
        #finding where the bit is not equal to zero and setting its corresponding error to a large value
        for k in range(len(bitmask)):
            if k != 0:
                error[k] = 1.e+10
def chop_chips(i):
    spec = all_fits[i][0]
    err = all_fits[i][1]
    wave = all_fits[i][3]

    
    spec_chip = chips(spec)
    err_chip = chips(err)
#     bit_chip = chips(bit)
    wave_chip = chips(wave)
    
#     return spec_chip, err_chip, bit_chip, wave_chip
    return wave_chip, spec_chip, err_chip

def chips(data):
    
    chip1 = data[15150: 15800]
    chip2 = data[15890: 16430]
    chip3 = data[16490: 16950]
    
    return np.array(np.concatenate([chip1, chip2, chip3]))
npz = np.load('cannon_continuum_apogee.npz')
trusted = npz['trusted']
wavelength = npz['wavelengths'] 
trusted_waves = wavelength[trusted == True]
def get_npz_data(num):
    apstar_spec = all_fits[num][0]
    apstar_err = all_fits[num][1]
    apstar_wave = all_fits[num][3]
    
    index = np.searchsorted(apstar_wave,trusted_waves)
    wave_npz = []
    spec_npz = []
    err_npz = []
    for i in range(len(index)):
        if abs(apstar_wave[index[i]+1]-apstar_wave[index[i]]) > abs(apstar_wave[index[i]-1]-apstar_wave[index[i]]):
            data_wave = apstar_wave[index[i]-1]
            data_spec = apstar_spec[index[i]-1]
            data_err = apstar_err[index[i]-1]
            wave_npz.append(data_wave)
            spec_npz.append(data_spec)
            err_npz.append(data_err)
        else:
            data_wave = apstar_wave[index[i]+1]
            data_spec = apstar_spec[index[i]+1]
            data_err = apstar_err[index[i]+1]
            wave_npz.append(data_wave)
            spec_npz.append(data_spec)
            err_npz.append(data_err)
        
    a = {'npz waves':wave_npz, 'npz specs':spec_npz, 'npz errs':err_npz}
    b = Table(a)
    data_table = np.lib.recfunctions.structured_to_unstructured(np.array(Table(b)))
    nans = np.any(np.isnan(data_table), axis=1)
    nanless_table = b[~nans]

    new_wave = np.array(nanless_table['npz waves'])
    new_spec = np.array(nanless_table['npz specs'])
    new_err = np.array(nanless_table['npz errs'])
    

    return new_wave, new_spec, new_err
def poly_and_chunk(wave_arr, spec_arr, err_arr, mininmum, maximum):

    chunk1_wave = []
    chunk1_spec = []
    chunk1_err = []
    
    for i in wave_arr:
        if i >= mininmum and i <= maximum:
            chunk1_wave.append(i)
            if i >= min(chunk1_wave) and i <= max(chunk1_wave):
                chunk1_spec.append(spec_arr[wave_arr == i])
                chunk1_err.append(err_arr[wave_arr == i])
 
    chunk1_wave = np.array(chunk1_wave)
    
    n = []
    m = []
    for i in chunk1_spec:
        for j in i:
            n.append(j)
    chunk1_spec = np.array(n)
    for i in chunk1_err:
        for j in i:
            m.append(j)
    chunk1_err = np.array(m)

    chunk1_coeffs = np.polyfit(chunk1_wave, chunk1_spec, 2, w = 1 / (chunk1_err**2))
    return chunk1_coeffs, chunk1_wave, chunk1_spec, chunk1_err
def plot_polyfit(wave_chunk, coeffs, poly_color): #ax^2 + bx +c
    
    return plt.plot(wave_chunk, (coeffs[0] * (wave_chunk**2)) + (coeffs[1] * (wave_chunk)) + (coeffs[2]), color=poly_color)

def polyfit(wave,spec,err,full_wave):
    coeffs = np.polyfit(wave, spec, 2, w = 1 / (err**2))
    return (coeffs[0] * (full_wave**2)) + (coeffs[1] * (full_wave)) + (coeffs[2])
def chunkifier(i):
    new_wave, new_spec, new_err = get_npz_data(i)
    
#     wave_chip, spec_chip, err_chip = chop_chips(i)
#     full_wave = all_fits[i][3]
    
    c1_wave = new_wave[:241]
    c2_wave = new_wave[241:349]
    c3_wave = new_wave[349:]
    
    c1_spec = new_spec[:241]
    c2_spec = new_spec[241:349]
    c3_spec = new_spec[349:]
    
    c1_err = new_err[:241]
    c2_err = new_err[241:349]
    c3_err = new_err[349:]
    
    c1,w1,s1,e1 = poly_and_chunk(c1_wave, c1_spec, c1_err, 15150, 15800)
    c2,w2,s2,e2 = poly_and_chunk(c2_wave, c2_spec, c2_err, 15890, 16430)
    c3,w3,s3,e3 = poly_and_chunk(c3_wave,c3_spec, c3_err, 16490, 16950)
    
    plt.figure(figsize=(10,7))
    plot_polyfit(w1,c1, 'blue')
    plot_polyfit(w2,c2, 'green')
    plot_polyfit(w3,c3, 'red')
      
    plt.plot(w1,s1,marker='.',markersize=1,ls='none',color='k',label='blue chip data')
    plt.plot(w2,s2,marker='.',markersize=1,ls='none',color='k',label='green chip data')
    plt.plot(w3,s3,marker='.',markersize=1,ls='none',color='k',label='red chip data')
        
    plt.title(f'Spectrum of {names[i]} with polyfitted chunks of spectra chopped by chip data.', fontsize= 14)
    plt.xlabel('Wavelength (Angstrom)', fontsize= 14)
    plt.ylabel('Flux ($ergs/cm^2$)', fontsize= 14)
    plt.legend()
def plot_normalize_spectra(i,plot=False):
    new_wave, new_spec, new_err = get_npz_data(i)
    
    #indices found based on wavelnegths from paper
    c1_wave = new_wave[:241]
    c2_wave = new_wave[241:349]
    c3_wave = new_wave[349:]
    
    c1_spec = new_spec[:241]
    c2_spec = new_spec[241:349]
    c3_spec = new_spec[349:]
    
    c1_err = new_err[:241]
    c2_err = new_err[241:349]
    c3_err = new_err[349:]
    
    wave_chip, spec_chip, err_chip = chop_chips(i)
    full_wave = all_fits[i][3]

    c1= polyfit(c1_wave, c1_spec, c1_err,full_wave[236:3277])
    c2= polyfit(c2_wave, c2_spec, c2_err,full_wave[3688:6107])
    c3= polyfit(c3_wave, c3_spec, c3_err,full_wave[6371:8362])

    line = np.concatenate([c1,c2,c3])
    normal = np.divide(spec_chip[:6303],line)
    error = np.divide(err_chip[:6303], line)

    if plot == True:
        plt.figure(figsize=(10,5))
        plt.plot(wave_chip,normal,marker='.',markersize=1,ls='none',color='k')
        plt.title(f'Normalized Spectrum of {names[i]}', fontsize= 14)
        plt.xlabel('Wavelength (Angstrom)', fontsize= 14)
        plt.ylabel('Flux ($ergs/cm^2$)', fontsize= 14)
        
    else:
        return normal,wave_chip,error
norm_train = []
clean_train = []
err_train = []
for i in randomlist:
    norm_train.append(norms_all[i])
    err_train.append(norms_err[i])
    clean_train.append(clean_all_table[i])
    
norm_train = np.array(norm_train)
clean_train = np.array(clean_train)

len(norm_train),len(clean_train)
cval_set = []
cval_clean = []
cval_err = []
cval_index = []
for i in range(len(norms_all)):
    if i not in randomlist:
        cval_set.append(norms_all[i])
        cval_err.append(norms_err[i])
        cval_clean.append(clean_all_table[i])
        cval_index.append(i)
cval_set = np.array(cval_set)
cval_err = np.array(cval_err)

len(cval_set),len(cval_err)
norm_train = norm_train.tolist()
df_train = pd.DataFrame(np.array(clean_train))
df_train
df_train.insert(14, "NORM_SPEC", norm_train, True) #run once!!
df_train

teff = np.array(df_train['TEFF'])
logg = np.array(df_train['LOGG'])
feH = np.array(df_train['FE_H'])
mgFE = np.array(df_train['MG_FE'])
siFE = np.array(df_train['SI_FE'])

X = np.vstack([np.ones(len(teff)), teff,logg,feH,mgFE,siFE,
              teff**2, teff*logg, teff*feH, teff*mgFE, teff*siFE,
              logg**2, logg*feH, logg*mgFE, logg*siFE,
              feH**2, feH*mgFE, feH*siFE,
              mgFE**2, mgFE*siFE, siFE**2])

def predict_normalized_spectra(teff, logg, feH, mgFE, siFE):
    X = np.vstack([np.ones(len(teff)), teff,logg,feH,mgFE,siFE,
              teff**2, teff*logg, teff*feH, teff*mgFE, teff*siFE,
              logg**2, logg*feH, logg*mgFE, logg*siFE,
              feH**2, feH*mgFE, feH*siFE,
              mgFE**2, mgFE*siFE, siFE**2])
    X = X.T
    f = X @ ptrn
    return f

theta = np.linalg.lstsq(X,norm_train) #is the model that train on for data
star_spectra_normal = np.array(star['NORM_SPEC'])
star_teff = np.array(star['TEFF'])
star_logg = np.array(star['LOGG'])
star_feH = np.array(star['FE_H'])
star_mgFE = np.array(star['MG_FE'])
star_siFE = np.array(star['SI_FE'])
star_spectra_predicted = predict_normalized_spectra(star_teff,star_logg,star_feH,star_mgFE,star_siFE)
def chips2(data):
    
    chip1 = data[236:3277]
    chip2 = data[3688:6107]
    chip3 = data[6371:8362]
    
    return np.array(np.concatenate([chip1, chip2, chip3]))
mystery_spec = fits.open('mystery_spec_wiped.fits')
spec_mys = mystery_spec[1].data
err_mys = mystery_spec[2].data
mask_mys = mystery_spec[3].data

CRVAL = f1[1].header['CRVAL1']
CRDELT1 = f1[1].header['CDELT1'] 
wave_mys = 10**np.arange(CRVAL, CRVAL + len(spec_mys)*CRDELT1, CRDELT1)

index = np.searchsorted(wave_mys,trusted_waves)
wave_npz = []
spec_npz = []
err_npz = []
for i in range(len(index)):
    if abs(wave_mys[index[i]+1]-wave_mys[index[i]]) > abs(wave_mys[index[i]-1]-wave_mys[index[i]]):
        data_wave = wave_mys[index[i]-1]
        data_spec = spec_mys[index[i]-1]
        data_err = err_mys[index[i]-1]
        wave_npz.append(data_wave)
        spec_npz.append(data_spec)
        err_npz.append(data_err)
    else:
        data_wave = wave_mys[index[i]+1]
        data_spec = spec_mys[index[i]+1]
        data_err = err_mys[index[i]+1]
        wave_npz.append(data_wave)
        spec_npz.append(data_spec)
        err_npz.append(data_err)

a = {'npz waves':wave_npz, 'npz specs':spec_npz, 'npz errs':err_npz}
b = Table(a)
data_table = np.lib.recfunctions.structured_to_unstructured(np.array(Table(b)))
nans = np.any(np.isnan(data_table), axis=1)
nanless_table = b[~nans]

new_wave = np.array(nanless_table['npz waves'])
new_spec = np.array(nanless_table['npz specs'])
new_err = np.array(nanless_table['npz errs'])

c1_wave = new_wave[:241]
c2_wave = new_wave[241:349]
c3_wave = new_wave[349:]

c1_spec = new_spec[:241]
c2_spec = new_spec[241:349]
c3_spec = new_spec[349:]

c1_err = new_err[:241]
c2_err = new_err[241:349]
c3_err = new_err[349:]

wave_chip = chips2(wave_mys)
spec_chip = chips2(spec_mys)
err_chip = chips2(err_mys)

c1= polyfit(c1_wave, c1_spec, c1_err,wave_mys[236:3277])
c2= polyfit(c2_wave, c2_spec, c2_err,wave_mys[3688:6107])
c3= polyfit(c3_wave, c3_spec, c3_err,wave_mys[6371:8362])

line = np.concatenate([c1,c2,c3])
normal = np.divide(spec_chip,line)
error_norm = np.divide(err_chip, line)


