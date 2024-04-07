#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
from astropy import constants as const
import numpy as np
from astropy import units as u
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
import os
import subprocess

def ReadParams(input_file):
    """
    Function for reading parameters from the input json-file.

    Parameters
    ----------
    input_file : string
        Name of the input json-file.

    Returns
    -------
    ds : float
        Source distant in kpc.
    s_flux : float
        Baseline flux.
    s_vel : float
        Source relative velocity in km/s.
    l_mass : float
        Lens mass in Solar masses.
    dl : float
        Lens distant in kpc.
    u0 : float
        Impact parameter.
    t0 : float
        Time of magnification peak.

    """
    with open(input_file, 'r') as f:
        s = json.load(f)
    
    ds = s["point source"]["Ds"]
    s_flux = s["point source"]["flux"]
    s_vel = s["point source"]["velocity"]
    
    l_mass = s["point lens"]["mass"]
    dl = s["point lens"]["Dl"]
    u0 = s["point lens"]["impact parameter"]
    t0 = s["point lens"]["t0"]
    
    return ds, s_flux, s_vel, l_mass, dl, u0, t0

def CreateDir(save_dir):
    """
    Function for creating directory for saving the results.

    Parameters
    ----------
    save_dir : string
        Directory for saving the results.

    Returns
    -------
    None.

    """
    if not(os.path.exists(save_dir)):
        key = True
        while key:
            ans = input("No path %s. Do you want to create path? (y/n) \n" %save_dir)
            if (ans == "y"):
                subprocess.run("mkdir -p %s" %save_dir, shell=True)  
                key = not(key)
            elif(ans == "n"):
                print("Exit")
                key = not(key)
                return
    return

def EinsteinRadius(l_mass, ds, dl):
    """
    Function for calculating the Einstein radius.

    Parameters
    ----------
    l_mass : float
        Lens mass in Solar masses.
    ds : float
        Source distant in kpc.
    dl : float
        Lens distant in kpc.

    Returns
    -------
    theta_E : float
        Einstein radius in arcsec.

    """
    mass = l_mass * const.M_sun
    G = const.G
    c = const.c
    aconv = np.rad2deg(1.0) * 3600.0 * u.arcsecond
    theta_E = np.sqrt(4.0 * (G * mass/c/c).to('kpc') * (ds - dl)/dl/ds/u.kpc) * aconv
    return(theta_E)
    
def EinsteinCrossTime(theta_E, dl, s_vel):
    """
    Function for calculating the Einstein radius crossing time.

    Parameters
    ----------
    theta_E : float
        Einstein radius in arcsec.
    dl : float
        Lens distant in kpc.
    s_vel : float
        Source relative velocity in km/s.

    Returns
    -------
    float
        Einstein radius crossing time.

    """
    return(((theta_E.to('radian').value * dl * u.kpc).to('km')/s_vel/u.km * u.s).to('day'))

def Y(t, t0, tE, u0):
    """
    Function calculating the coordinates of the unlensed source at time t

    Parameters
    ----------
    t : float
        Time.
    t0 : float
        Time of magnification peak.
    tE : float
        Einstein radius crossing time.
    u0 : float
        Impact parameter.

    Returns
    -------
    tuple
        Coordinates of the unlensed source.

    """
    y1 = (t - t0)/tE.value
    y2 = np.ones(len(t)) * u0
    return(y1, y2)

def X(y1, y2, plus=True):
    """
    Function calculating the coordinates of the x+ and x- image at time t

    Parameters
    ----------
    y1, y2 : float
        Coordinates of the unlensed source.
    plus : bool
        Specifies for which image, x+ or x-, the coordinates are calculated. The default is True.

    Returns
    -------
    tuple
        Coordinates of the x+ and x- image.

    """
    Q = np.sqrt(y1**2 + y2**2 + 4) / np.sqrt(y1**2 + y2**2)
    if plus:
        x1 = 0.5 * (1 + Q) * y1
        x2 = 0.5 * (1 + Q) * y2
    else:
        x1 = 0.5 * (1 - Q) * y1
        x2 = 0.5 * (1 - Q) * y2
    return(x1, x2)
       
def Magnification(y1, y2, s_flux):
    """
    Function calculating the magnification as a function of time

    Parameters
    ----------
    y1, y2 : float
        Coordinates of the unlensed source.
    s_flux : float
        Baseline flux.

    Returns
    -------
    float
        Magnification.

    """
    y = np.sqrt(y1**2 + y2**2)
    return(s_flux * (y**2 + 2) / y / np.sqrt(y**2 + 4))
   
def X_ext_source(y1, y2, r, plus=True):
    """
    Function calculating the coordinates of the x+ and x- image at time t, when the source has an extended size.
    We assign to the source a circular shape.

    Parameters
    ----------
    y1, y2 : float
        Coordinates of the unlensed source.
    r : float
        Source radius.
    plus : bool
        Specifies for which image, x+ or x-, the amplification is calculated. The default is True.

    Returns
    -------
    tuple
        Coordinates of the x+ and x- image.

    """
    phi = np.linspace(0.0, 2*np.pi, 360)
    dy1 = r * np.cos(phi)
    dy2 = r * np.sin(phi)
    yy1 = y1 + dy1
    yy2 = y2 + dy2
    Q =  np.sqrt(yy1**2 + yy2**2 + 4)/(np.sqrt(yy1**2 + yy2**2))
    if plus:
        x1 = 0.5 * (1 + Q) * yy1
        x2 = 0.5 * (1 + Q) * yy2
    else:
        x1 = 0.5 * (1 - Q) * yy1
        x2 = 0.5 * (1 - Q) * yy2
    return(x1, x2)

def Plot(t, magn, xp1_e, xp2_e, xm1_e, xm2_e, tt, magn_ext, n):
    """
    Function for plotting images at time t.

    Parameters
    ----------
    t : float
        Time.
    magn : float
        Magnification.
    xp1_e, xp2_e : float
        Coordinates of the x+ image, when the source has an extended size.
    xm1_e, xm2_e : float
        Coordinates of the x- image, when the source has an extended size.
    tt : float
        In this moment of time we calculate magnification of the extended size source.
    magn_ext : float
        Magnification of the extended size source.
    n : int
        Number of iteration.

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(1, 2, figsize=(20,10))
    font = {'size'   : 20, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    
    ax[1].set_ylabel(r'$A(t)$')
    ax[1].set_xlabel(r'$(t-t_0)/t_E$')
    
    ax[0].plot([0.0], [0.0], '*', markersize=20, color='gold')
    ax[0].plot(xp1_e, xp2_e, lw=2, color='cornflowerblue')
    ax[0].plot(xm1_e, xm2_e, lw=2, color='cornflowerblue')
    
    ax[1].plot(t, magn, 'k-')
    ax[1].plot(tt, magn_ext, 'o', markersize=10, color='cornflowerblue')
    
    ax[0].set_xlim([-2.5,2.5])
    ax[0].set_ylim([-2.3,2.7])
    circle=plt.Circle((0,0), 1, color='black', fill=False)
    ax[0].add_artist(circle)
    
    plt.savefig('images/' + '{:03d}'.format(n) + '.png')
    plt.clf()
    plt.cla()
    plt.close()
    return   

def Animation():
    """
    Function for creating an animation, using plots from the folder 'images/'

    Returns
    -------
    None.

    """
    if os.path.exists('animation.mp4'):
        subprocess.run('rm animation.mp4', shell=True)
    cmd = "ffmpeg -framerate 20 -pattern_type glob -i 'images/*.png' -c:v libx264 -r 30 -pix_fmt yuv420p animation.mp4"
    subprocess.run(cmd, shell=True)
    return

def run(input_file):
    """
    The main function

    """
    ds, s_flux, s_vel, l_mass, dl, y0, t0 = ReadParams(input_file)
    theta_E = EinsteinRadius(l_mass, ds, dl)
    tE = EinsteinCrossTime(theta_E, dl, s_vel)

    t = t0 + np.linspace(-8, 8, 2000) * tE.value

    y1, y2 = Y(t, t0, tE, y0)
    xp1, xp2 = X(y1, y2, plus=True)
    xm1, xm2 = X(y1, y2, plus=False)
    magn = Magnification(y1, y2, s_flux)

    fig, ax = plt.subplots(1, 2, figsize=(20,10))
    font = {'size'   : 20, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    
    ax[0].plot(y1, y2, '--', label='source traj.')
    ax[0].plot(xp1, xp2, '--', label='image $x_+$')
    ax[0].plot(xm1, xm2, '--', label='image $x_-$')
    ax[1].plot(t, magn, 'k-')
    
    CreateDir('images')
    
    print('Start counting and plotting coordinates and magnification at different moments of time')
    t_sparse = t0 + np.linspace(-2, 2, 100) * tE.value
    color = iter(cm.rainbow(np.linspace(0, 1, t_sparse.size)))
    n = 0
    for tt in tqdm(t_sparse):
        n += 1
        c = next(color)
        y1_ext, y2_ext = Y(np.array([tt]), t0, tE, y0)
        xp1_e, xp2_e = X_ext_source(y1_ext, y2_ext, 0.05, plus=True)
        ax[0].plot(xp1_e, xp2_e, color=c, lw=2)
        xm1_e, xm2_e = X_ext_source(y1_ext, y2_ext, 0.05, plus=False)
        ax[0].plot(xm1_e, xm2_e, color=c, lw=2)
        magn_ext = Magnification(y1_ext, y2_ext, s_flux)
        ax[1].plot([tt], [magn_ext], 'o', markersize=10, color=c)
        Plot(t, magn, xp1_e, xp2_e, xm1_e, xm2_e, tt, magn_ext, n)
    
    print('\n' + 'Start creating animation')
    Animation()    
    
    ax[0].set_xlim([-2,2])
    ax[0].set_ylim([-1.8,2.2])
    ax[0].plot([0.0], [0.0], '*', markersize=20, color='gold')
    circle=plt.Circle((0,0), 1, color='black', fill=False)
    ax[0].add_artist(circle)
    ax[0].legend()
    ax[1].set_ylabel(r'$A(t)$')
    ax[1].set_xlabel(r'$(t-t_0)/t_E$')
    plt.savefig('Microlensing_full.png')
    plt.clf()
    plt.cla()
    plt.close()    

    print('Done!')   
    return

if __name__ == "__main__":
    run("INPUT.json")
