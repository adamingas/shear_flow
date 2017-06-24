from __future__ import division
import numpy as np
import sys
import os
import math
import config
def walk():
    """
    Performs the random walk. Creates the Matrix with the coefficients of the velocity profile created when
    one of the spheres moves through the fluid.
    """
    global hydro
    if str.upper(hydro) in ("YES","Y","HYDRO"):
        hydro = "yes"
        for j in xrange(runs):
            out = open("Run{}".format(j), "w")

            x1, y1, x2 = 0,0,0
            y2 = init_separation
            out.write("{} {} {} {} {}\n".format(0,x1,y1,x2,y2))
            for i in xrange(steps):
                #S is the stokeslet tensor
                distance_sq = (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1)
                Sxy = (3/(4*math.sqrt(distance_sq)))*(
                (x2-x1)*(y2 - y1)/
                             distance_sq)
                Sxx = (3/(4*math.sqrt(distance_sq)))*(1+
                (x2-x1)*(x2-x1)/distance_sq)
                Syy =  (3/(4*math.sqrt(distance_sq)))*(1+
                (y2 - y1)*(y2-y1)/distance_sq)

                Matrix = np.array([[1,0,-Sxx,-Sxy],[0,1,-Sxy,-Syy],[-Sxx,-Sxy,1,0],[-Sxy,-Syy,0,1]])
                Minv = np.linalg.inv(Matrix)
                random_array = np.random.randn(4)
                x1new= x1 + np.dot(Minv,random_array)[0]*np.sqrt(constant * time_step)
                y1new = y1+ np.dot(Minv, random_array)[1] * np.sqrt(
                    constant * time_step)
                x2new = x2 + np.dot(Minv, random_array)[2] * np.sqrt(
                    constant * time_step)
                y2new = y2 + np.dot(Minv, random_array)[3] * np.sqrt(
                    constant * time_step)
                x1,y1,x2,y2 = x1new,y1new,x2new,y2new
                out.write("{} {} {} {} {}\n".format(time_step*(i+1), x1, y1, x2, y2))
            update_progress(j / (runs))
            out.close()
    elif str.upper(hydro) in ("NO","N","NOHYDRO") :
        hydro = "no"
        for j in xrange(runs):
            out = open("Run{}".format(j), "w")

            x1, y1, x2 = 0,0,0
            y2 = init_separation
            out.write("{} {} {} {} {}\n".format(0,x1,y1,x2,y2))
            for i in xrange(steps):
                #S is the stokeslet tensor

                random_array = np.random.randn(4)
                x1new= x1 + random_array[0]*np.sqrt(constant * time_step)
                y1new = y1+ random_array[1] * np.sqrt(constant * time_step)
                x2new = x2 + random_array[2] * np.sqrt(constant * time_step)
                y2new = y2 + random_array[3] * np.sqrt(constant * time_step)
                x1,y1,x2,y2 = x1new,y1new,x2new,y2new
                out.write("{} {} {} {} {}\n".format(time_step*(i+1), x1, y1, x2, y2))
            update_progress(j / (runs))
            out.close()

# update_progress() : Displays or updates a console progress bar
# Accepts a float between 0 and 1. Any int will be converted to a float.
# A value under 0 represents a 'halt'.
# A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

if __name__ == '__main__':

    '''
    To initialise the program pass the following arguments after the name of the .py file
    
    steps: Number of steps the random walk will take
    runs: How many times the walk will be carried out
    constant: Constant used which is equal to kT/6pi*viscocity*radius
    time step: How big of a time step the algorithm will take
    repeats: How many times will the program run (Used because there is a limitation on the size of an array)
    initial separation: The initial separation of the two spheres
    '''
    steps = int(config.steps)
    runs = int(config.runs)
    constant = float(config.constant)
    time_step = float(config.time_step)
    init_separation = float(config.init_separation)
    hydro = str(config.hydro)
    file = open("Results.out","a")
    try:
        os.mkdir("Walk hydro: {} init_sep: {} constant: {} time_step: {} steps: {} runs: {}".format(hydro,init_separation,
            constant,time_step,steps,runs))
    except OSError, e:
        if e.errno != 17:
            raise
    os.chdir("Walk hydro: {} init_sep: {} constant: {} time_step: {} steps: {} runs: {}".format(hydro,init_separation,
            constant,time_step,steps,runs))
    walk()
