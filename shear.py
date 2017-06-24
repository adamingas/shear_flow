from __future__ import division

import sys
import os
import math
import time
import linecache
import config
import argparse
#import numpy as np


# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def walk(k):
    """
    Performs the random walk. Creates the Matrix with the coefficients of the velocity profile created when
    one of the spheres moves through the fluid.
    """
    # if str.upper(hydro) in ("YES","Y","HYDRO"):
    for j in xrange(runs):
        out = open("Run{}".format(j), "w")

        x = 0
        y = init_separation
        rseparation = x * x + y * y
        out.write("{} {}\n".format(0, rseparation))
        for i in xrange(steps):
            # S is the stokeslet tensor
            rseparation = x * x + y * y
            # Sxy = (3/(4*math.sqrt(distance_sq)))*(
            # (x2-x1)*(y2 - y1)/
            #              distance_sq)
            # Sxx = (3/(4*math.sqrt(distance_sq)))*(1+
            # (x2-x1)*(x2-x1)/distance_sq)
            # Syy =  (3/(4*math.sqrt(distance_sq)))*(1+
            # (y2 - y1)*(y2-y1)/distance_sq)
            #
            # Matrix = np.array([[1,0,-Sxx,-Sxy],[0,1,-Sxy,-Syy],[-Sxx,-Sxy,1,0],[-Sxy,-Syy,0,1]])
            # Minv = np.linalg.inv(Matrix)
            # random_array = np.random.randn(4)
            xnew = x + constant[k] * y * time_step - (x * time_step) / (1 - rseparation)
            ynew = y - (y * time_step) / (1 - rseparation)
            # x2new = x2 + np.dot(Minv, random_array)[2] * np.sqrt(
            #     constant * time_step)
            # y2new = y2 + np.dot(Minv, random_array)[3] * np.sqrt(
            #     constant * time_step)
            x, y = xnew, ynew
            out.write("{} {} {} {}\n".format(time_step * (i + 1), (x * x + y * y), x, y))
        update_progress(j / (runs))
        out.close()
    # elif str.upper(hydro) in ("NO","N","NOHYDRO") :
    #
    # for j in xrange(runs):
    #     out = open("Run{}".format(j), "w")
    #
    #     x1, y1, x2 = 0,0,0
    #     y2 = init_separation
    #     out.write("{} {} {} {} {}\n".format(0,x1,y1,x2,y2))
    #     for i in xrange(steps):
    #         #S is the stokeslet tensor
    #
    #         random_array = np.random.randn(4)
    #         x1new= x1 + random_array[0]*np.sqrt(constant * time_step)
    #         y1new = y1+ random_array[1] * np.sqrt(constant * time_step)
    #         x2new = x2 + random_array[2] * np.sqrt(constant * time_step)
    #         y2new = y2 + random_array[3] * np.sqrt(constant * time_step)
    #         x1,y1,x2,y2 = x1new,y1new,x2new,y2new
    #         out.write("{} {} {} {} {}\n".format(time_step*(i+1), x1, y1, x2, y2))
    #     update_progress(j / (runs))
    #     out.close()


def update_progress(progress):
    barLength = 10  # Modify this to change the length of the progress bar
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
    block = int(round(barLength * progress))
    text = "\rPercent: [{0}] {1}% {2}".format("#" * block + "-" * (barLength - block), progress * 100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def analyse(j):
    '''
    To initialise the program pass the following arguments after the name of the .py file
    
    steps: Number of steps the random walk will take
    runs: How many times the walk will be carried out
    constant: Constant used which is equal to kT/6pi*viscocity*radius
    time step: How big of a time step the algorithm will take
    repeats: How many times will the program run (Used because there is a limitation on the size of an array)
    initial separation: The initial separation of the two spheres
    '''
    while j<len(constant):
        Rfile = open(os.getcwd() + "/Results.out", "a")
        os.chdir(
            "Shear_Walk hydro: {} init_sep: {} constant: {} time_step: {} steps: {} runs: {}".format(hydro, init_separation,
                                                                                                     constant[j],
                                                                                                     time_step, steps,
                                                                                                     runs))
        out = open("MSS_{}_{}_{}_{}_{}_{}.out".format(hydro, init_separation, constant[j], time_step, steps, runs), "w")
        nout = open("MS_{}_{}_{}_{}_{}_{}.out".format(hydro, init_separation, constant[j], time_step, steps, runs), "w")
        ssq = 0
        s = 0
        # run_files = ["Run{}".format(x) for x in range(runs)]
        # filedata = {run: open(run, "r") for run in run_files}
        # filedata2 = {run: open(run, "r") for run in run_files}
        thousands_of_runs = int(round(runs / 1000))
        ms_list = []
        mss_list = []
        for k in range(thousands_of_runs):
            filedata = {i: open("Run{}".format(i), "r") for i in xrange(k * 1000, min(runs, (k + 1) * 1000))}
            ms_list.append(open("ms_{}th_thousand.tmp".format(k), "w"))
            mss_list.append(open("mss_{}th_thousand.tmp".format(k), "w"))
            for j in xrange(steps + 1):
                s = 0
                ssq = 0

                for file in filedata.values():
                    token = str.split(file.readline())

                    t = float(token[0])
                    x1 = float(token[1])
                    y1 = float(token[2])
                    x2 = float(token[3])
                    y2 = float(token[4])
                    s += (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)
                    ssq += math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
                mss_list[k].write("{} {}\n".format(t, s / runs))
                ms_list[k].write("{} {}\n".format(t, (ssq / runs)))
                update_progress(j / (steps))
            for fruns in filedata.values():
                fruns.close()
            ms_list[k].close()
            mss_list[k].close()
            ms_list[k] = open("ms_{}th_thousand.tmp".format(k), "r")
            mss_list[k] = open("mss_{}th_thousand.tmp".format(k), "r")
        for j in xrange(steps + 1):
            mean_mss = 0
            mean_ms = 0
            for k in range(thousands_of_runs):
                mstoken = str.split(ms_list[k].readline())
                msstoken = str.split(mss_list[k].readline())
                t = float(mstoken[0])
                mssn = float(msstoken[1])
                msn = float(mstoken[1])
                mean_mss += mssn
                mean_ms += msn
            out.write("{} {}\n".format(t, mean_mss))
            nout.write("{} {}\n".format(t, mean_ms))
        for k in range(thousands_of_runs):
            os.remove(mss_list[k].name)
            os.remove(ms_list[k].name)
        out.close()
        nout.close()
        meansqsep = float(str.split(linecache.getline(out.name, steps + 1))[1])
        meansep = float(str.split(linecache.getline(nout.name, steps + 1))[1])
        print("Mean squared separation over {} runs: {} ".format(runs, meansqsep))
        print("Root Mean squared separation over {} runs: {} ".format(runs, math.sqrt(meansqsep)))
        print ("Mean separation over {} runs : {}".format(runs, meansep))
        # Appending the results at the end of the Results.out file
        Rfile.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~R E S U L T S~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
        Rfile.write(
            "Time-Step:{} Steps:{} runs:{} constant:{} Initial separation:{} hydro: {} time&date: {} \n".format(time_step,
                                                                                                                steps, runs,
                                                                                                                constant[j],
                                                                                                                init_separation,
                                                                                                                hydro,
                                                                                                                time.strftime(
                                                                                                                    "%c")))
        Rfile.write("Mean Squared separation {}\n".format(meansqsep))
        Rfile.write("Root Mean squared separation {}\n".format(math.sqrt(meansqsep)))
        Rfile.write("Mean separation {}\n".format(meansep))
        Rfile.close()
        # Mean squared displacement. Each row has a colour of the rainbow.
        os.chdir('..')
        j = j +1

def simulate():
    j = 0
    while j < len(constant):
        try:
            os.mkdir(
                "Shear_Walk hydro: {} init_sep: {} constant: {} time_step: {} steps: {} runs: {}".format(hydro,
                                                                                                         init_separation,
                                                                                                         constant[j],
                                                                                                         time_step,
                                                                                                         steps, runs))
        except OSError, e:
            if e.errno != 17:
                raise
        os.chdir("Shear_Walk hydro: {} init_sep: {} constant: {} time_step: {} steps: {} runs: {}".format(hydro,
                                                                                                          init_separation,
                                                                                                          constant[j],
                                                                                                          time_step,
                                                                                                          steps, runs))
        walk(j)
        os.chdir('..')
        j = j + 1


if __name__ == "__main__":

    steps = int(config.steps)
    runs = int(config.runs)
    constant = config.constant
    time_step = float(config.time_step)
    init_separation = float(config.init_separation)
    hydro = str(config.hydro)
    parser = argparse.ArgumentParser(description="Program version 1"
                                                 "The program simulates te motion of a polymer in shear flow.\n "
                                                 "The model is of a finite extensibility non-linear elastic spring"
                                                 "(FENE).\n"
                                                 "Parameters of the simulation can be found in the config.py file.\n"
                                                 "\n Version 1: Has no thermal fluctuations so analyser doesn't do anything")
    # parser.add_argument("echo", help="echo the string you use here")
    parser.add_argument("-a", "--analyse", help="Run analyser. For version 1 this does nothing (v1) since there "
                                                "are no thermal fluctuations.",
                        action="store_true")
    parser.add_argument("-w", "--walk", help="Simulate a walk with parameters of config.py file", action="store_true")

    args = parser.parse_args()
    # print args.echo
    if args.walk:
        simulate()
    if args.analyse:
        analyse(0)
