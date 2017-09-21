from __future__ import division

import sys
import os
import math
import time
import linecache
import config
import argparse
import numpy as np

"""
Version 2
Simulate run for a number of constants and find maximum separation squared.
"""

def walk(k,max_file):
    """
    Performs the random walk.
    There are 4 possibilities of a walk, which are combinations of the binary states of hydrodynamic
    interactions and thermal noise.
    k is the weissenberg number
    """
    # No hydrodynamic interaction and no noise
    if str.upper(hydro) in ("NO","N","NO HYDRO") and str.upper(noise) in ("NO","N","NO NOISE"):
        for j in xrange(runs):
            out = open("Run{}_c{}".format(j,k), "w")
            z = zinitial
            x = xinitial
            y = yinitial
            rseparation = x * x + y * y + z*z
            out.write("{} {} {} {}\n".format(0, x, y, z))
            maxr = rseparation
            tmax=0
            for i in xrange(steps):
                # S is the stokeslet tensor
                rseparation = x * x + y * y +z*z
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
                xnew = x + k * y * time_step - (x * time_step) / (1 - rseparation)
                ynew = y - (y * time_step) / (1 - rseparation)
                znew = z - (z * time_step)/(1-rseparation)
                # x2new = x2 + np.dot(Minv, random_array)[2] * np.sqrt(
                #     constant * time_step)
                # y2new = y2 + np.dot(Minv, random_array)[3] * np.sqrt(
                #     constant * time_step)
                x, y,z = xnew, ynew , znew
                #Finding out largest separation squared
                if maxr <= x*x + y*y + z*z and args.max:
                    maxr = x*x+y*y +z*z
                    tmax = (i+1)*time_step
                out.write("{} {} {} {}\n".format(time_step * (i + 1),x, y,z))

            max_file.write("{} {} {}\n".format(k,maxr/(init_separation**2),tmax))
            update_progress(j / (runs))
            out.close()

    # Hydrodynamic interactions but not noise
    elif str.upper(hydro) in ("YES","Y","HYDRO") and str.upper(noise) in ("NO","N","NO NOISE"):
        for j in xrange(runs):
            out = open("Run{}_c{}".format(j, k), "w")
            z = zinitial
            x = xinitial
            y = yinitial
            rseparation = x * x + y * y + z*z
            out.write("{} {} {} {}\n".format(0, x, y, z))
            maxr = rseparation
            tmax = 0
            for i in xrange(steps):
                rseparation = x * x + y * y +z*z
                ynew = y + time_step*(1/(1+3/(4*math.sqrt(rseparation)*ra_ratio)))*(-y/(1-rseparation) +
                    (1/(ra_ratio*math.sqrt(rseparation)*(4/3) + 2))*(-y*y*x*k/rseparation +
                    (y*y*y + y*x*x + y*z*z)/((1-rseparation)*rseparation)))

                xnew = x + time_step * (1 / (1 + 3 / (4 * math.sqrt(rseparation) * ra_ratio))) * (k*y
                -x / (1 - rseparation) +(1 / (ra_ratio * math.sqrt(rseparation)*(4/3) + 2)) * (-y * x * x * k / rseparation
                + (x * x * x + x * y * y + x*z*z) / ((1 - rseparation) * rseparation)))

                znew = z + time_step * (1 / (1 + 3 / (4 * math.sqrt(rseparation) * ra_ratio))) * (
                        -z / (1 - rseparation) +(1 / (ra_ratio * math.sqrt(rseparation) * (4 / 3) + 2)) * (-z*y*x*k/rseparation +
                        (z*y*y + z*x*x + z*z*z) / ((1 - rseparation) * rseparation)))
                x, y, z = xnew, ynew, znew
                # Finding out largest separation squared
                if maxr <= x * x + y * y + z*z and args.max:
                    maxr = x * x + y * y +z*z
                    tmax = (i+1) * time_step
                out.write("{} {} {} {}\n".format(time_step * (i + 1), x, y, z))
            if args.max:
                max_file.write("{} {} {}\n".format(k, maxr / (init_separation ** 2), tmax))
            update_progress(j / (runs))
            out.close()

    elif str.upper(hydro) in ("NO","N","NO HYDRO") and str.upper(noise) in ("YES","Y","NOISE"):
        pass

    elif str.upper(hydro) in ("YES","Y","HYDRO") and str.upper(noise) in ("YES","Y","NOISE"):
        for j in xrange(runs):
            out = open("Run{}_c{}".format(j, k), "w")
            z = zinitial
            x = xinitial
            y = yinitial
            out.write("{} {} {} {}\n".format(0, x, y, z))

            for i in xrange(steps):
                rseparation = x * x + y * y +z*z
                rvector = np.array([x,y,z])
                #Matrix of stokeslet is created
                Mxx = 1- (3/(4*math.sqrt(rseparation) * ra_ratio))*(1 + x*x/rseparation)
                Mxy = - 3/(4*math.sqrt(rseparation) * ra_ratio)*(x*y/rseparation)
                Mxz = - 3/(4*math.sqrt(rseparation) * ra_ratio)*(x*z/rseparation)
                Myy = 1- (3/(4*math.sqrt(rseparation) * ra_ratio))*(1 + y*y/rseparation)
                Myz = - 3/(4*math.sqrt(rseparation) * ra_ratio)*(y*z/rseparation)
                Mzz = 1- 3/(4*math.sqrt(rseparation) * ra_ratio)*(1 + z*z/rseparation)
                Mmatrix = np.array([[Mxx, Mxy,Mxz],[Mxy,Myy,Myz],[Mxz,Myz,Mzz]])
                noise_vector = noise_producer(Mmatrix)
                xnew = x + k*y*time_step - time_step*Mmatrix[0].dot(rvector)/(1-rseparation) + noise_vector[0]
                ynew = y - time_step * Mmatrix[1].dot(rvector) / (1 - rseparation) + \
                       noise_vector[1]
                znew = z - time_step * Mmatrix[2].dot(rvector) / (1 - rseparation) + \
                       noise_vector[2]
                x,y,z = xnew,ynew,znew
                out.write("{} {} {} {}\n".format(time_step * (i + 1), x, y, z))
            update_progress(j / (runs))
            out.close()

def noise_producer(matrix):
    """
    This method calculates the noise exerted on the 2 spheres and returns a 3 dimensional vector.
    :param matrix: Stokeslet matrix
    :return: Random noise vector
    """
    randgaussnumbers = math.sqrt(chi* time_step)*np.random.randn(3)
    psi1 = math.sqrt(matrix[0,0])*randgaussnumbers[0]

    print (matrix[0,0]*matrix[1,1] -matrix[0,1]**2)
    psi2 = matrix[0,1]/math.sqrt(matrix[0,0])*randgaussnumbers[0] +(
            (math.sqrt(matrix[0,0]*matrix[1,1] -matrix[0,1]**2)/matrix[0,0])*randgaussnumbers[1] )
    psi3 = (matrix[0,2]/math.sqrt(matrix[0,0])) *randgaussnumbers[0] +(
        (matrix[1,2]*matrix[0,0] - matrix[0,1]*matrix[0,2])/math.sqrt(matrix[1,1]*matrix[0,0] - matrix[0,1]**2))*(
        randgaussnumbers[1]) +math.sqrt(
        matrix[2,2] -(matrix[0,2]/math.sqrt(matrix[0,0]))**2 -
        ((matrix[1, 2] * matrix[0, 0] - matrix[0, 1] * matrix[0, 2]) / math.sqrt(
                matrix[1, 1] * matrix[0, 0] - matrix[0, 1] ** 2))**2) * randgaussnumbers[2]

    return np.array([psi1,psi2,psi3])
def analyse():
    """
    This function goes to the folder of the previously run simulation, and averages over
    all the runs (which, because of thermal noise, are different). It writes the mean
    squared separation and mean separation of each constant simulated in a file. It also
    creates a files with the maximum separation of every constant and at which time it occured
    """

    os.chdir("Shear_constants:{}-{}_numc:{}_hydro:{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0],
    constant[-1],len(constant),hydro, steps,time_step, ra_ratio,noise))
    max_file = open("AMax_File_con"
                    ":{}-{}_numc:{}_h{}_s{}_ts{}_ra{}_n{}".format(
        constant[0], constant[-1], len(constant), hydro, steps, time_step, ra_ratio, noise), "w")
    max_file.write("#Constant, Maxseparation over initial separation, time of max separation\n")

    for v in constant:
        #Rfile = open(os.getcwd() + "/Results.out", "a")

        out = open("MSS_{}_{}_{}_{}_{}_{}.out".format(v,hydro, steps,time_step,ra_ratio,noise), "w")
        nout = open("MS_{}_{}_{}_{}_{}_{}.out".format(v,hydro, steps,time_step,ra_ratio,noise), "w")

        #The following way of reading files is because of a limitation
        #in the number of files a computer can have open at the same time
        thousands_of_runs = int(math.ceil(runs / 1000))
        ms_list = []
        mss_list = []
        # Reads every thousand runs of a simulation
        for k in range(thousands_of_runs):
            # Opens the first 1000 runs in a dictionary, then opens the next 1000 and so on.
            filedata = {i: open("Run{}_c{}".format(i,v), "r") for i in xrange(k * 1000, min(runs, (k + 1) * 1000))}
            # Mean separation and Mean square separation lists that contain temporary files
            # with the respective values for every thousand runs. They are deleted afterwards
            ms_list.append(open("ms_{}th_thousand.tmp".format(k), "w"))
            mss_list.append(open("mss_{}th_thousand.tmp".format(k), "w"))

            # Adding squared separation and separation together
            # to average noise
            for lines in xrange(steps + 1):
                s = 0
                ssq = 0

                for file in filedata.values():
                    token = str.split(file.readline())
                    # This convenion will most likely change in the 3rd version of the program
                    t = float(token[0])
                    x = float(token[1])
                    y = float(token[2])
                    z = float(token[3])
                    rsepparation = x*x + y*y + z*z

                    s += rsepparation
                    ssq += math.sqrt(rsepparation)
                mss_list[k].write("{} {}\n".format(t, s / runs))
                ms_list[k].write("{} {}\n".format(t, (ssq / runs)))
                update_progress(lines / (steps))
            for fruns in filedata.values():
                fruns.close()
            ms_list[k].close()
            mss_list[k].close()
            ms_list[k] = open("ms_{}th_thousand.tmp".format(k), "r")
            mss_list[k] = open("mss_{}th_thousand.tmp".format(k), "r")

        # This loop goes through the temporary file in ms_list and mss_list and finds the
        # largest sepparation. It also finds the mean separation and separation squared if
        # the number of runs was more than 1000. If its under 1000 runs then this loop will
        # slow down the computation by a bit.
        # ~~~~~~~~~ NOTE: If computation time is an issue then modify this ~~~~~~~~~~~~~~~~~~~~~~~~~~
        print "~~~~~~~Merging and finding Max value~~~~~~~~~"
        maxr = 0
        tmax = 0
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
            if maxr <= mean_mss:
                maxr = mean_mss
                tmax = t
        # Max separation squared over initial separation squared is stored in a max file for
        # every constant
        # The loop deletes the unnecessary temporary files
        max_file.write("{} {} {}\n".format(v, maxr / (init_separation ** 2), tmax))
        for k in range(thousands_of_runs):
            os.remove(mss_list[k].name)
            os.remove(ms_list[k].name)
        out.close()
        nout.close()
        meansqsep = float(str.split(linecache.getline(out.name, steps + 1))[1])
        meansep = float(str.split(linecache.getline(nout.name, steps + 1))[1])
        #print("Mean squared separation over {} runs: {} ".format(runs, meansqsep))
        #print("Root Mean squared separation over {} runs: {} ".format(runs, math.sqrt(meansqsep)))
        #print ("Mean separation over {} runs : {}".format(runs, meansep))
        # Appending the results at the end of the Results.out file
        # Rfile.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~R E S U L T S~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
        # Rfile.write(
        #     "Time-Step:{} Steps:{} runs:{} constant:{} Initial separation:{} hydro: {} time&date: {} \n".format(time_step,
        #                                                                                                         steps, runs,
        #                                                                                                         constant[j],
        #                                                                                                         init_separation,
        #                                                                                                         hydro,
        #                                                                                                         time.strftime(
        #                                                                                                             "%c")))
        # Rfile.write("Mean Squared separation {}\n".format(meansqsep))
        # Rfile.write("Root Mean squared separation {}\n".format(math.sqrt(meansqsep)))
        # Rfile.write("Mean separation {}\n".format(meansep))
        # Rfile.close()
        # Mean squared displacement. Each row has a colour of the rainbow.
        # if args.walk:
        #     os.chdir("Max_Separation_constants:{}-{}_numc:{}".format(constant[0],constant[-1],len(constant)))
        #     max_file = open("Max_Separation_constants:{}-{}_numc:{}".format(constant[0],constant[-1],len(constant)),"w")
        #     for n,j in enumerate(constant):


def simulate():
    """
    For version 1 the creation of folder is inside the for loop.
    For version 2 all simulations are created in a single folder but separate files.
    :return:
    """


    #This code simulates the walk of the polymer and stores the max separation in a file

    try:
        os.mkdir("Shear_constants:{}-{}_numc:{}_hydro:{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0],constant[-1],
        len(constant),hydro,steps,time_step,ra_ratio,noise))
    except OSError as e:
        if e.errno != 17:
            raise
    os.chdir("Shear_constants:{}-{}_numc:{}_hydro:{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0],constant[-1],
        len(constant),hydro,steps,time_step,ra_ratio,noise))
    if args.max:
        max_file = open("Max_File_con"
                        ":{}-{}_numc:{}_h{}_s{}_ts{}_ra{}_n{}".format(
            constant[0],constant[-1],len(constant),hydro,steps,time_step,ra_ratio,noise),"w")
        max_file.write("#Constant, Maxseparation over initial separation, time of max separation\n")
    else:
        max_file = None
    for n,j in enumerate(constant):

        walk(j,max_file)
        print("wi number {} done, {} left".format(j,len(constant)-n))

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

if __name__ == "__main__":
    """
    Program gets variables from config file. 
    To run the program, arguments have to be passed while initialising it in the form of flags.
    The -walk argument runs the random walk simulation using the parameters of the config file
    The -analyse argument finds the folder that matches the parameters of the config file
    and subsequently runs through the files to produce a mean separation squared and max separation
    files for every constant.
    The -max argument is used only when there is no noise in the simulation so as to avoid running
    the computationally expensive -analyse method 
    """
    steps = int(config.steps)
    runs = int(config.runs)
    constant = config.constant #Weissenberg numbers in a numpy array
    time_step = float(config.time_step)
    yinitial= float(config.yinitial)
    xinitial = float(config.xinitial)
    zinitial = float(config.zinitial)
    init_separation = math.sqrt(xinitial**2 + yinitial**2 + zinitial**2)
    hydro = str(config.hydro)
    ra_ratio = float(config.ra_ratio)
    noise = str(config.noise)
    chi = float(config.chi)
    parser = argparse.ArgumentParser(description="Program version 2"
                                                 "The program simulates the motion of a polymer in shear flow.\n "
                                                 "The model is of a finite extensibility non-linear elastic spring"
                                                 "(FENE).\n"
                                                 "Parameters of the simulation can be found in the config.py file.\n"
                                                 "\nVersion 1: Has no thermal fluctuations so analyser doesn't do anything\n"
                                     "Version 2: There are hydrodynamic interactions between the two ends of the polymer chain.\n"
                                     "Version 3: There is thermal noise acting on the polymer")
    # parser.add_argument("echo", help="echo the string you use here")
    parser.add_argument("-a", "--analyse", help="Run analyser. For version 1 this does nothing (v1) since there "
                                                "are no thermal fluctuations.",
                        action="store_true")
    parser.add_argument("-w", "--walk", help="Simulate walks with parameters of config.py file", action="store_true")

    parser.add_argument("-m", "--max", help="Store max separation squared in a separate file", action="store_true")


    args = parser.parse_args()
    # print args.echo
    if args.walk:
        simulate()
    if args.analyse:
        analyse()
