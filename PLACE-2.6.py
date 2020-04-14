#!/usr/bin/env python2.7

#Project start 10/11/11 - Last updated 2/15/12 
#Place Correction in Python V2.6 - By Nicholas Mei
#Goal: Does Place Correction for specified scan days

# 2/15/12: Major overhaul of interface and usage of AFNIpyIO instead of nibabel
# 2/16/12: Added some error handling and warnings
# 2/29/12: Added command line callability, if no options are passed then PLACE code defaults to using Tkinter gui

#======================== required modules to be installed =================
# Requires numpy, and scipy to be installed as Python modules.
# Num.py - http://numpy.scipy.org/
# Sci.py - http://www.scipy.org/
# Also requires AFNIpyIO module to load in .HEAD/.BRIK files

#======================== command line options ========================
#usage: PLACE-2.6.py [-h] -d DSETS [DSETS ...] -p PARSCAN -m DMAP
#                     [-s [SAVE_PATH]]
#Do PLACE correction
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -d DSETS [DSETS ...], -dset DSETS [DSETS ...], --dset DSETS [DSETS ...]
#                        path(s) to .HEAD files (separated by spaces)
#  -p PARSCAN, -pscan PARSCAN, --pscan PARSCAN
#                        path to the ParScan file
#  -m DMAP, -dmap DMAP, --dmap DMAP
#                        path to the Dmap file
#  -s [SAVE_PATH], -save [SAVE_PATH], --save [SAVE_PATH]
#                        optional save path (default: same as first -dset path)

# example call:
#PLACE-2.6.py -d /Users/mein/Desktop/AFNI_Files/LY/111209/EPI_LY_111209_E07+orig.HEAD\
#		/Users/mein/Desktop/AFNI_Files/LY/111209/EPI_LY_111209_E11+orig.HEAD\
#	     -p /Users/mein/Desktop/AFNI_Files/LY/place_test2/ParScan\
#             -m /Users/mein/Desktop/AFNI_Files/LY/place_test2/Dmap\
#	     -s /Users/mein/Desktop/AFNI_Files/Ly/place_test2/
# (the extra empty line after the final argument is important!! don't forget it or calling PLACE-2.6.py won't work!)

import sys
import os
import time
import string
import copy

import Tkinter as tk
import tkFileDialog

import subprocess as sp
import numpy as np
#custom python module to read in and write out AFNI files
from afnipyio import AFNIPyIO as afni

from scipy.sparse import *
from operator import itemgetter
from contextlib import contextmanager

import argparse

root = tk.Tk()

#Version so generated logs and scripts are traceable to a specific version of the code (especially important if there are errors)
version = 'V2.6'

#==================================== CUSTOM USER SETTINGS CHANGE BEFORE YOU USE FOR YOURSELF ========================

defaultdir = '/nfs/data42/seidlitzjm'

#====================================================================================================================

#this will change once you've selected a Dmap or Parscan file
dmap_parscan_defaultdir = defaultdir
loglist = []

#Error Handler
class Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

#========================= makeunwarpmatrix function translation from matlab code ===============
#All credit for the logic of the following code goes to GKA and was written July 2006
#Questions about the overall logic of this function should go to (GKA) geoff.adams@gmail.com
#Questions about the num.py and sci.py conversion or reports for bugs should go to njmei11@gmail.com

def makeunwarpmatrix(dmap, expan):
    
    expan = float(expan)

    #transpose dmap to put phase in columns and read in rows
    #results in a 3200x64x28 matrix (for 11/02/25 CS scan)
    #remember!!! axes in python = matlabaxes-1
    dmap3 = np.transpose(dmap,[1,0,2])

    #print 'Transposed Dmap is: \n' + str(dmap3)

    print 'Dimensions of transposed dmap are: ' + str(dmap3.shape)

    #print dmap3.dtype
    nslice = np.size(dmap3, axis=2)             #Number of slices
    nphasevol = np.size(dmap3, axis=0)/expan    #Number of phase points in volume
    nphasedmap = np.size(dmap3, axis=0)         #Number of phase points in dmap
    nread = np.size(dmap3, axis=1)              #Number of read points

    print 'nslice is: ' + str(nslice)
    print 'nphasevol is: ' + str(nphasevol)
    print 'nphasedmap is: ' + str(nphasedmap)
    print 'nread is: ' + str(nread)

    crop = np.array([1,1,nread,nphasevol])

    #print str(crop)
    
    #------------- Convert from displacement mapping to absolute indexing ---------------------  
    #Matlab code: phasePoints = repmat([1:nPhaseDmap]', [1,nRead,nSlice]);
    #Python indexes start at 0 but we will do the "offset" later to preserve parity with matlab code
    a = np.arange(1,nphasedmap+1)                               # a = [1:nPhaseDmap]
    #print a

    #tile promotes 'a' to the dimension of the tile argument by prepending new axes as necessary
    #matlab seems to work in the opposite way...
    #fix is to append [:,np.newaxis,np.newaxis] to force the matlab tiling method
    phasepoints = np.tile(a[:,np.newaxis,np.newaxis], [1,nread,nslice])

    #print str('Phasepoints are: ' + str(phasepoints))    
    print 'Dimension of phasepoints variable is: ' + str(np.shape(phasepoints))

    phaseorigin = np.floor(phasepoints-dmap3, dtype=np.double)

    #print np.max(phaseorigin)
    #print np.min(phaseorigin)
    #print str('Phase origins are: ' + str(phaseorigin))
    #print str(np.size(phaseorigin))
    print 'Dimensions of phaseorigin variable is: ' + str(np.shape(phaseorigin)) + '\n'

    #-------------------- Correct wrapping errors --------------------------------------
    #create a boolean index 'outofbound' where phaseorigin values are > nphasedmap
    outofbound = phaseorigin > nphasedmap
    #print 'First out of bound' + str(outofbound)
    phaseorigin[outofbound] = phaseorigin[outofbound] - nphasedmap
    
    outofbound = phaseorigin < 1
    #print 'Second out of bound' + str(outofbound)
    phaseorigin[outofbound] = phaseorigin[outofbound] + nphasedmap

    #----------------------- Cropping --------------------------------------
    x1 = crop[0]                        #Read start     x1 evaluates as 1 when crop = [1, 1, 64, 32]
    x2 = x1 + crop[2]-1                 #Read end       
    y1 = (crop[1]-1)*expan+1            #Phase start    
    y2 = (crop[1] + crop[3]-1)*expan    #Phase end

    #phaseOrigin = phaseOrigin(y1:y2,x1:x2,:);
    #nRead = crop(3);
    #nPhaseVol = crop(4);
    #nPhaseDmap = nPhaseVol*expan;

    #cropping phaseorigin to "correct" dimensions. Without user input crop value, it's kind of useless...
    #slicing in python is "inclusive" i.e. cropping anything below 0 and above 3200 requires [0:3200]
    phaseorigin = phaseorigin[int(y1-1):int(y2), int(x1-1):int(x2)]
    nread = crop[2]
    nphasevol = crop[3]
    nphasedmap = nphasevol*expan

    phaseorigin = phaseorigin -y1 + 1
    
    outofbound = phaseorigin > nphasedmap
    phaseorigin[outofbound] = nphasedmap
    
    outofbound = phaseorigin < 1
    phaseorigin[outofbound] = 1

    #print np.max(phaseorigin)
    #print np.min(phaseorigin)

    #--------------------------- Reshaping phaseorigin into a vector to construct warping matrix ---------------
    #reshape should be in 'Fortran' mode because that's how Matlab does it
    phaseorigin = np.reshape(phaseorigin, [nphasedmap,-1], order='F')

    #Matlab code: phaseOrigin = phaseOrigin + repmat([0:nRead*nSlice-1]*nPhaseDmap, [nPhaseDmap, 1]);

    b = np.arange(0,nread*nslice)
    b = b*nphasedmap

    rep = np.tile(b, [nphasedmap, 1])

    phaseorigin = phaseorigin + rep

    #print np.max(phaseorigin)
    #print np.min(phaseorigin)

    #Matlab: phaseOrigin = ceil(phaseOrigin(:)/expan);

    phaseorigin = np.ceil(phaseorigin.flatten(1)/expan)
    #Cludgy fix for problem of max(phaseorigin) being 57344 instead of 57343, not sure if this fix is "correct"
    #update (10/20/11): resulting unwarpmatrix matches output of matlab program (after accounting for 0 vs. 1 index differences)
    phaseorigin = phaseorigin -1

    #print np.max(phaseorigin)
    #print np.min(phaseorigin)

    numpoints = nread*nslice*nphasevol
    #this code in Matlab doesn't do anything....
    #Matlab: points = range(1,numpoints+1)

    #------------------------------------ Final Step! ---------------------------------------
    #Matlab sparse takes: S = sparse(i,j,s,m,n)
    #Generates an m by n sparse matrix such that: S(i(k),j(k)) = s(k)
    #Matlab matrix type is CSC matrix

    #scipy csc matrix takes: csc_matrix((data, ij), [shape=(M, N)])
    #a[ij[0,k],ij[1,k]] = data[k]
    #i would be phaseOrigin variable
    
    #Matlab code is: unwarpMatrix = sparse(phaseOrigin, ceil([1:nRead*nSlice*nPhaseDmap]/expan), 1, numPoints, numPoints)/expan;
    size = nread*nslice*nphasedmap

    #should use this instead because csc_matrix requires largest values in phaseorigin and j to be <= numpoints
    #essentially j and phaseorigin are indexes and because python starts indexes from 0 we should correct for that
    #instead of doing something like: ceil([1:nRead*nSlice*nPhaseDmap]/expan) 
    j = np.floor(np.arange(0,size, dtype=np.double)/expan)

    #Matlab apparently treats '1' as a scalar so I should be tiling 1 to the same size as j and phaseorigin
    s = np.tile(1,size)
    
    unwarpmatrix = csc_matrix((s,(phaseorigin, j)), shape=(int(numpoints),int(numpoints)))/expan
    return  unwarpmatrix

#========================= Unwarp a volume given an unwarp matrix ======================
def placeunwarp(vol, unwarpmatrix):
    #Written by GKA
    #Extracted from Blockana by dbtm (matlab comment)

    #transpose matrix such that columns 1,2,3 (x, y, z) --> 2,1,3 (y, x, z)
    vol = np.transpose(vol, [1,0,2])

    origshape = np.array(np.shape(vol))

    #turn volume into linear vector (array)
    flatvol = vol.flatten(1)
    #flatvol = vol.flatten(order='F')

    #heart of the correction! Dot product of the unwarpmatrix vector and flatvol vector
    #sparsematrix.dot() is very difference from np.dot() or np.inner() the first one works, the second two don't
    dot = unwarpmatrix.dot(flatvol)
    
    #reshape back into original volume shape
    corvol = np.reshape(dot, origshape, order='F')

    #tranpose matrix back to 1,2,3
    corvol = np.transpose(corvol, [1,0,2])

    return corvol

def printandlog(some_string):
    global loglist

    some_string = str(some_string)
    
    print some_string
    loglist.append(some_string + '\n')

#code for benchmarking purposes
@contextmanager  
def measureTime(title):
    t1 = time.time()
    yield
    t2 = time.time()
    printandlog('%s: Generating unwarp sparse matrix took %0.2f seconds ' % (title, t2-t1))

#function to write out a logfile
def writelog(logfilepath, loglist):
    print '\n' + 'Writing out analysis log!' + '\n'
    
    from time import gmtime, strftime
    time = strftime("%Y-%m-%d %H:%M:%S", gmtime())

    #Opening our file!
    LOG = open(logfilepath, 'w')

    #Writing the header and time
    LOG.write('#PLACE correction ' + version + ' by Nicholas Mei\n')
    LOG.write('#This log was automatically generated on ' + time + ' GMT' + '\n' + '\n')

    #Dump everything collected in the log into our text file!
    for line in loglist:
        LOG.write(line)

    #Close our log!
    LOG.close()

    print 'log file written to' + logfilepath + '\n' + '\n' + 'Finally done with this condition!!' + '\n'

    #empty loglist and reinitialize it for next run
    del loglist[:]

def main():

    #--------------------------- PLACE correction function ----------------------

    def placecor(gui = True, final_dset_list = None, dmap_path = None, parscan_path = None, save_dir = None):
        if gui:
            
            final_dset_list = list(dset_list.get(0, tk.END))
            dmap_path = dmap_list.get(0,tk.END)
            parscan_path = parscan_list.get(0, tk.END)
            save_dir = save_list.get(0, tk.END)

        if not final_dset_list:
            printandlog("You did not select any datasets to PLACE correct!")
        if not dmap_path:
            printandlog("You did not select a Dmap file to use!")
        if not parscan_path:
            printandlog("You did not select a Parscan file to use!")
        if not save_dir:
            printandlog("You did not select a directoy to save PLACE corrected files to!")

        if (final_dset_list and dmap_path and parscan_path and save_dir):
            
            dmap_path = dmap_path[0]
            parscan_path = parscan_path[0]
            if gui:
                save_dir = save_dir[0]
            
            printandlog("Everything checks out, you selected the following: \n")
            printandlog("Dataset list: " + str(final_dset_list) + "\n")
            printandlog("Dmap Location: " + dmap_path)
            printandlog("Parscan Location: " + parscan_path)
            printandlog("Save Directory: " + save_dir + "\n")

            #===================================== Parscan preparation ========================================
            #let's parse our Parscan!
            try:
                rawparscan = open(parscan_path, 'r')
            except:
                raise Error('ParScan file was unable to be opened!')

            rawparscanstring = rawparscan.read()
            rawparscan.close()

            #separate all values from the raw string by ' ' (space)
            placepars = rawparscanstring.split()
            #take first 5 ParScan parameters (xres, yres, zres, reps, expansion)
            placepars = placepars[:5]

            placepars = [int(i) for i in placepars]

            printandlog('ParScan values of interest are (xres, yres, zres, reps, expansion): ' + str(placepars))

            unwarpnumread = placepars[0]
            unwarpnumphase = placepars[1]
            unwarpnumslice = placepars[2]
            unwarpnumreps = placepars[3]
            unwarpexpan = placepars[4]

            #===================================== Dmap preparation ===========================================
            #let's open up our Dmap!  
            try:
                binarydmap = open(dmap_path,'r')
            except:
                raise Error('Dmap file was unable to be opened!')

            binarydmapstring = binarydmap.read()
            binarydmap.close()

            #use the numpy fromstring() function to convert the binary file
            dmap = np.fromstring(binarydmapstring, dtype='int16')

            #loglist.append('Dmap is:' + str(dmap) + '\n')
            #print 'Dmap is: ' + str(dmap)

            #size of Dmap should be xres*yres*expan*nslices
            printandlog('Dimensions of dmap are: ' + str(np.size(dmap)))

            #reshape should be in 'Fortran' read order mode because that's how Matlab does it...
            dmap = np.reshape(dmap, [unwarpnumread, unwarpnumphase*unwarpexpan, unwarpnumslice], order='F')

            #print 'Reshaped Dmap is: \n' + str(dmap2)
            printandlog('Dimensions of reshaped dmap are: ' + str(dmap.shape))

            with measureTime('makeunwarpmatrix'):

                unwarpmatrix = makeunwarpmatrix(dmap, unwarpexpan)

            printandlog('\nSample of unwarp sparse matrix: \n \n' + str(unwarpmatrix) + '\n' + '\n')

            #========================================= Reading in data with nibabel and doing actual volume correction =========================

            time1 = time.time()
            skipped = 0

            #go through all epis in our list
            for epi in final_dset_list:

                printandlog('\nProcessing scan: ' + str(epi))
                
                #load img1 using nibabel
                img = afni.load(epi)
                
                #Warning: data seems to be a pointer to image data so modifying it is not recommended
                #I.E. only reference and read values from it. Don't operate on it...
                data = copy.copy(img.brik.volume)

                printandlog('Data matrix shape is: ' + str(np.shape(data)))

                #ensure that the unwarp matrix and data matrix have the same shape otherwise bad stuff will happen...
                if tuple(np.shape(data)[0:3]) == (unwarpnumread, unwarpnumphase, unwarpnumslice):

                    #initialize our new img matrix (arrays)
                    tempimg = np.empty(np.shape(data), dtype=img.head.dtype)

                    #go through each timepoint: np.size(data, axis=3) could also be replaced by numtimes
                    for t in range(np.size(data, axis=3)):

                        tempimg[:,:,:,t] = placeunwarp(data[:,:,:,t], unwarpmatrix)

                    tempimg = tempimg[np.arange(0,unwarpnumread),:,:,:]

                    tempimg = tempimg[:,np.arange(0,unwarpnumphase),:,:]

                    #print(np.shape(tempimg))

                    img.brik.volume = tempimg

                    if epi.endswith('+orig.HEAD'):
                        correctedname = epi.rstrip('+orig.HEAD') + '_pc+orig.HEAD'
                    if epi.endswith('+acpc.HEAD'):
                        correctedname = epi.rstrip('+acpc.HEAD') + '_pc+acpc.HEAD'
                    if epi.endswith('+tlrc.HEAD'):
                        correctedname = epi.rstrip('+tlrc.HEAD') + '_pc+tlrc.HEAD'

                    corrected_base = os.path.basename(correctedname)

                    save_path = os.path.join(save_dir, corrected_base)

                    if os.path.exists(save_path):
                        printandlog("You already have a file of the same name as: " + save_path)
                        printandlog("PLACE correction was NOT done for: " + epi)
                        skipped += 1
                    else:
                        img.save(save_path)
                        printandlog('New PLACE corrected file will be saved to: ' + save_path)
                        printandlog('Alright!! Finished place correcting: ' + str(epi) + '\n')        
                else:
                    printandlog("Your scan: " + epi + " Does not match the dimensions of your Parscan parameters!!!")
                    printandlog("Skipping PLACE correction of: " + epi)
                    skipped += 1

            time2 = time.time()

            printandlog('\nPhew finished converting ' + str(len(final_dset_list)-skipped) + ' scans\n')
            printandlog('Entire conversion took: ' + str((time2-time1)) + ' seconds\n')

            log_name = 'PLACE_log.txt'
            log_path = os.path.join(save_dir, log_name)

            writelog(log_path, loglist)


    #check how many arguments are being passed to the PLACE-2.6.py if it is just 1 (the default) start the GUI
    if len(sys.argv) == 1:

        #function to hide Tkinter console
        def hideTkConsole(root):
            try:
                root.tk.call('console', 'hide')
            except tk.TclError:
                # Some versions of the Tk framework don't have a console object
                pass

        #------------------------- GUI functions --------------------
        def dsetchoose():
            print "Select the scans you want to run PLACE correction on, .BRIK files will be ignored"
            filelist = tkFileDialog.askopenfilenames(parent=root, title='Select the scans you want to run PLACE correction on', initialdir=defaultdir)

            if filelist:
                #dsetdir = str(os.path.dirname(filelist[0]))
                #print dsetdir   
                headnamelist = [i for i in filelist if '.HEAD' in i]
                temp_dset_list = list(dset_list.get(0, tk.END))    
                for item in headnamelist:
                    if item not in temp_dset_list:
                        dset_list.insert(tk.END, item)

        def dmapchoose():
            global dmap_parscan_defaultdir
            print "Please select the appropriate Dmap file"
            dmappath = tkFileDialog.askopenfilename(parent=root, title="Please select the appropriate Dmap file", initialdir=dmap_parscan_defaultdir)

            if dmappath:
                if dmap_list.size() > 0:
                    dmap_list.delete(0,tk.END)
                dmap_list.insert(0, dmappath)

                dmap_parscan_defaultdir = str(os.path.dirname(dmappath))

        def parscanchoose():
            global dmap_parscan_defaultdir
            print "Please select the appropriate Parscan file"
            parscanpath = tkFileDialog.askopenfilename(parent=root, title="Please select the appropriate Parscan file", initialdir=dmap_parscan_defaultdir)

            if parscanpath:
                if parscan_list.size() > 0:
                    parscan_list.delete(0,tk.END)
                parscan_list.insert(0, parscanpath)

                dmap_parscan_defaultdir = str(os.path.dirname(parscanpath))

        def savechoose():
            print "Select the directory where you wish to save PLACE corrected files to"
            savepath = tkFileDialog.askdirectory(parent=root, title='Select the directory where you wish to save PLACE corrected files to:', initialdir=defaultdir, mustexist=True)

            if savepath:
                if save_list.size() > 0:
                    save_list.delete(0,tk.END)
                save_list.insert(0, savepath)

        def deleteselection() :
                items = dset_list.curselection()
                pos = 0
                for i in items :
                    idx = int(i) - pos
                    dset_list.delete( idx,idx )
                    pos = pos + 1

        #========================= PLACE Correction Tkinter GUI =========================
        root.wm_title("Place Correction " + version)
        hideTkConsole(root)

        #------------------ Dataset and scrollbar frame -----------------
        tk.Label(root, text="Datasets to PLACE correct:").pack(side=tk.TOP)

        dset_frame = tk.Frame(root, bd=0, relief=tk.SUNKEN, padx=30)
        dset_frame.pack(fill=tk.X)

        dset_scrollbar = tk.Scrollbar(dset_frame, orient=tk.VERTICAL)
        dset_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        dset_list = tk.Listbox(dset_frame, width=80, height=8, selectmode = tk.EXTENDED, exportselection = False)
        dset_list.pack(side=tk.RIGHT, fill=tk.X, expand=True)

        dset_list.config(yscrollcommand=dset_scrollbar.set)
        dset_scrollbar.config(command=dset_list.yview)

        #------------------ Dataset buttons frame ----------------------
        dset_button_frame = tk.Frame(root, padx=30)
        dset_button_frame.pack(side=tk.TOP, fill=tk.X)

        choose_dset = tk.Button(dset_button_frame, text='Choose', command = dsetchoose)
        choose_dset.pack(side=tk.RIGHT)

        remove_dset = tk.Button(dset_button_frame, text='Remove Selected', command = deleteselection)
        remove_dset.pack(side=tk.RIGHT)

        #------------------ Dmap path selection -------------------------
        tk.Label(root, text="Select the Dmap file:").pack(side=tk.TOP)

        dmap_frame = tk.Frame(root, padx=30)
        dmap_frame.pack(side=tk.TOP, fill=tk.X)

        dmap_list = tk.Listbox(dmap_frame, width=50, height=1)
        dmap_list.pack(side=tk.LEFT, fill=tk.X, expand=True)

        choose_dmap = tk.Button(dmap_frame, text='Choose', command = dmapchoose)
        choose_dmap.pack(side=tk.RIGHT)

        #------------------ Parscan path selection ----------------------
        tk.Label(root, text="Select the Parscan file:").pack(side=tk.TOP)

        parscan_frame = tk.Frame(root, padx=30)
        parscan_frame.pack(side=tk.TOP, fill=tk.X)

        parscan_list = tk.Listbox(parscan_frame, width=50, height=1)
        parscan_list.pack(side=tk.LEFT, fill=tk.X, expand=True)

        choose_parscan = tk.Button(parscan_frame, text='Choose', command = parscanchoose)
        choose_parscan.pack(side=tk.RIGHT)

        #------------------ Save directory selection --------------------
        tk.Label(root, text="Directory to save PLACE corrected files to:").pack(side=tk.TOP)

        save_frame = tk.Frame(root, padx=30)
        save_frame.pack(side=tk.TOP, fill=tk.X)

        save_list = tk.Listbox(save_frame, width=50, height=1)
        save_list.pack(side=tk.LEFT, fill=tk.X, expand=True)

        choose_save = tk.Button(save_frame, text='Choose', command = savechoose)
        choose_save.pack(side=tk.RIGHT)

        #----------------- run button --------------------

        runit = tk.Button(root, text="Do the PLACE correction!", command = placecor)
        runit.pack(side=tk.BOTTOM, fill=tk.X, pady=15, padx=30)

        root.mainloop()

    #===================================== Command Line argument parsing =====================================

    #if the number of arguments passed is != 1 then user entered some arguments
    else:
        parser = argparse.ArgumentParser(description = 'Do PLACE correction')
        parser.add_argument('-d', '-dset', '--dset', dest = 'dsets', nargs = '+', required = True, help = 'path(s) to .HEAD files (separated by spaces)')
        parser.add_argument('-p', '-pscan', '--pscan', dest = 'parscan', nargs = 1, required = True, help = 'path to the ParScan file')
        parser.add_argument('-m', '-dmap', '--dmap', dest = 'dmap', nargs = 1, required = True, help = 'path to the Dmap file')
        parser.add_argument('-s', '-save', '--save', dest = 'save_path', nargs = '?', help = 'optional save path (default: same as first -dset path)')

        args = parser.parse_args()

        print args.dsets
        print args.parscan
        print args.dmap
        print args.save_path

        dsets_good = True
        parscan_good = True
        dmap_good = True
        save_path_good = True

        #Do some rudimentary argument checking to make sure nothing is amiss...
        for dset_path in args.dsets:
            if not os.path.exists(dset_path):
                print "The dset path you entered: " + str(dset_path) + " does not exist!!"
                dsets_good = False
            if not dset_path.endswith(".HEAD"):
                print "The dset path you entered: " + str(dset_path) + " does not end with .HEAD!!"
                dsets_good = False

        if not os.path.exists(args.parscan[0]):
            print "The ParScan path you entered: " + str(args.parscan[0]) + " does not exist!!"
            parscan_good = False
        if not os.path.basename(args.parscan[0]).startswith("ParScan"):
            print "The ParScan path you entered: " + str(args.parscan[0]) + " does not appear to point to a ParScan file!!"
            parscan_good = False

        if not os.path.exists(args.dmap[0]):
            print "The Dmap path you entered: " + str(args.dmap[0]) + " does not exist!!"
            dmap_good = False
        if not os.path.basename(args.dmap[0]).startswith("Dmap"):
            print "The Dmap path you entered: " + str(args.dmap[0]) + " does not appear to point to a Dmap file!!"
            dmap_good = False

        if args.save_path:
            if not os.path.exists(args.save_path):
                print "The save path directory you entered: " + str(args.save_path) + " does not exist!!"
                save_path_good = False
        else:
            if os.path.exists(os.path.dirname(args.dsets[0])):
                print "You did not enter a save directory path, program will default to using: " + str(os.path.dirname(args.dsets[0])) + " as the save path."
                args.save_path = os.path.dirname(args.dsets[0])
            else:
                print "Program tried to default to using: " + str(os.path.dirname(args.dsets[0])) + " as the save path but failed!!"
                save_path_good = False

        if dsets_good and parscan_good and dmap_good and save_path_good:
            placecor(gui = False, final_dset_list = args.dsets, dmap_path = args.dmap, parscan_path = args.parscan, save_dir = args.save_path)

if __name__ == "__main__":
    main()
