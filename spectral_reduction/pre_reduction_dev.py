#!/usr/bin/env python
from __future__ import print_function
import sys, os, pdb, argparse, datetime
from optparse import OptionParser
import util
import quick_reduc
import time
import glob
import matplotlib
import instruments
from astropy.io import fits, ascii
from pyraf import iraf

matplotlib.use('TkAgg')
from flat_utils import combine_flats,inspect_flat

def pre_reduction_dev():


    return 0


def pre_reduction_orig():

    description = "> Performs pre-reduction steps"
    usage = "%prog    \t [option] \n Recommended syntax: %prog -i -c"
  
    parser = OptionParser(usage=usage, description=description, version="0.1" )
    parser.add_option('-n','--nflats',type=int,default=50,
                        help='max number of flats to use (to prevent IRAF crash)')
    option, args = parser.parse_args()
    MAX_N_FLATS = option.nflats
    
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.specred(_doprint=0)

    iraf.ccdred.verbose = 'no'
    iraf.specred.verbose = 'no'
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''

    iraf.longslit.mode = 'h'
    iraf.specred.mode = 'h'
    iraf.noao.mode = 'h'
    iraf.ccdred.instrument = "ccddb$kpno/camera.dat"

    mkarc = raw_input("Make arc? ([y]/n): ")
    mkflat = raw_input("Make flat? ([y]/n): ")

    if len(args) > 1:
        files=[]
        sys.argv.append('--help')
        option, args = parser.parse_args()
        sys.exit()
    elif len(args) == 1:
        files = util.readlist(args[0])
        sys.exit()
    else:
        listfile = glob.glob('*.fits')
        files_science = []
        files_arc = []
        files_dflat = []
        #print 'checking your files ...'
        for img in listfile:
            _type = ''
            hdr0 = util.readhdr(img)
            _type=util.readkey3(hdr0, 'object')
            if 'flat' in _type.lower():
                files_dflat.append(img)
            elif 'arc' not in _type.lower() and 'arc' not in img.lower():
                files_science.append(img)
        if mkarc != 'n':
            mkarc_b = raw_input("List blue arc files to combine (.fits will be added): ").split()
            mkarc_r = raw_input("List red arc files to combine (.fits will be added): ").split()
            for arc in mkarc_b:
                files_arc.append(arc + '.fits')
            for arc in mkarc_r:
                files_arc.append(arc + '.fits')

    if mkarc != 'n':
        list_arc_b = []
        list_arc_r = []
        for arcs in files_arc:
            if instruments.blue_or_red(arcs)[0] == 'blue':
                list_arc_b.append(arcs)
            elif instruments.blue_or_red(arcs)[0] == 'red':
                list_arc_r.append(arcs)
            else:
                sys.exit()

    if mkflat != 'n':
        list_flat_b = []
        list_flat_r = []
        for dflats in files_dflat:
            if instruments.blue_or_red(dflats)[0] == 'blue':
                if len(list_flat_b) < MAX_N_FLATS:
                    list_flat_b.append(dflats)
            elif instruments.blue_or_red(dflats)[0] == 'red':
                if len(list_flat_r) < MAX_N_FLATS:
                    list_flat_r.append(dflats)
            else:
                sys.exit()
                
                
    # make pre_reduced if it doesn't exist
    if not os.path.isdir('pre_reduced/'):
        os.mkdir('pre_reduced/')
        
    # log the existing processed files (need to verify this works if pre_reduced is empty...)
    pfiles = []
    new_files = []
    for root, dirnames, filenames in os.walk('pre_reduced'):
        for file in filenames:
            if file.startswith('to'):
                pfiles.append(file)
    #print(pfiles)

    # loop over each image in pre_reduced
    for img in listfile:
        hdr = util.readhdr(img)
        targ=util.readkey3(hdr, 'object')
        
        # if file is not not a processed file, run the overscan+trim code
        if 'to'+ img not in pfiles:
            
            # if the file is a science file, grab the name for later
            if 'arc' not in targ.lower() and 'flat' not in targ.lower():
                new_files.append(img)
                print ('Adding data for: ' + targ)
                
            inst = instruments.blue_or_red(img)[1]

            iraf.specred.dispaxi = inst.get('dispaxis')
            iraf.longslit.dispaxi = inst.get('dispaxis')

            _biassec0 = inst.get('biassec')
            _trimsec0 = inst.get('trimsec')
            
            ######################################################################
            #
            # JB: this chunk of code needs attention
            # It seems hacky for anything but Kast...
            #
            # overscan
            if not img.startswith('o') and inst.get('observatory')=='lick':
                if os.path.isfile('pre_reduced/o'+img):
                    os.remove('pre_reduced/o'+img)
                util.kastbias(img,'pre_reduced/o'+img)
            elif not img.startswith('o') and inst.get('observatory')!='lick':
                if os.path.isfile('pre_reduced/o'+img):
                    os.remove('pre_reduced/o'+img)
                os.system('cp ' +  img + ' ' + 'pre_reduced/' + img)

                
            # trim
            if not img.startswith('t')and inst.get('observatory')=='lick':
                if os.path.isfile('pre_reduced/to'+img):
                    os.remove('pre_reduced/to'+img)
                iraf.ccdproc('pre_reduced/o'+img, output='pre_reduced/to'+img, 
                             overscan='no', trim='yes', zerocor="no", flatcor="no", 
                             readaxi='line',trimsec=str(_trimsec0), Stdout=1)

            elif not img.startswith('t')and inst.get('observatory')!='lick':
                if os.path.isfile('pre_reduced/to'+img):
                    os.remove('pre_reduced/to'+img)
                # iraf.ccdproc('pre_reduced/'+img, output='pre_reduced/to'+img, 
                #              overscan='yes', trim='yes', zerocor="no", flatcor="no", 
                #              readaxi='line',trimsec=str(_trimsec0), biassec=str(_biassec0), Stdout=1)
                iraf.ccdproc('pre_reduced/'+img, output='pre_reduced/to'+img, 
                             overscan='no', trim='yes', zerocor="no", flatcor="no", 
                             readaxi='line',trimsec=str(_trimsec0), biassec=str(_biassec0), Stdout=1)

        else:
            pass 
            # 'to'+img exists, so do nothing, 

    # combine the arcs
    if mkarc != 'n':
        
        # blue arcs
        if len(list_arc_b) > 0:
            if len(list_arc_b) == 1:
                arc_blue = list_arc_b[0]
                os.system('cp ' + 'pre_reduced/to'+ arc_blue + ' ' + 'pre_reduced/ARC_blue.fits')
            else:
                arc_str = ''
                for arc in list_arc_b:
                    arc_str = arc_str + 'pre_reduced/to'+ arc + ','
                if os.path.isfile('pre_reduced/ARC_blue.fits'):
                    os.remove('pre_reduced/ARC_blue.fits')
                iraf.imcombine(arc_str, output='pre_reduced/ARC_blue.fits')

        # red arcs
        if len(list_arc_r) > 0:
            if len(list_arc_r) == 1:
                arc_red = list_arc_r[0]
                os.system('cp ' + 'pre_reduced/to'+ arc_red + ' ' + 'pre_reduced/ARC_red.fits')
            else:
                arc_str = ''
                for arc in list_arc_r:
                    arc_str = arc_str + 'pre_reduced/to'+ arc + ','
                if os.path.isfile('pre_reduced/ARC_red.fits'):
                    os.remove('pre_reduced/ARC_red.fits')
                iraf.imcombine(arc_str, output='pre_reduced/ARC_red.fits')

    # combine the flats
    if mkflat != 'n':
        inter = 'yes'
        
        # blue flats
        if len(list_flat_b) > 0:
            br, inst = instruments.blue_or_red(list_flat_b[0])
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_blue = 'pre_reduced/toFlat_blue.fits'

            flat_list = []
            for flat in list_flat_b:
                flat_list.append('pre_reduced/to'+ flat)
            if os.path.isfile(Flat_blue):
                os.remove(Flat_blue)

            # first, combine all the flat files into a master flat
            res = combine_flats(flat_list,OUTFILE=Flat_blue,MEDIAN_COMBINE=True)
            
            # run iraf response
            iraf.specred.response(Flat_blue, 
                                   normaliz=Flat_blue, 
                                   response='pre_reduced/RESP_blue', 
                                   interac=inter, thresho='INDEF',
                                   sample='*', naverage=2, function='legendre', 
                                   low_rej=3,high_rej=3, order=60, niterat=20, 
                                   grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['pre_reduced/RESP_blue.fits'],DISPAXIS=dispaxis)

        # red flats
        if len(list_flat_r) > 0:
            br, inst = instruments.blue_or_red(list_flat_r[0])
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_red = 'pre_reduced/toFlat_red.fits'


            flat_list = []
            for flat in list_flat_r:
                flat_list.append('pre_reduced/to'+ flat)
            if os.path.isfile(Flat_red):
                os.remove(Flat_red)

            # first, combine all the flat files into a master flat
            res = combine_flats(flat_list,OUTFILE=Flat_red,MEDIAN_COMBINE=True)

            #What is the output here? Check for overwrite
            iraf.specred.response(Flat_red, 
                                  normaliz=Flat_red, 
                                  response='pre_reduced/RESP_red', 
                                  interac=inter, thresho='INDEF',
                                  sample='*', naverage=2, function='legendre', 
                                  low_rej=3,high_rej=3, order=80, niterat=20, 
                                  grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat(['pre_reduced/RESP_red.fits'],DISPAXIS=dispaxis)
    

    # science files should have 't' in front now
    # this just gets the base name, to prefix assumed below
    if new_files is not None:
        files_science = new_files

    # get all the science objects for the night
    science_targets = []
    for obj in files_science:
        hdr = util.readhdr(obj)
        _type=util.readkey3(hdr, 'object')
        science_targets.append(_type)

    # make a dir for each sci object
    science_targets = set(science_targets)
    for targ in science_targets:
        if not os.path.isdir('pre_reduced/' + targ + '/'):
            os.mkdir('pre_reduced/'+ targ + '/')

    # copy the files into the obj dir
    for obj in files_science:
        hdr = util.readhdr(obj)
        targ=util.readkey3(hdr, 'object')
        if not obj.startswith('to'):
            os.system('cp ' + 'pre_reduced/to'+ obj + ' ' + 'pre_reduced/' + targ + '/')
        else:
            os.system('cp ' +  'pre_reduced/'+ obj + ' ' + 'pre_reduced/' + targ + '/')

    rawfiles = glob.glob('*.fits')
    ofiles = glob.glob('pre_reduced/o'+ '*.fits')
    tfiles = glob.glob('pre_reduced/to'+ '*.fits')
    
    # delete raw files from the pre_reduced dir
    # there shouldn't be any there though?
    # maybe if the overscan isn't implemented for that detector
    for img in rawfiles:
        util.delete('pre_reduced/' + img)
        
    # delete the ofiles from pre_reduced dir
    for img in ofiles:
        util.delete(img)
    
    
    
def parse_cmd_args():
    ''' Parse the command line options '''

    # init parser
    descStr = 'Pre-reduction for the UCSC Spectroscopic Pipeline'
    parser = argparse.ArgumentParser(description=descStr)
    conflictGroup = parser.add_mutually_exclusive_group()

    # optional
    parser.add_argument('-v','--verbose',
                        help='print diagnostic info',action='store_true')

    parser.add_argument('-o','--original',
                        help='print diagnostic info',action='store_true')
    parser.add_argument('-c','--config-file',
                        help='config file to use for pre-reduction')

    # parse
    cmdArgs = parser.parse_args()

    #print cmdArgs to determine what defaults are
    pdb.set_trace()

    kwargs['VERBOSE'] = cmdArgs.verbose
    kwargs['ORIGINAL'] = cmdArgs.original
    kwargs['CONFIG_FILE'] = cmdArgs.config_file

    args = ()
    kwargs = {}

    return (args,kwargs)

def main(*args,**kwargs):
    '''
    Main driver for running pre-reduction versions

    Parameters
    ----------
    data : list
        List to parse for unique values

    IS_TUPLE : bool, optional
        If True, data is assumed to be a list of tuples

    KEY_INDEX : int, optional
        If IS_TUPLE, KEY_INDEX is used as the tuple index which
        holds the values to construct the unique set from

    Returns
    -------
    set : A set of unique values in data

    '''

    if kwargs.get('ORIGINAL'):
        res = pre_reduction_orig()

    else:
        if kwargs.get('CONFIG_FILE'):
            res = pre_reduction_dev(config=config)
        else:
            res = pre_reduction_dev(config=False)

    return 0

if __name__=='__main__':
    ''' Run parsing, then main '''
    args,kwargs = parse_cmd_args()
    main(*args,**kwargs)