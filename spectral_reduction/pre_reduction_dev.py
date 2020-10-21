#!/usr/bin/env python
from __future__ import print_function
import sys, os, shutil, pdb, argparse, datetime, json
from optparse import OptionParser
import util
import quick_reduc
import time
import glob
import matplotlib
import instruments
from astropy.io import fits, ascii
from pyraf import iraf
import pyds9 as pyds9
import keck_basic_2d
import manifest_utils as mu
import host_galaxies as host_gals
import numpy as np

matplotlib.use('TkAgg')
from flat_utils import combine_flats,inspect_flat


def show_ds9_listfile(listFile,instanceName='default'):
    ''' reads list of images from a file, displays in a DS9 instance '''

    imgList = []

    # read in the filenames
    with open(listFile,'r') as fin:
        for line in fin:
            if line[0] != '#':
                fileName = line.split()[0]
                imgList.append(fileName)
    disp = show_ds9_list(imgList,instanceName=instanceName)
    return disp

def show_ds9_list(imgList,instanceName='default'):
    ''' display a list of images in a DS9 instance '''

    #Setup the DS9 instance
    disp = pyds9.DS9(instanceName)
    disp.set("frame delete all")
    disp.set("view image no")
    disp.set("view colorbar no")
    disp.set("view panner no")
    disp.set("view info yes")
    disp.set("view magnifier no")
    disp.set("view buttons no")
    disp.set("tile yes")
    disp.set("tile column")
    disp.set("width 1200")
    disp.set("height 275")

    #Display the images
    for i,img in enumerate(imgList):
        disp.set("frame {}".format(i))
        ds9cmd = "file fits {}".format(img)
        disp.set(ds9cmd)
        disp.set("scale mode minmax")
        disp.set("regions delete all")
        disp.set("scale log")
        disp.set("cmap invert yes")
        disp.set(ds9cmd)
    return disp


def config_to_pretty_str(configDict):
    ''' Dummy function for reformatting configDicts '''
    return json.dumps(configDict,indent=4)



def user_adjust_config(configDict,operation=''):
    ''' Allow user to make adjustments to their reduction config interactively '''

    if operation == 'REMOVE':

        usrResp = ''
        while usrResp.upper() not in ['C','Q']:
            configStr = config_to_pretty_str(configDict)

            promptStr = '\n\n\n You\'ve chosen to remove some file(s) from the config.\n'
            promptStr += 'Here\'s the current state:\n'
            promptStr += configStr
            promptStr += '\n\nAt this point you may:\n\n'
            promptStr += '  Enter the filename you wish to remove, or\n'
            promptStr += '  (C)ontinue with these files as is\nCommand: '
            usrRespOrig = raw_input(promptStr)


            try:
                usrResp = usrRespOrig.strip()
            except Exception as e:
                errStr = 'I had some problem with that input: {}'.format(e)
                print(errStr)

            if usrResp.upper() not in ['C','Q']:
                try:
                    removeFilename = usrResp
                    REMOVED = False
                    for imgType,typeDict in configDict.items():
                        for chan,objDict in typeDict.items():
                            for obj,fileList in objDict.items():
                                if removeFilename in fileList:
                                    configDict[imgType][chan][obj].remove(removeFilename)
                                    REMOVED = True
                                    print('i removed {}'.format(removeFilename))
                    if not REMOVED:
                        raise ValueError('Couldn\'t locate {} in config'.format(removeFilename))

                except (Exception,ValueError) as e:
                    errStr = 'I had some problem: {}'.format(str(e))
                    print(errStr)



    if operation == 'ADD':

        usrResp = ''
        while usrResp.upper() not in ['C','Q']:
            configStr = config_to_pretty_str(configDict)

            promptStr = '\n\n\n You\'ve chosen to add some file(s) from the config.\n'
            promptStr += 'Here\'s the current state:\n'
            promptStr += configStr
            promptStr += '\n\nAt this point you may:\n\n'
            promptStr += '  (C)ontinue with these files as is, or\n'
            promptStr += '  Enter the filename you wish to add, \n'
            promptStr += '  according to the format TYPE CHAN OBJECT FILENAME\n\n'
            promptStr += '  (e.g., CAL_ARC BLUE CALIBRATION_ARC r1021.fits)\n'
            promptStr += '  (e.g., CAL_FLAT RED CALIBRATION_FLAT r1022.fits)\n'
            promptStr += '  (e.g., STD BLUE BD28411 b1054.fits)\n'
            promptStr += '  (e.g., SCI BLUE 2019ltw b1055.fits)\n\n'
            promptStr += '  Make sure your input is accurately formatted; this program\n'
            promptStr += '  performs conservative validation and may not add your\n'
            promptStr += '  file if the input is not formatted correctly. In that case\n'
            promptStr += '  quit and just add the file(s) the config file by hand.\nCommand: '
            usrRespOrig = raw_input(promptStr)


            try:
                usrResp = usrRespOrig.strip()
            except Exception as e:
                errStr = 'I had some problem with that input: {}'.format(e)
                print(errStr)

            if usrResp.upper() not in ['C','Q']:
                try:
                    ADDED = False
                    inputList = usrResp.split()
                    imgType = inputList[0]
                    chan = inputList[1]
                    obj = inputList[2]
                    filename = inputList[3]

                    # check if file is available
                    availableFiles = [os.path.basename(file) for file in glob.glob('*.fits')]
                    if filename not in availableFiles:
                        raise ValueError('Couldn\'t find file {}'.format(filename))

                    # try to add
                    if imgType not in ['CAL_ARC','CAL_FLAT','STD','SCI']:
                        raise ValueError('Bad image type supplied!')
                    if chan not in ['BLUE','RED']:
                        raise ValueError('Bad channel supplied!')

                    if obj not in configDict[imgType][chan].keys():
                        configDict[imgType][chan][obj] = [filename]
                        ADDED = True
                    else:
                        if filename not in configDict[imgType][chan][obj]:
                            configDict[imgType][chan][obj].append(filename)
                            ADDED = True

                    if not ADDED:
                        raise ValueError('Couldn\'t add {} to config'.format(filename))

                except Exception as e:
                    errStr = 'I had some problem: {}'.format(str(e))
                    print(errStr)


    return configDict


def basic_2d_proc(rawFile,imgType=None,CLOBBER=False):

    # set up file names based on our convention
    oScanFile = 'pre_reduced/o{}'.format(rawFile)
    toScanFile = 'pre_reduced/to{}'.format(rawFile)
    toScanFlat = 'pre_reduced/to{}'.format(rawFile.replace('.fits','_norm.fits'))

    # run the basic 2D stuff on this image if necessary
    if (not os.path.isfile(oScanFile)) or CLOBBER:

        # get the instrument configuration
        inst = instruments.blue_or_red(rawFile)[1]
        iraf.specred.dispaxi = inst.get('dispaxis')
        iraf.longslit.dispaxi = inst.get('dispaxis')
        _biassec0 = inst.get('biassec')
        _trimsec0 = inst.get('trimsec')
        _flatsec0 = inst.get('flatsec')

        # remove destination files
        for file in [oScanFile,toScanFile]:
            if os.path.isfile(file) and CLOBBER:
                os.remove(file)

        # check the ultimate desination file, since intermediates get deleted
        if not os.path.isfile(toScanFile):
            if inst.get('instrument')=='kast':

                # do Lick specific bias operations
                util.kastbias(rawFile,oScanFile)
            elif inst.get('instrument')=='binospec':

                util.binospecbias(rawFile,oScanFile,inst)

            # do general (IRAF is in such a sorry state I'm not even sure if this works)
            else:
                iraf.ccdproc(rawFile, output=oScanFile,
                             overscan='yes', trim='yes',
                             zerocor="no", flatcor="no",readaxi='line',
                             trimsec=str(_trimsec0), biassec=str(_biassec0),
                             Stdout=1)

            # trim (same trimming operation for all telescopes)
            iraf.ccdproc(oScanFile, output=toScanFile,
                        overscan='no', trim='yes', zerocor="no", flatcor="no",
                        readaxi='line',trimsec=str(_trimsec0), Stdout=1)

            #create trimmed flats for norm region
            if imgType == 'CAL_FLAT' and inst.get('instrument')=='kast':
                iraf.ccdproc(oScanFile, output=toScanFlat,
                        overscan='no', trim='yes', zerocor="no", flatcor="no",
                        readaxi='line',trimsec=str(_flatsec0), Stdout=1)
            if imgType == 'CAL_FLAT' and inst.get('instrument')=='binospec':
                iraf.ccdproc(oScanFile, output=toScanFlat,
                        overscan='no', trim='yes', zerocor="no", flatcor="no",
                        readaxi='line',trimsec=str(_flatsec0), Stdout=1)
            os.remove(oScanFile)

    return 0


def fix_lris_header(rawFile):
    hdu = fits.open(rawFile)
    hdu[0].header['OBJECT'] = hdu[0].header['TARGNAME']
    hdu.writeto(rawFile, overwrite=True)

def reorg_files(configDict,CLOBBER=False):
    ''' Move files into their appropriate subdirectories '''

    for imgType,typeDict in configDict.items():
        if imgType in ['STD','SCI']:
            for chan,objDict in typeDict.items():
                for obj,fileList in objDict.items():

                    destDir = 'pre_reduced/{}'.format(obj)
                    if not os.path.isdir(destDir):
                        os.mkdir(destDir)

                    for rawFile in fileList:
                        procFile = 'pre_reduced/to{}'.format(rawFile)
                        destFile = '{}/{}'.format(destDir,os.path.basename(procFile))

                        if not os.path.isfile(destFile):
                            shutil.copy(procFile,destFile)

    return 0




def pre_reduction_dev(*args,**kwargs):

    # parse kwargs
    VERBOSE = kwargs.get('VERBOSE')
    CLOBBER = kwargs.get('CLOBBER')
    FAKE_BASIC_2D = kwargs.get('FAKE_BASIC_2D')
    FULL_CLEAN = kwargs.get('FULL_CLEAN')
    FAST = kwargs.get('FAST')
    CONFIG_FILE = kwargs.get('CONFIG_FILE')
    MAKE_ARCS = kwargs.get('MAKE_ARCS')
    MAKE_FLATS = kwargs.get('MAKE_FLATS')
    QUICK = kwargs.get('QUICK')
    RED_AMP_BAD = kwargs.get('RED_AMP_BAD')
    BLUE_AMP_BAD = kwargs.get('BLUE_AMP_BAD')
    HOST = kwargs.get('HOST')
    TRIM = kwargs.get('TRIM')
    FIX_AMP_OFFSET = kwargs.get('FIX_AMP_OFFSET')
    RED_DIR = kwargs.get('RED_DIR')

    # init iraf stuff
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

    prereddir= os.path.join(RED_DIR, 'pre_reduced/')

    # set up config
    if CONFIG_FILE:
        with open(CONFIG_FILE,'r') as fin:
            configDict = json.load(fin)
    else:
        STANDARD_STAR_LIBRARY = mu.construct_standard_star_library()
        observations = sorted(glob.glob('*.fits'))

        #TODO: Better first pass at config file, std have exptime < 250s(?)
        # CDK - updated this and placed it after arm/inst_name definition to
        # better generalize configDict
        configDict = {}
        for obsfile in observations:

            arm, inst_dict = instruments.blue_or_red(obsfile)
            inst_name = inst_dict.get('instrument').upper()
            for key in ['SCI','STD','CAL_ARC','CAL_FLAT','CAL']:
                if key not in configDict.keys():
                    configDict[key]={}
                if arm.upper() not in configDict[key].keys():
                    configDict[key][arm.upper()]={}
                    if 'CAL_' in key:
                        full_key = key.replace('CAL_','CALIBRATION_')
                        configDict[key][arm.upper()]={full_key: []}

            # CDK - added this to explicitly skip files that dont need to be
            # reduced.
            if not mu.needs_to_be_reduced(obsfile): continue

            hdu = fits.open(obsfile)
            if 'lris' in inst_dict.get('name'):
                hdu[0].header['OBJECT']=hdu[0].header['TARGNAME']
                hdu.writeto(obsfile, overwrite=True)

            use_ext = inst_dict['use_ext']
            header = hdu[use_ext].header

            imageType = mu.determine_image_type(header, inst_name, STANDARD_STAR_LIBRARY)

            channel, inst_dict = instruments.blue_or_red(obsfile)
            obj = header.get('OBJECT').strip()
            if imageType == 'SCI' or imageType == 'STD':
                if obj in configDict[imageType][channel.upper()].keys():

                    configDict[imageType][channel.upper()][obj].append(obsfile)
                else:
                    configDict[imageType][channel.upper()][obj] = [obsfile]
            if imageType == 'CAL_ARC' and 'foc' not in obsfile:
                configDict[imageType][channel.upper()]['CALIBRATION_ARC'].append(obsfile)
            if imageType == 'CAL_FLAT':
                configDict[imageType][channel.upper()]['CALIBRATION_FLAT'].append(obsfile)


        with open('custom_config.json','w') as fout:
            fout.write(json.dumps(configDict,indent=4))

        outStr = '\n\nOk, not config supplied, so I wrote a first pass custom_config.json\n'
        outStr += 'Use at your own risk! Manually edit if needed and run again with -c custom_config.json\n'
        outStr += 'You can manually add files to the appropriate lists and rerun.\n'
        outStr += 'WARNING: make sure you rename your config file, or it could get overwritten!\n\n'
        print(outStr)

        sys.exit(1)

    if not FAST:
        # do visual inspection of frames via ds9 windows
        usrResp = ''
        while usrResp != 'C':
            promptStr = '\nYou\'ve opted to display images before kicking off the reduction.\n'
            promptStr += 'At this point you may:\n'
            promptStr += '  (D)isplay the current state of the reduction config\n'
            promptStr += '  (C)ontinue with these files as is\n'
            promptStr += '  (R)emove a file from the current config\n'
            promptStr += '  (A)dd a file to the current config\n'
            promptStr += '  (Q)uit the whole thing. \n'
            promptStr += 'I recommend you (D)isplay and remove unwanted frames from your config file,\n'
            promptStr += '(Q)uit, and then rerun with the updated config file.\nCommand: '
            usrRespOrig = raw_input(promptStr)


            try:
                usrResp = usrRespOrig.strip().upper()
            except Exception as e:
                usrResp = 'nothing'

            # (D)isplay all the images in the lists
            if usrResp == 'D':

                blueArcList = configDict['CAL_ARC']['BLUE']['CALIBRATION_ARC']
                redArcList = configDict['CAL_ARC']['RED']['CALIBRATION_ARC']

                blueFlatList = configDict['CAL_FLAT']['BLUE']['CALIBRATION_FLAT']
                redFlatList = configDict['CAL_FLAT']['RED']['CALIBRATION_FLAT']

                blueStdList = []
                redStdList = []

                blueSciList = []
                redSciList = []

                for targ,imgList in configDict['STD']['BLUE'].items():
                    for file in imgList:
                        blueStdList.append(file)
                for targ,imgList in configDict['STD']['RED'].items():
                    for file in imgList:
                        redStdList.append(file)
                for targ,imgList in configDict['SCI']['BLUE'].items():
                    for file in imgList:
                        blueSciList.append(file)
                for targ,imgList in configDict['SCI']['RED'].items():
                    for file in imgList:
                        redSciList.append(file)

                blueArcDS9 = show_ds9_list(blueArcList,instanceName='BlueArcs')
                redArcDS9 = show_ds9_list(redArcList,instanceName='RedArcs')

                blueFlatDS9 = show_ds9_list(blueFlatList,instanceName='BlueFlats')
                redFlatDS9 = show_ds9_list(redFlatList,instanceName='RedFlats')

                blueStdDS9 = show_ds9_list(blueStdList,instanceName='BlueStandards')
                redStdDS9 = show_ds9_list(redStdList,instanceName='RedStandards')

                blueSciDS9 = show_ds9_list(blueSciList,instanceName='BlueScience')
                redSciDS9 = show_ds9_list(redSciList,instanceName='RedScience')

            if usrResp == 'R':
                configDict = user_adjust_config(configDict,operation='REMOVE')
            if usrResp == 'A':
                configDict = user_adjust_config(configDict,operation='ADD')
            if usrResp == 'Q':
                print('Okay, quitting pre_reduction...')
                sys.exit(1)


    # pre_reduced does not exist, needs to be made
    if not os.path.isdir(prereddir):
        os.mkdir(prereddir)

    if QUICK:
        file =glob.glob('*.fits')[0]
        inst = instruments.blue_or_red(file)[1]
        if 'kast' in inst['name']:
            b_inst = instruments.kast_blue
            r_inst = instruments.kast_red
        if 'lris' in inst['name']:
            b_inst = instruments.lris_blue
            r_inst = instruments.lris_red
        if 'goodman' in inst['name']:
            b_inst = instruments.goodman_blue
            r_inst = instruments.goodman_red
        if not os.path.isdir(prereddir+'master_files/'):
            os.mkdir(prereddir+'master_files/')
        b_arcsol = b_inst.get('archive_arc_extracted_id')
        b_resp = b_inst.get('archive_flat_file')
        r_arcsol = r_inst.get('archive_arc_extracted_id')
        r_resp = r_inst.get('archive_flat_file')
        if os.path.isdir(prereddir+'master_files/'):
            os.system('cp ' + b_arcsol + ' ' + 'pre_reduced/master_files/')
            os.system('cp ' + b_resp + ' ' + 'pre_reduced/')
            os.system('cp ' + r_arcsol + ' ' + 'pre_reduced/master_files/')
            os.system('cp ' + r_resp + ' ' + 'pre_reduced/')


    # pre_reduced exists, but we want to clobber/do a clean reduction
    elif FULL_CLEAN:

        promptStr = 'Do you really want to wipe pre_reduced? [y/n]: '
        usrRespOrig = raw_input(promptStr)
        if usrRespOrig and usrRespOrig[0].strip().upper() == 'Y':

            # remove all pre_reduced files
            shutil.rmtree(prereddir)
            os.mkdir(prereddir)

    # pre_reduced exists, need to document what is there
    else:

        # get existing pre_reduced files
        preRedFiles = glob.glob(prereddir+'*.fits')

    # loop over raw files in configDict, if the destination exists, do nothing
    # # otherwise, do the bias/reorient/trim/output/etc
    for imgType,typeDict in configDict.items():
        for chan,objDict in typeDict.items():
            for obj,fileList in objDict.items():
                for rawFile in fileList:
                    # try:
                    #     res = basic_2d_proc(rawFile,CLOBBER=CLOBBER)
                    #     if res != 0:
                    #         raise ValueError('Something bad happened in basic_2d_proc on {}'.format(rawFile))
                    # except (Exception,ValueError) as e:
                    #     print('Exception (basic_2d): {}'.format(e))

                    if not FAKE_BASIC_2D:
                        inst = instruments.blue_or_red(rawFile)[1]
                        if inst['name'] == 'lris_blue' or inst['name'] == 'lris_red':
                            fix_lris_header(rawFile)
                            # res = keck_basic_2d.main([rawFile])
                            if imgType != 'CAL_FLAT':
                                print (imgType)
                                res = keck_basic_2d.main([rawFile], TRIM=TRIM, ISDFLAT=False, RED_AMP_BAD=RED_AMP_BAD, MASK_MIDDLE_RED=False, MASK_MIDDLE_BLUE=False, FIX_AMP_OFFSET=FIX_AMP_OFFSET, BLUE_AMP_BAD=BLUE_AMP_BAD)
                            else:
                                print (imgType)
                                res = keck_basic_2d.main([rawFile], TRIM=TRIM, ISDFLAT = True, RED_AMP_BAD=RED_AMP_BAD, MASK_MIDDLE_RED=False, MASK_MIDDLE_BLUE=False, FIX_AMP_OFFSET=FIX_AMP_OFFSET, BLUE_AMP_BAD=BLUE_AMP_BAD)
                        else:
                            res = basic_2d_proc(rawFile,imgType=imgType,CLOBBER=CLOBBER)
                    else:
                        # here we're faking the basic 2D reduction because we've done
                        # specialized 2D reduction (e.g., keck_basic_2d)
                        res = 0
                    if res != 0:
                        raise ValueError('Something bad happened in basic_2d_proc on {}'.format(rawFile))

    # move the std and sci files into their appropriate directories
    try:
        res = reorg_files(configDict,CLOBBER=CLOBBER)
        if res != 0:
            raise ValueError('Something bad happened in reorg_files')
    except (Exception,ValueError) as e:
        print('Exception (reorg): {}'.format(e))


    ### some blocks of code from the original pre_reduction ###
    # combine the arcs
    if MAKE_ARCS:
        # CDK - generalized arc creation method
        for key in configDict['CAL_ARC'].keys():
            list_arc = configDict['CAL_ARC'][key]['CALIBRATION_ARC']
            if len(list_arc)>0:
                first = prereddir+'to{}'.format(list_arc[0])
                br, inst = instruments.blue_or_red(first)
                destFile = prereddir+'ARC_{0}.fits'.format(key.lower())

                util.make_arc(list_arc, destFile, inst, iraf)


    # combine the flats
    if MAKE_FLATS and 'lris' in inst['name']:

        list_flat_b = configDict['CAL_FLAT']['BLUE']['CALIBRATION_FLAT']
        list_flat_r = configDict['CAL_FLAT']['RED']['CALIBRATION_FLAT']
        inter = 'yes'

        b_amp1_list = []
        b_amp2_list = []
        r_amp1_list = []
        r_amp2_list = []
        predir = prereddir+'to'
        for flat in list_flat_b:
            suffix = '.'.join(flat.split('.')[1:])
            b_amp1_file = flat.split('.')[0]+'_amp1.'+suffix
            b_amp2_file = flat.split('.')[0]+'_amp2.'+suffix
            if os.path.exists(predir+b_amp1_file): b_amp1_list.append(b_amp1_file)
            if os.path.exists(predir+b_amp2_file): b_amp2_list.append(b_amp2_file)
        for flat in list_flat_r:
            suffix = '.'.join(flat.split('.')[1:])
            r_amp1_file = flat.split('.')[0]+'_amp1.'+suffix
            r_amp2_file = flat.split('.')[0]+'_amp2.'+suffix
            print(r_amp1_file)
            if os.path.exists(predir+r_amp1_file): r_amp1_list.append(r_amp1_file)
            if os.path.exists(predir+r_amp2_file): r_amp2_list.append(r_amp2_file)

        # blue flats
        if len(list_flat_b) > 0:
            # br, inst = instruments.blue_or_red(list_flat_b[0])
            br, inst = instruments.blue_or_red(prereddir+'to{}'.format(b_amp1_list[0]))
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_blue_amp1 = prereddir+'toFlat_blue_amp1.fits'
            Flat_blue_amp2 = prereddir+'toFlat_blue_amp2.fits'

            flat_list_amp1 = []
            for flat in b_amp1_list:
                flat_list_amp1.append(prereddir+'to'+ flat)
            if os.path.isfile(Flat_blue_amp1):
                os.remove(Flat_blue_amp1)

            # first, combine all the flat files into a master flat
            res = combine_flats(flat_list_amp1,OUTFILE=Flat_blue_amp1,MEDIAN_COMBINE=True)

            # run iraf response
            iraf.specred.response(Flat_blue_amp1,
                                   normaliz=Flat_blue_amp1,
                                   response=prereddir+'RESP_blue_amp1',
                                   interac=inter, thresho='INDEF',
                                   sample='*', naverage=2, function='legendre',
                                   low_rej=5,high_rej=5, order=60, niterat=20,
                                   grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            res = inspect_flat([prereddir+'RESP_blue_amp1.fits'],
                OUTFILE=prereddir+'RESP_blue_amp1.fits', DISPAXIS=dispaxis)

            hdu_amp1 = fits.open(prereddir+'RESP_blue_amp1.fits')
            amp1_flatten = np.asarray(hdu_amp1[0].data).flatten()
            shape1=hdu_amp1[0].data.shape
            dim1 = shape1[0]
            dim2 = shape1[1]

            concat_amps = amp1_flatten

            if not BLUE_AMP_BAD:
                flat_list_amp2 = []
                for flat in b_amp2_list:
                    flat_list_amp2.append(prereddir+'to'+ flat)
                if os.path.isfile(Flat_blue_amp2):
                    os.remove(Flat_blue_amp2)

                res = combine_flats(flat_list_amp2,OUTFILE=Flat_blue_amp2,MEDIAN_COMBINE=True)
                iraf.specred.response(Flat_blue_amp2,
                                       normaliz=Flat_blue_amp2,
                                       response=prereddir+'RESP_blue_amp2',
                                       interac=inter, thresho='INDEF',
                                       sample='*', naverage=2, function='legendre',
                                       low_rej=5,high_rej=5, order=60, niterat=20,
                                       grow=0, graphic='stdgraph')

                res = inspect_flat([prereddir+'RESP_blue_amp2.fits'],
                    OUTFILE=prereddir+'RESP_blue_amp2.fits', DISPAXIS=dispaxis)
                hdu_amp2 = fits.open(prereddir+'RESP_blue_amp2.fits')
                amp2_flatten = np.asarray(hdu_amp2[0].data).flatten()
                shape2=hdu_amp2[0].data.shape

                dim1+=shape2[0]
                concat_amps = np.concatenate([concat_amps, amp2_flatten])

            print('Output blue response dimensions:',dim1,dim2)
            resp_blue_data = np.reshape(concat_amps, (dim1, dim2))

            header = hdu_amp1[0].header
            if os.path.exists(prereddir+'RESP_blue.fits'):
                os.remove(prereddir+'RESP_blue.fits')

            hdu = fits.PrimaryHDU(resp_blue_data,header)
            hdu.writeto(prereddir+'RESP_blue.fits',output_verify='ignore')

            resp_files = ['RESP_blue_amp1.fits','RESP_blue_amp2.fits']
            for file in resp_files:
                if os.path.exists(prereddir+file):
                    os.remove(prereddir+file)

        # red flats
        if len(list_flat_r) > 0:
            # br, inst = instruments.blue_or_red(list_flat_r[0])
            br, inst = instruments.blue_or_red(prereddir+'to{}'.format(r_amp1_list[0]))
            dispaxis = inst.get('dispaxis')
            iraf.specred.dispaxi = dispaxis
            Flat_red_amp1 = prereddir+'toFlat_red_amp1.fits'
            Flat_red_amp2 = prereddir+'toFlat_red_amp2.fits'


            flat_list_amp1 = []
            for flat in r_amp1_list:
                flat_list_amp1.append(prereddir+'to'+ flat)
            if os.path.isfile(Flat_red_amp1):
                os.remove(Flat_red_amp1)

            flat_list_amp2 = []
            for flat in r_amp2_list:
                flat_list_amp2.append(prereddir+'to'+ flat)
            if os.path.isfile(Flat_red_amp2):
                os.remove(Flat_red_amp2)

            amp2flag = len(flat_list_amp2)>0

            # first, combine all the flat files into a master flat
            if amp2flag:
                res = combine_flats(flat_list_amp1,OUTFILE=Flat_red_amp1,MEDIAN_COMBINE=True)
                res = combine_flats(flat_list_amp2,OUTFILE=Flat_red_amp2,MEDIAN_COMBINE=True)
            else:
                res = combine_flats(flat_list_amp1,OUTFILE=Flat_red_amp1,MEDIAN_COMBINE=True)

            #What is the output here? Check for overwrite
            iraf.specred.response(Flat_red_amp1,
                                  normaliz=Flat_red_amp1,
                                  response='pre_reduced/RESP_red_amp1',
                                  interac=inter, thresho='INDEF',
                                  sample='*', naverage=2, function='legendre',
                                  low_rej=5,high_rej=5, order=80, niterat=20,
                                  grow=0, graphic='stdgraph')
            if amp2flag:
                iraf.specred.response(Flat_red_amp2,
                                      normaliz=Flat_red_amp2,
                                      response='pre_reduced/RESP_red_amp2',
                                      interac=inter, thresho='INDEF',
                                      sample='*', naverage=2, function='legendre',
                                      low_rej=5,high_rej=5, order=80, niterat=20,
                                      grow=0, graphic='stdgraph')

            # finally, inspect the flat and mask bad regions
            if amp2flag:
                res = inspect_flat([prereddir+'RESP_red_amp1.fits'],
                    OUTFILE=prereddir+'RESP_red_amp1.fits', DISPAXIS=dispaxis)
                res = inspect_flat([prereddir+'RESP_red_amp2.fits'],
                    OUTFILE=prereddir+'RESP_red_amp2.fits', DISPAXIS=dispaxis)
            else:
                res = inspect_flat([prereddir+'RESP_red_amp1.fits'],
                    OUTFILE=prereddir+'RESP_red.fits', DISPAXIS=dispaxis)

            if amp2flag:
                hdu_amp1 = fits.open(prereddir+'RESP_red_amp1.fits')
                hdu_amp2 = fits.open(prereddir+'RESP_red_amp2.fits')
                head = hdu_amp1[0].header
                shape1 = hdu_amp1[0].data.shape
                shape2 = hdu_amp2[0].data.shape
                amp1_flatten = np.asarray(hdu_amp1[0].data).flatten()
                amp2_flatten = np.asarray(hdu_amp2[0].data).flatten()
                concat_amps = np.concatenate([amp2_flatten, amp1_flatten])
                xbin, ybin = [int(ibin) for ibin in head['BINNING'].split(',')]
                resp_red_data = np.reshape(concat_amps, (shape1[0]+shape2[0],shape1[1]))

                resp_red_data[278:294,:] = 1.

                header = hdu_amp1[0].header
                if os.path.isfile(prereddir+'RESP_red.fits'):
                    os.remove(prereddir+'RESP_red.fits')

                hdu = fits.PrimaryHDU(resp_red_data,header)
                hdu.writeto(prereddir+'RESP_red.fits',output_verify='ignore')

                os.remove(prereddir+'RESP_red_amp1.fits')
                os.remove(prereddir+'RESP_red_amp2.fits')
            else:
                os.remove(prereddir+'RESP_red_amp1.fits')

    elif MAKE_FLATS:
        # CDK - generalized arc creation method
        for key in configDict['CAL_FLAT'].keys():
            list_flat = configDict['CAL_FLAT'][key]['CALIBRATION_FLAT']
            if len(list_flat)>0:

                first = prereddir+'to{}'.format(list_flat[0])
                br, inst = instruments.blue_or_red(first)
                dispaxis = inst.get('dispaxis')

                iraf.specred.dispaxi = dispaxis
                inter = True
                Flat_out = prereddir+'toFlat_{0}.fits'.format(br.lower())
                dum = prereddir+'dummy_{0}.fits'.format(br.lower())
                resp = prereddir+'RESP_{0}.fits'.format(br.lower())

                flat_list = []
                norm_list = []
                for flat in list_flat:
                    flat_list.append(prereddir+'to'+flat)
                    norm = flat.replace('.fits','_norm.fits')
                    norm_list.append(prereddir+'to'+norm)
                for file in [Flat_out, dum, resp]:
                    if os.path.exists(file):
                        os.remove(Flat_out)

                # first, combine all the flat files into a master flat
                res = combine_flats(flat_list,OUTFILE=Flat_out,
                    MEDIAN_COMBINE=True)
                # combine all the flat files for norm region
                res = combine_flats(norm_list,OUTFILE=dum, MEDIAN_COMBINE=True)

                iraf.specred.response(Flat_out, normaliz=dum, response=resp,
                    interac=inter, thresho='INDEF', sample='*', naverage=2,
                    function='legendre', low_rej=5,high_rej=5, order=60,
                    niterat=20, grow=0, graphic='stdgraph')

                # finally, inspect the flat and mask bad regions
                res = inspect_flat([resp], OUTFILE=resp, DISPAXIS=dispaxis)

                for flat in norm_list:
                    os.remove(flat)

    if HOST:
        host_gals.make_host_metadata(configDict)

    return 0

def parse_cmd_args():
    ''' Parse the command line options '''

    # init parser
    descStr = 'Pre-reduction for the UCSC Spectroscopic Pipeline'
    parser = argparse.ArgumentParser(description=descStr)
    configCG = parser.add_mutually_exclusive_group()
    basicProcCG = parser.add_mutually_exclusive_group()

    # optional
    parser.add_argument('-q','--quicklook',
                        help='Move archival calibrations to enable quicker reductions', action='store_true')
    parser.add_argument('-v','--verbose',
                        help='Print diagnostic info',action='store_true')
    parser.add_argument('-k','--clobber',
                        help='Clobber existing files/dirs',action='store_true')
    parser.add_argument('-f','--fast',
                        help='Use config as is, don\'t prompt user for anything',action='store_true')


    parser.add_argument('--make-arcs',
                        help='Combine arcs into master arc images', action='store_true')
    parser.add_argument('--make-flats',
                        help='Combine flats into master flat images', action='store_true')
    parser.add_argument('--host',
                        help='Obtain relevant host galaxy metadata', action='store_true')
    parser.add_argument('--red-amp-bad',
                        help='Red side amplifier is bad so trim', action='store_true')
    parser.add_argument('--blue-amp-bad',
                        help='Red side amplifier is bad so trim', action='store_true')
    parser.add_argument('--fix-amp-offset',
                        help='Fix amplifier offset', action='store_true')
    parser.add_argument('--no-trim','--nt',
        help='Do not trim images during 2D processing',action='store_true')
    parser.add_argument('--red-dir', default='',
        help='Reduction directory if different from current')

    basicProcCG.add_argument('--fake-basic-2d',
                        help='Fake the basic 2D reductions', action='store_true')
    basicProcCG.add_argument('--full-clean',
                        help='Completely wipe pre_reduction (USE WITH CAUTION)', action='store_true')

    configCG.add_argument('-o','--original',
                    help='Run the original pre_reduction',action='store_true')
    configCG.add_argument('-c','--config-file',
                    help='Config file to use for pre-reduction')




    # parse
    cmdArgs = parser.parse_args()


    args = ()
    kwargs = {}

    kwargs['VERBOSE'] = cmdArgs.verbose
    kwargs['CLOBBER'] = cmdArgs.clobber
    kwargs['FULL_CLEAN'] = cmdArgs.full_clean
    kwargs['FAKE_BASIC_2D'] = cmdArgs.fake_basic_2d
    kwargs['FAST'] = cmdArgs.fast
    kwargs['ORIGINAL'] = cmdArgs.original
    kwargs['CONFIG_FILE'] = cmdArgs.config_file
    kwargs['MAKE_ARCS'] = cmdArgs.make_arcs
    kwargs['MAKE_FLATS'] = cmdArgs.make_flats
    kwargs['QUICK'] = cmdArgs.quicklook
    kwargs['HOST'] = cmdArgs.host
    kwargs['RED_AMP_BAD'] = cmdArgs.red_amp_bad
    kwargs['BLUE_AMP_BAD'] = cmdArgs.blue_amp_bad
    kwargs['TRIM'] = not cmdArgs.no_trim
    kwargs['FIX_AMP_OFFSET'] = cmdArgs.fix_amp_offset
    kwargs['RED_DIR'] = cmdArgs.red_dir

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
        res = pre_reduction_dev(**kwargs)

    return 0

if __name__=='__main__':
    ''' Run parsing, then main '''
    args,kwargs = parse_cmd_args()
    main(*args,**kwargs)
