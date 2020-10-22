#!/usr/bin/env python
from __future__ import print_function
import sys
import os
from optparse import OptionParser
import util
import quick_reduc
import time
import glob
import matplotlib
matplotlib.use('TkAgg')

# Returns
def get_files(directory, args):

    science = [] ; arc = [] ; flat = []
    if len(args)==1:
        files = util.readlist(args[0])
    else:
        file_pattern = os.path.join(directory, 't*.fits')
        arc_pattern = os.path.join(directory, '../ARC*.fits')
        resp_pattern = os.path.join(directory, '../RESP*.fits')
        files = glob.glob(file_pattern)
        files += glob.glob(arc_pattern)
        files += glob.glob(resp_pattern)

    for file in files:
        hdr0 = util.readhdr(file)
        _type=util.readkey3(hdr0, 'object')
        if _type.startswith('arc') or os.path.basename(file).startswith('ARC'):
            arc.append(file)
        elif 'RESP' in file:
            flat.append(file)
        else:
            science.append(file)

    files = []
    if len(arc)==0 and len(args)==1:
        q = '\n# No arc available.  Do you want to glob ../ARC*.fits ([y]/n)?'
        inp = raw_input(q)
        if q=='' or q=='y' or q=='yes':
            arc_pattern = os.path.join(directory, '../ARC*.fits')
            files = glob.glob(arc_pattern)
    if len(flat)==0 and len(args)==1:
        q = '\n# No arc available.  Do you want to glob ../RESP*.fits ([y]/n)?'
        inp = raw_input(q)
        if q=='' or q=='y' or q=='yes':
            flat_pattern = os.path.join(directory, '../RESP*.fits')
            files = glob.glob(flat_pattern)

    for file in files:
        hdr0 = util.readhdr(file)
        _type=util.readkey3(hdr0, 'object')
        if _type.startswith('arc'):
            arc.append(file)
        elif 'RESP' in file:
            flat.append(file)
        else:
            science.append(file)

    return(science, arc, flat)

if __name__ == '__main__':
    description = '> Fast reduction of spectra '
    usage = '%prog    \t [option] \n Recommended syntax: %prog -i -c'
    parser = OptionParser(usage=usage, description=description, version='0.1' )
    parser.add_option('-a', '--arc', dest='arc', action='store_true',
        help='identify arcs of the same night')
    parser.add_option('-i', '--interactive_extraction', dest='interactive',
        action='store_true',help='extract spectrum interactively')
    parser.add_option('-c', '--cosmic', dest='cosmic',
        action='store_true',help='cosmic ray removal')
    parser.add_option('-f', '--fast', dest='fast',
        action='store_true',help='fast reduction')
    parser.add_option('-b', '--back', '--background', dest='back',
        action='store_true',help='separately fit sky background')
    parser.add_option('-s','--separate','--sep', dest='separate',
        action='store_true',help='reduce 2D spectra separately then stack 1D')
    parser.add_option('--host', dest='host',
        action='store_true',help='host reduction')
    parser.add_option('--dir','-d', dest='directory', default='',
        help='directory to reduce if not current directory')

    opt, args = parser.parse_args()

    kwargs = {}
    # Pass options as argument to quick_reduc
    for key in vars(opt).keys():
        kwargs['_'+key] = vars(opt)[key]

    science, arc, flat = get_files(opt.directory, args)

    starttime = time.time()

    if len(science) > 0:
        print('\n#####################################\n### start of reduction')
        outputfile = quick_reduc.reduce(science, arc, flat, **kwargs)
        stoptime = time.time()
    else:
        print('\n### No science files! exiting...')

    stoptime = time.time()
    print('\n### wow, only ' + str(stoptime - starttime) + ' seconds')
    print('\n### end of reduction')
    print('\n#######################################')
