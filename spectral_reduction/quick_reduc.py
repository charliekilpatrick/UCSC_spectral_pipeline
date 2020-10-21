from __future__ import print_function

def reduce(imglist, files_arc, files_flat, _cosmic=False,
    _interactive=False, _arc=False, _fast=False, _host=False,
    _back=False, **kwargs):

    import string
    import os
    import re
    import sys
    import pdb
    import shutil
    os.environ["PYRAF_BETA_STATUS"] = "1"
    try:      from astropy.io import fits
    except:      import   pyfits as fits
    import numpy as np
    import glob
    import util
    import instruments
    import combine_sides as cs
    import cosmics
    from pyraf import iraf
    import pyzapspec
    import host_galaxies as host_gals

    dv = util.dvex()
    scal = np.pi / 180.

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.specred(_doprint=0)
    iraf.disp(inlist='1', reference='1')

    toforget = ['ccdproc', 'imcopy',
                'specred.apall', 'longslit.identify',
                'longslit.reidentify', 'specred.standard',
                'longslit.fitcoords', 'onedspec.wspectext']
    for t in toforget:
        iraf.unlearn(t)
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


    list_arc = {}
    for arcs in files_arc:
        hdr = util.readhdr(arcs)
        arm, inst = instruments.blue_or_red(arcs)
        if arm not in list_arc.keys():
            list_arc[arm]=[arcs]
        else:
            list_arc[arm].append(arcs)

    asci_files = []
    newlist = {}

    print('\n### images to reduce :',imglist)
    for img in imglist:
        arm, inst = instruments.blue_or_red(img)
        if arm not in newlist.keys():
            newlist[arm] = [img]
        else:
            newlist[arm].append(img)

    if len(newlist.keys())==0:
        print('ERROR: no input files')
        sys.exit()
    elif len(newlist.keys())==1:
        key = newlist.keys()[0]
        print('Only one side, reducing side={0}'.format(key))
        newlist = [newlist[key] for key in newlist.keys()]
    else:
        query_str = "Reduce which side? ([all]"
        q = '/'+'/'.join([str(k)[0] for k in newlist.keys()])
        query_str += q + '):'
        sides = raw_input(query_str)
        if not sides: sides='all'

        key_map = {'b':'blue','blue':'blue','r':'red','red':'red','o':'only',
            'only':'only'}
        if sides=='all' or sides=='a':
            newlist = [newlist[key] for key in newlist.keys()]
        elif sides in key_map.keys():
            key = key_map[sides]
            if key not in newlist.keys():
                print('ERROR: side={0} not in file list'.format(key))
                sys.exit()
            else:
                newlist = [newlist[key]]

    for imgs in newlist:
        print (imgs)
        hdr = util.readhdr(imgs[0])
        br, inst = instruments.blue_or_red(imgs[0])
        flat_file = '../RESP_{0}'.format(br.lower())

        iraf.specred.dispaxi = inst.get('dispaxis')
        iraf.longslit.dispaxi = inst.get('dispaxis')

        _gain = inst.get('gain')
        _ron = inst.get('read_noise')
        iraf.specred.apall.readnoi = _ron
        iraf.specred.apall.gain = _gain

        _object0 = util.readkey3(hdr, 'OBJECT')
        _date0 = util.readkey3(hdr, 'DATE-OBS')

        _object0 = re.sub(' ', '', _object0)
        _object0 = re.sub('/', '_', _object0)
        nameout0 = str(_object0) + '_' + inst.get('name')

        exdir = _object0+'_ex/'

        nameout0 = util.name_duplicate(imgs[0], nameout0, '')
        print ('NAMEOUT:', nameout0)
        timg = nameout0
        print('\n### now processing :',timg,' for -> ',inst.get('name'))

        print ('IMAGES: ', imgs)
        if len(imgs) > 1:
            img_str = ''
            tmphdu = None
            for kk,i in enumerate(imgs):
                if _cosmic:
                    print('\n### starting cosmic removal')
                    files = glob.glob('*.fits')
                    if 'cosmic_{}'.format(i) not in glob.glob('*.fits'):

                        outimg,outmask,header = pyzapspec.pyzapspec(i,
                            outfile='cosmic_{}'.format(i),
                            WRITE_OUTFILE = True,
                            boxsize=inst.get('pyzap_boxsize',7),
                            nsigma=inst.get('pyzap_nsigma',16),
                            subsigma=inst.get('pyzap_subsigma',3))
                    img = 'cosmic_{}'.format(i)
                    img_str = img_str + img + ','

                    print('\n### cosmic removal finished')
                else:
                    print('\n### No cosmic removal, saving normalized image for inspection???')

                    img_str = img_str + i + ','
            print (img_str)
            combine_type = 'average'
            scale_type = None
            if len(imgs)>2:
                combine_type='median'
                scale_type = '!EXPTIME'
            elif len(imgs)==2:
                combine_type='average'
                scale_type='!EXPTIME'
            # CDK - keep getting a bunch of weird downstream effects when timg
            # already exists, so adding this
            if os.path.exists(timg):
                os.remove(timg)
            iraf.imcombine(img_str, output=timg, combine=combine_type,
                scale=scale_type)
        else:
            i = imgs[0]
            if _cosmic:
                print('\n### starting cosmic removal')
                files = glob.glob('*.fits')
                if 'cosmic_{}'.format(i) not in glob.glob('*.fits'):

                    outimg,outmask,header = pyzapspec.pyzapspec(i,
                        outfile='cosmic_{}'.format(i),
                        WRITE_OUTFILE = True,
                        boxsize=inst.get('pyzap_boxsize',7),
                        nsigma=inst.get('pyzap_nsigma',16),
                        subsigma=inst.get('pyzap_subsigma',3))
                img = 'cosmic_{}'.format(i)

                print('\n### cosmic removal finished')

            else:
                img = i #TH: Needs to redefine img other the script will make a imcopy of the last item in imglist.
                print('\n### No cosmic removal, saving normalized image for inspection???')

            if os.path.isfile(timg):
                os.system('rm -rf ' + timg)
            iraf.imcopy(img, output=timg)

        iraf.ccdproc(timg, output='',
                           overscan='no',
                           trim='no',
                           zerocor="no",
                           flatcor="yes",
                           readaxi='line',
                           flat=flat_file,
                           Stdout=1)

        if 'kast_red' in inst.get('name'):
            tfits = fits.open(timg)
            tdata = tfits[0].data
            theader = tfits[0].header
            if theader.get('Mirrored', None) is None:
                print ('Mirroring Image!')
                flip_tdata = np.fliplr(tdata)
                theader.set('Mirrored',  'True')
                hdu = fits.PrimaryHDU(flip_tdata, theader)
                hdu.writeto(timg,output_verify='ignore', clobber=True)

        img = timg
        arcfile = list_arc[br][0]

        if arcfile is not None and not arcfile.endswith('.fits'):
            arcfile=arcfile+'.fits'

        if not os.path.isdir('database/'):
            os.mkdir('database/')

        #There is a bug in identify when the path to the coordlist is too long
        #This is my hacky workaround for now, the file  is deleted later
        os.system('cp ' + inst.get('line_list') + ' .')

        if _arc:
            if not os.path.isdir('../master_files/'):
                os.mkdir('../master_files/')

            arcref = None
            arcfile = 'ARC_{0}.fits'.format(br.lower())
            wave_sol_file = 'idARC_{0}.ms'.format(br.lower())

            masters = [os.path.basename(x) for x in glob.glob('../master_files/*')]
            if wave_sol_file in masters:
                wave_sol= raw_input("Use your master wavelength solution? [y]/n: ") or 'y'
                if wave_sol.strip().lower() == 'y':
                    print ('Copying master file')
                    arc_ex=re.sub('.fits', '.ms.fits', arcfile)
                    os.system('cp ' + '../master_files/' + wave_sol_file + ' ./database/')
                else:
                    os.system('cp ' + '../' + arcfile + ' .')
                    arc_ex=re.sub('.fits', '.ms.fits', arcfile)
                    print('\n### arcfile : ',arcfile)
                    print('\n### arcfile extraction : ',arc_ex)
                    print(inst.get('line_list'))
                    print ('here')
                    # remove the file from the extract destination
                    if os.path.isfile(arc_ex):
                        os.remove(arc_ex)

                    iraf.specred.apall(arcfile,
                                        output=arc_ex,
                                        line = inst.get('approx_extract_line'),
                                        nsum=10,
                                        interactive='no',
                                        extract='yes',
                                        find='yes',
                                        nfind=1 ,
                                        format='multispec',
                                        trace='no',
                                        back='no',
                                        recen='no')
                    iraf.longslit.identify(images=arc_ex,
                                            section='line {}'.format(inst.get('approx_extract_line')),
                                            coordli='lines.dat',
                                            function = 'spline3',
                                            order=3,
                                            mode='h')
            else:
                os.system('cp ' + '../' + arcfile + ' .')
                arc_ex=re.sub('.fits', '.ms.fits', arcfile)
                print('\n### arcfile : ',arcfile)
                print('\n### arcfile extraction : ',arc_ex)
                print(inst.get('line_list'))
                # remove the file from the extract destination
                if os.path.isfile(arc_ex):
                    os.remove(arc_ex)

                iraf.specred.apall(arcfile,
                                    output=arc_ex,
                                    line = inst.get('approx_extract_line'),
                                    nsum=10,
                                    interactive='no',
                                    extract='yes',
                                    find='yes',
                                    nfind=1 ,
                                    format='multispec',
                                    trace='no',
                                    back='no',
                                    recen='no')
                iraf.longslit.identify(images=arc_ex,
                    section='line {}'.format(inst.get('approx_extract_line')),
                    coordli='lines.dat', function='spline3', order=3, mode='h')

                os.system('cp ' + 'database/' + wave_sol_file + ' ../master_files/')

        else:
            arcref = None
            arcfile = 'ARC_{0}.fits'.format(br.lower())
            wave_sol_file = 'idARC_{0}.ms'.format(br.lower())

            os.system('cp ' + '../' + arcfile + ' .')
            arc_ex=re.sub('.fits', '.ms.fits', arcfile)

            arcref = inst.get('archive_arc_extracted')
            arcref_img = string.split(arcref, '/')[-1]
            arcref_img = arcref_img.replace('.ms.fits', '')
            arcrefid = inst.get('archive_arc_extracted_id')
            os.system('cp ' + arcref + ' .')
            arcref = string.split(arcref, '/')[-1]
            os.system('cp ' + arcrefid + ' ./database')

            aperture = inst.get('archive_arc_aperture')
            os.system('cp ' + aperture + ' ./database')

            print('\n###  arcfile : ',arcfile)
            print('\n###  arcfile extraction : ',arc_ex)
            print('\n###  arc reference : ',arcref)

            # read for some meta data to get the row right
            tmpHDU = pyfits.open(arcfile)
            header = tmpHDU[0].header
            try:
                spatialBin = int(header['binning'].split(',')[0])
            except KeyError:
                spatialBin = 1
            apLine = 700//spatialBin

            iraf.specred.apall(arcfile,
                               output=arc_ex,
                               repef=arcref_img,
                               line = apLine,
                               nsum=10,
                               interactive='no',
                               extract='yes',
                               find='yes',
                               nfind=1 ,
                               format='multispec',
                               trace='no',
                               back='no',
                               recen='no')

            iraf.longslit.reidentify(referenc=arcref,
                                     images=arc_ex,
                                     interac='NO',
                                     section='middle line',
                                     coordli='lines.dat',
                                     shift='INDEF',
                                     search='INDEF',
                                     mode='h',
                                     verbose='YES',
                                     step=1,
                                     nsum=10,
                                     nlost=2,
                                     cradius=5,
                                     refit='yes',
                                     overrid='yes',
                                     newaps='no')

            os.system('cp ' + 'database/' + wave_sol_file + ' ../master_files/')

        util.delete('lines.dat')

        skyimg = None
        if _back:
            subimg, skyimg = util.subtract_background(img, inst,
                _inter=_interactive)
            img = subimg


        print('\n### extraction using apall')
        result = []
        hdr_image = util.readhdr(img)
        _type=util.readkey3(hdr_image, 'object')

        if (_type.startswith("arc") or
            _type.startswith("dflat") or
            _type.startswith("Dflat") or
            _type.startswith("Dbias") or
            _type.startswith("Bias")):
            print('\n### warning problem \n exit ')
            sys.exit()

        match_aperture = raw_input('Match aperture? y/[n]: ') or 'n'
        match_flag = match_aperture!='n'
        if _host:
            objname = _object0.lower().split('_')[0]
            host_dat = host_gals.calculate_ap_data(objname, inst)
            ap_pixs_phys, ap_pixs_sky, ap_width_kron, sep_pix,\
                ap_widths_arcsec, ap_widths_kpc, r_kron_rad = host_dat
            host_gals.write_host_ap(ap_pixs_phys, ap_pixs_sky, ap_width_kron,
                sep_pix, nameout0.split('.')[0], ap_widths_arcsec,
                ap_widths_kpc, r_kron_rad)

        imgex = util.extractspectrum(img, dv, inst, _interactive,
            'obj', host_ex=_host, match_flag=match_flag, skyimg=skyimg)

        save_ap = raw_input('Save as a master aperture ? y/[n]: ')
        outap = 'ap'+img.replace('.fits','')
        if save_ap!='n':
            shutil.copyfile('database/'+outap, '../master_files/'+outap)

        if os.path.exists('database/'+outap):
            with open('database/'+outap) as apfile:
                ap_data = apfile.readlines()
                count = 0
                for line in ap_data:
                    if line.startswith('begin'):
                        count+=1
                if count > 1:
                    with open('database/' + wave_sol_file) as arc:
                        with open('database/' + wave_sol_file + '_new','w') as arcex_new:
                            lines = arc.readlines()
                            for i in range(count):
                                for line in lines:
                                    if 'identify' in line:
                                        arcex_new.write(line.replace('Ap 1', 'Ap '+ str(i+1)))
                                    elif 'image' in line:
                                        arcex_new.write(line.replace('Ap 1', 'Ap '+ str(i+1)))
                                    elif 'aperture' in line:
                                        arcex_new.write(line.replace('1', str(i+1)))
                                    else:
                                        arcex_new.write(line)
                                arcex_new.write('\n')

                    shutil.copyfile('database/'+wave_sol_file+'_new',
                        'database/'+wave_sol_file)
                    util.delete('database/'+wave_sol_file+'_new')

        print('\n### applying wavelength solution:',arc_ex)
        iraf.disp(inlist=imgex, reference=arc_ex)

        result = result + [imgex] + [timg]

        if not os.path.isdir(exdir):
            os.mkdir(exdir)

        if not _arc:
            util.delete(arcref)
        else:
            util.delete(arcfile)

        util.delete(arc_ex)
        util.delete(imgex)
        if arcref is not None:
            util.delete(arcref)
        util.delete('logfile')

        shutil.move('d'+imgex, exdir+'d'+imgex)

        use_master = raw_input('Use master flux calibration? [y]/n ')
        if use_master != 'n':
            fluxstar='fluxstar'+inst.get('arm')+'.fits'
            bstar='bstar' + inst.get('arm') +  '.fits'
            shutil.copyfile('../master_files/'+fluxstar,exdir+fluxstar)
            shutil.copyfile('../master_files/'+bstar,exdir+bstar)
        else:
            use_sens = raw_input('Use archival flux calibration? [y]/n ')
            if use_sens != 'n':
                sensfile = inst.get('archive_sens')
                sensbase = os.path.split(sensfile)[1]
                bstarfile = inst.get('archive_bstar')
                bstarbase = os.path.split(bstarbase)

                if os.path.exists(sensfile):
                    shutil.copyfile(sensfile, exdir+sensbase)
                else:
                    print('WARNING: no sensitivity file', sensfile)
                if os.path.exists(bstarfile):
                    shutil.copyfile(bstarfile, exdir+bstarbase)
                else:
                    print('WARNING: no bstar file',bstarfile)


    return result


