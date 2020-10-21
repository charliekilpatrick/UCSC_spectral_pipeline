from __future__    import print_function
import os
import util
from astropy.io import fits, ascii

path_to_trunk = os.path.expandvars('$UCSC_SPECPIPE/spectral_reduction/trunk/')
if not os.path.exists(path_to_trunk):
    raise RuntimeError('Error : $UCSC_SPECPIPE variable is undefined or incorrect!')

kast_blue = {'name': 'kast_blue',
             'telescope': 'shane',
             'instrument': 'kast',
             'arm': 'blue',
             'use_ext': 0,
             'read_noise': 3.7,
             'gain': 1.2,
             'grism': 'temp',
             'filter': 'temp',
             'slit': 'temp',
             'dispaxis': 1,
             'biassec': '[1:2048,318:350]',
             # 'trimsec': '[50:2010,7:287]',
             # 'trimsec': '[80:1850,7:287]',
             # 'flatsec': '[80:1850,50:287]',
             'trimsec': '[80:1850,40:330]', #8-31-19
             'flatsec': '[80:1850,80:330]', #8-31-19
             # 'trimsec': '[26:1976,22:302]', # temporary
             'archive_zero_file': path_to_trunk + 'KAST_cals/Zero_blue_20180206.fits',
             'archive_flat_file': path_to_trunk + 'KAST_cals/RESP_blue.fits',
             'archive_sens': path_to_trunk + 'KAST_cals/fluxstarblue.fits',
             'archive_bstar': path_to_trunk + 'KAST_cals/bstarblue.fits',
             # 'archive_arc_extracted': path_to_trunk + 'KAST_cals/Blue_Arc_Ref.ms.fits',
             'archive_arc_extracted_id': path_to_trunk + 'KAST_cals/idARC_blue.ms',
             # 'archive_arc_aperture': path_to_trunk + 'KAST_cals/apBlue_Arc_Ref',
             'line_list': path_to_trunk+'KAST_cals/lines.dat',
             'extinction_file': path_to_trunk + 'KAST_cals/lick_extinction.dat',
             'observatory': 'lick',
             'sky_file': path_to_trunk + 'KAST_cals/licksky.fits',
             'section': 'middle line',
             'pyzap_boxsize': 5, # approximate seeing in pixels
             'pyzap_nsigma': 16,
             'pyzap_subsigma': 2.8,
             'approx_extract_line': 145,
             'approx_extract_column': 1000,
             'pixel_scale': .43, #arcsec/pix
             'spatial_binning': 1.
             }


kast_red = { 'name': 'kast_red',
             'telescope': 'shane',
             'instrument': 'kast',
             'arm': 'red',
             'use_ext': 0,
             'read_noise': 3.8,
             'gain': 1.9,
             'grism': 'temp',
             'filter': 'temp',
             'slit': 'temp',
             'dispaxis': 2,
             'biassec': '[360:405,1:2725]',
             # 'trimsec': '[66:346,80:2296]',
             # 'flatsec': '[110:346,80:2296]',
             'trimsec': '[70:370,80:2296]', #8-31-19
             'flatsec': '[110:370,80:2296]', #8-31-19
             'archive_zero_file': path_to_trunk + 'KAST_cals/Zero_red_20180206.fits',
             'archive_flat_file': path_to_trunk + 'KAST_cals/RESP_red.fits',
             'archive_sens': path_to_trunk + 'KAST_cals/fluxstarred.fits',
             'archive_bstar': path_to_trunk + 'KAST_cals/bstarred.fits',
             # 'archive_arc_extracted': path_to_trunk + 'KAST_cals/Red_Arc_Ref.ms.fits',
             'archive_arc_extracted_id': path_to_trunk + 'KAST_cals/idARC_red.ms',
             # 'archive_arc_aperture': path_to_trunk + 'KAST_cals/apRed_Arc_Ref',
             'line_list': path_to_trunk+'KAST_cals/lines.dat',
             'extinction_file': path_to_trunk + 'KAST_cals/lick_extinction.dat',
             'observatory': 'lick',
             'sky_file': path_to_trunk + 'KAST_cals/licksky.fits',
             'section': 'middle column',
             'pyzap_boxsize': 5,
             'pyzap_nsigma': 16,
             'pyzap_subsigma': 2.8,
             'approx_extract_line': 150,
             'approx_extract_column': 1000,
             'pixel_scale': .43, #arcsec/pix
             'spatial_binning': 1.
             }

lris_blue = {'name': 'lris_blue',
             'telescope': 'keck',
             'instrument': 'lris',
             'arm': 'blue',
             'use_ext': 0,
             'read_noise': 3.7,
             'gain': 1.5,
             'grism': 'temp',
             'filter': 'temp',
             'slit': 'temp',
             'dispaxis': 1,
             'archive_flat_file': path_to_trunk + 'LRIS_cals/RESP_blue.fits',
             'archive_sens': path_to_trunk + 'LRIS_cals/fluxstarblue.fits',
             'archive_bstar': path_to_trunk + 'LRIS_cals/bstarblue.fits',
             'archive_arc_extracted': path_to_trunk + 'LRIS_cals/LRIS_Blue_Arc_Ref.ms.fits',
             'archive_arc_extracted_id': path_to_trunk + 'LRIS_cals/id/idLRIS_Blue_Arc_Ref.ms',
             'archive_arc_aperture': path_to_trunk + 'LRIS_cals/apLRIS_Blue_Arc_Ref',
             'line_list': path_to_trunk+'LRIS_cals/lines.dat', #this nearly hits pathname length limits
             'extinction_file': path_to_trunk + 'LRIS_cals/lick_extinction.dat', # JB: need to fix
             'observatory': 'keck',
             'sky_file': path_to_trunk + 'LRIS_cals/kecksky.fits',
             'section': 'middle line',
             'pyzap_boxsize': 5,
             'pyzap_nsigma': 16,
             'pyzap_subsigma': 8,
             # 'approx_extract_line': 740, jon old
             'approx_extract_line': 200,
             'pixel_scale': .135, #arcsec/pix
             'spatial_binning': 1.
             }

lris_red = { 'name': 'lris_red',
             'telescope': 'keck',
             'instrument': 'lris',
             'arm': 'red',
             'use_ext': 0,
             'read_noise': 4.7,
             'gain': 1.2,
             'grism': 'temp',
             'filter': 'temp',
             'slit': 'temp',
             'dispaxis': 1,
             'archive_flat_file': path_to_trunk + 'LRIS_cals/RESP_red.fits',
             'archive_sens': path_to_trunk + 'LRIS_cals/fluxstarred.fits',
             'archive_bstar': path_to_trunk + 'LRIS_cals/bstarred.fits',
             'archive_arc_extracted': path_to_trunk + 'LRIS_cals/LRIS_Red_Arc_Ref.ms.fits',
             'archive_arc_extracted_id': path_to_trunk + 'LRIS_cals/id/idLRIS_Red_Arc_Ref.ms',
             'archive_arc_aperture': path_to_trunk + 'LRIS_cals/apLRIS_Red_Arc_Ref',
             'line_list': path_to_trunk+'LRIS_cals/lines.dat', #this nearly hits pathname length limits
             'extinction_file': path_to_trunk + 'LRIS_cals/lick_extinction.dat', # JB: need to fix
             'observatory': 'keck',
             'sky_file': path_to_trunk + 'LRIS_cals/kecksky.fits',
             'section': 'middle line',
             'pyzap_boxsize': 5,
             'pyzap_nsigma': 16,
             'pyzap_subsigma': 8,
             # 'approx_extract_line': 370, jon old
             'approx_extract_line': 69,
             'pixel_scale': .135, #arcsec/pix
             'spatial_binning': 2.
             }

goodman_m1={ 'name': 'goodman_blue',
             'telescope': 'soar',
             'instrument': 'goodman',
             'arm': 'blue',
             'use_ext': 0,
             'read_noise': 3.99,
             'gain': 2.06,
             'grism': 'temp',
             'filter': 'temp',
             'slit': 'temp',
             'dispaxis': 1,
             'biassec': '[1:25,100:800]', # RED CAMERA
             'trimsec': '[400:2068,84:900]', # RED CAMERA
             # 'biassec': '[2060:2070,100:800]', # BLUE CAMERA
             # 'trimsec': '[80:2025,65:775]', # BLUE CAMERA
             # 'trimsec': '[400:2025,65:775]', # BLUE CAMERA (cuts out flat field artifact)
             'archive_zero_file': path_to_trunk + 'SOAR_cals/Zero_red_20180206.fits',
             'archive_flat_file': path_to_trunk + 'SOAR_cals/RESP_blue.fits',
             'archive_sens': path_to_trunk + 'SOAR_cals/fluxstarblue.fits',
             'archive_bstar': path_to_trunk + 'SOAR_cals/bstarblue.fits',
             'archive_arc_extracted': path_to_trunk + 'SOAR_cals/m1_Arc_Ref.ms.fits',
             'archive_arc_extracted_id': path_to_trunk + 'SOAR_cals/id/idm1_Arc_Ref.ms',
             'archive_arc_aperture': path_to_trunk + 'SOAR_cals/apm1_Arc_Ref',
             'apline': 'INDEF',
             'line_list': path_to_trunk+'SOAR_cals/lines.dat',
             'extinction_file': path_to_trunk + 'lick_extinction.dat', # JB: need to fix
             'observatory': 'soar',
             'sky_file': path_to_trunk + 'licksky.fits', # JB: need to fix
             'section': 'middle column',
             'pyzap_boxsize': 5,
             'pixel_scale': .15, #arcsec/pix
             'spatial_binning': 1.
             }

goodman_m2={ 'name': 'goodman_red',
             'telescope': 'soar',
             'instrument': 'goodman',
             'arm': 'red',
             'use_ext': 0,
             'read_noise': 3.99,
             'gain': 2.06,
             'grism': 'temp',
             'filter': 'temp',
             'slit': 'temp',
             'dispaxis': 1,
             #'biassec': '[1:25,100:800]', # RED CAMERA
             #'trimsec': '[400:2068,84:900]', # RED CAMERA
             'biassec': '[2060:2070,100:800]', # BLUE CAMERA
             'trimsec': '[80:2025,65:775]', # BLUE CAMERA
             'archive_zero_file': path_to_trunk + 'SOAR_cals/Zero_red_20180206.fits',
             'archive_flat_file': path_to_trunk + 'SOAR_cals/RESP_red.fits',
             'archive_sens': path_to_trunk + 'SOAR_cals/fluxstarred.fits',
             'archive_bstar': path_to_trunk + 'SOAR_cals/bstarred.fits',
             'archive_arc_extracted': path_to_trunk + 'SOAR_cals/m2_Arc_Ref.ms.fits',
             'archive_arc_extracted_id': path_to_trunk + 'SOAR_cals/id/idm2_Arc_Ref.ms',
             'archive_arc_aperture': path_to_trunk + 'SOAR_cals/apm2_Arc_Ref',
             'apline': 500,
             'line_list': path_to_trunk+'SOAR_cals/lines.dat',
             'extinction_file': path_to_trunk + 'SOAR_cals/lick_extinction.dat', # JB: need to fix
             'observatory': 'soar',
             'sky_file': path_to_trunk + 'SOAR_cals/licksky.fits', # JB: need to fix
             'section': 'middle column',
             'pyzap_boxsize': 5,
             'pixel_scale': .15, #arcsec/pix
             'spatial_binning': 1.
             }

binospec={ 'name': 'binospec',
             'telescope': 'mmt',
             'instrument': 'binospec',
             'arm': 'only',
             'use_ext': 2,
             'read_noise': 3.2,
             'gain': 1.0,
             'grism': 'temp',
             'filter': 'temp',
             'slit': 'temp',
             'dispaxis': 1,
             'datasec': '[1:4096,1:4112]',   # IM5
             'biassec': '[1:4096,1:820]',    # IM5
             'trimsec': '[1:4096,1000:2780]', # IM5
             'flatsec': '[1:4096,1000:2780]', # IM5
             'archive_flat_file': path_to_trunk + 'BINO_cals/RESP.fits',
             'archive_sens': path_to_trunk + 'BINO_cals/fluxstarblue.fits',
             'archive_bstar': path_to_trunk + 'BINO_cals/bstarblue.fits',
             'archive_arc_extracted': path_to_trunk + 'BINO_cals/Arc_Ref.ms.fits',
             'archive_arc_extracted_id': path_to_trunk + 'BINO_cals/id/id_Arc_Ref.ms',
             'archive_arc_aperture': path_to_trunk + 'BINO_cals/ap_Arc_Ref',
             'apline': '2083',
             'line_list': path_to_trunk+'BINO_cals/lines.dat',
             'extinction_file': path_to_trunk + 'lick_extinction.dat', # JB: need to fix
             'observatory': 'mmt',
             'sky_file': path_to_trunk + 'licksky.fits', # JB: need to fix
             'section': 'middle column',
             'pyzap_boxsize': 3,
             'pixel_scale': .24, #arcsec/pix
             'spatial_binning': 1.
             }


def blue_or_red(img):
    hdu = fits.open(img)
    hdr = hdu[0].header

    # Add keys as you need them for different instruments
    check_keys = ['VERSION','WAVMODE','INSTRUME']
    check = {}
    for key in check_keys:
        if key in hdr.keys() and key not in check.keys():
            check[key] = hdr[key]

    if len(hdu)>1:
        hdr = hdu[1].header
    for key in check_keys:
        if key in hdr.keys() and key not in check.keys():
            check[key] = hdr[key]


    # kast
    if util.readkey3(check, 'VERSION') == 'kastb':
        return 'blue', kast_blue
    elif util.readkey3(check, 'VERSION') == 'kastr':
        return 'red', kast_red
    # soar
    elif util.readkey3(check, 'WAVMODE') == '400_M1':
        return 'blue', goodman_m1
    elif util.readkey3(check, 'WAVMODE') == '400_M2':
        return 'red', goodman_m2
    # lris
    elif util.readkey3(check, 'INSTRUME') == 'LRISBLUE':
        return 'blue', lris_blue
    elif util.readkey3(check, 'INSTRUME') == 'LRIS':
        return 'red', lris_red
    elif util.readkey3(check, 'INSTRUME') == 'Binospec':
        return 'only', binospec
    else:
        print('Current instrument not in database!')
        return None, None




