def final(objectlist,gratcode,secondord,gratcode2,user):
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import numpy as np
    from datetime import datetime
    import os, pdb
    import inspect
    from astropy.io import fits
    from astropy import units as u
    from astropy.coordinates import SkyCoord
#    from astropy.time import Time
    from tmath.wombat.get_screen_size import get_screen_size
    from tmath.wombat.getmswave import getmswave
    from tmath.wombat.inputter_single import inputter_single
    from tmath.wombat.inputter import inputter
    from tmath.pydux.scipyrebinsky import scipyrebinsky
    from tmath.wombat.yesno import yesno
    from tmath.wombat.womashrebin import womashrebin
    from tmath.pydux.getfitsfile import getfitsfile
    from tmath.pydux.pacheck import pacheck
    from tmath.pydux.jdcnv import jdcnv
    from tmath.pydux.baryvel import baryvel
    from tmath.pydux.finalscaler import finalscaler
    from tmath.pydux.envelope import envelope
    from tmath.pydux.congrid import congrid
    from tmath.pydux.xcor import xcor
    from tmath.pydux.telluric_remove import telluric_remove
    from tmath.pydux.waveparse import waveparse
    from tmath.pydux.parsehourangle import parsehourangle
    from tmath.pydux import RADEG
    plt.ion()

    screen_width, screen_height=get_screen_size()
    screenpos='+{}+{}'.format(int(screen_width*0.2),int(screen_height*0.05))
    fig=plt.figure()

    fig.canvas.manager.window.wm_geometry(screenpos)
    fig.canvas.set_window_title('Final Reductions')
    fig.set_size_inches(18,12)
    # turns off key stroke interaction
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)

    # this should bring the window to the top, but it doesn't
    wm=plt.get_current_fig_manager()
    #fig.canvas.manager.window.attributes('-topmost', 1)
    blah=wm.window.attributes('-topmost', 1)

    #plt.pause(0.2)
    #fig.canvas.manager.window.attributes('-topmost', 0)
    blah=wm.window.attributes('-topmost', 0)


    bstarfits=fits.open('bstar'+gratcode+'.fits')
    bstarstar=bstarfits[0].data
    bstarhead=bstarfits[0].header
    bstarwavezero=float(bstarhead['CRVAL1'])
    bstarwavedelt=float(bstarhead['CDELT1'])
    bstarwave=np.arange(len(bstarstar))*bstarwavedelt+bstarwavezero
    bstarairmass=float(bstarhead['AIRMASS'])
    bstarname=bstarhead['OBJECT']
    try:
        bstarnum=int(bstarhead['OBSNUM'])
    except KeyError:
        bstarnum=0
    if (secondord):
        bstarfits2=fits.open('bstarstar'+gratcode2+'.fits')
        bstarstar2=bstarfits[0].data
        bstarhead2=bstarfits[0].header
        bstarwavezero2=float(bstarhead2['CRVAL1'])
        bstarwavedelt2=float(bstarhead2['CDELT1'])
        bstarwave2=np.arange(len(bstarstar2))*bstarwavedelt2+bstarwavezero2
        bstarairmass2=float(bstarhead2['AIRMASS'])
        bstarname2=bstarhead2['OBJECT']
        try:
            bstarnum2=int(bstarhead2['OBSNUM'])
        except KeyError:
            bstarnum2=0
    observat=bstarhead['OBSERVAT'].strip().lower()
    #reduxdir=os.environ['PY_REDUX']
    reduxdir=inspect.getfile(xcor)
    reduxdir=reduxdir.rsplit('/',1)[0]+'/'
    kecksky=['keck','gemini-north','gemini-n','gemini-south','gemini-s', \
             'soar','ctio','vlt','lco','lco-imacs','lapalma']
    licksky=['lick','kpno','flwo','mmto','sso']
    if (observat in kecksky):
        mskyfile=reduxdir+'kecksky.fits'
    elif (observat in licksky):
        mskyfile=reduxdir+'licksky.fits'
    else:
        print('\nCannot find mastersky file and observatory unknown\n')
        mskyfile=getfitsfile('master sky','.fits')
    print('Using {} as master sky'.format(mskyfile))
    mskyfits=fits.open(mskyfile)
    mskydata=mskyfits[0].data
    mskyhead=mskyfits[0].header
    mskywavezero=float(mskyhead['CRVAL1'])
    mskywavedelt=float(mskyhead['CDELT1'])
    mskywave=np.arange(len(mskydata))*mskywavedelt+mskywavezero
    if (np.abs(np.mean(mskydata)) < 1e-7):
        mskydata=mskydata*1e15
    print('\nHello {}\n'.format(user))
    lat_dict={'keck': 19.8283,
              'gemini-north': 19.8238,
              'gemini-south': -30.228,
              'gemini-n': 19.8238,
              'gemini-s': -30.228,
              'soar': -30.228,
              'kpno': 31.9633,
              'lick': 37.3414,
              'palomar': 33.35611,
              'mcdonald': 30.6717,
              'flwo': 31.681,
              'mmto': 31.688,
              'sso':  -31.2734,
              'vlt':  -24.6254,
              'lco':  -29.01,
              'lco-imacs':  -29.01,
              'lapalma': 28.75833,
              'ctio': -30.16894}
    infile=open(objectlist,'r')
    objname=''
    secondtime = False
    seconddone = False
    for inputfile in infile:
        inputfile=inputfile.strip()
        if ('.fits' not in inputfile):
            inputfile=inputfile+'.fits'
        multifits=fits.open('c'+gratcode+inputfile)
        multispec=multifits[0].data
        mshead=multifits[0].header
        num_apertures=multispec.shape[1]
        num_bands=multispec.shape[0]
        wavearr=np.zeros((multispec.shape[2],multispec.shape[1]))
        if (secondord):
            multifits2=fits.open('c'+gratcode2+inputfile)
            multispec2=multifits2[0].data
            mshead2=multifits2[0].header
        pacheck(mshead)
        objectname=mshead['OBJECT']
        print('The object is: {}'.format(objectname))
        observat=mshead['OBSERVAT'].strip().lower()
        airmass=float(mshead['AIRMASS'])
        if (airmass < 1):
            airmass=1.0
        stringra=mshead['RA']
        stringdec=mshead['DEC']
        geminilist=['gemini-north','gemini-south','gemini-n','gemini-s']
        if (observat) in geminilist:
            ra=float(stringra)
            dec=float(stringdec)
        else:
            coords=SkyCoord(stringra,stringdec,unit=(u.hourangle,u.deg))
            ra=coords.ra.deg
            dec=coords.dec.deg
        ra=ra/RADEG
        dec=dec/RADEG
        haheader=mshead['HA']
        ha=parsehourangle(haheader)
        ha=ha*(360./24.)/RADEG

        if (observat in lat_dict):
            latitude=lat_dict[observat]
        else:
            print('Observatory Unknown!!!!!!!')
            latitude=0.0
        latitude=latitude/RADEG
        # Get Julian date and Earth's velocity
        #epoch=float(mshead['EPOCH'])
        date=mshead['DATE-OBS'].strip()
        # get year,month,day from ISO Y2K format YYYY-MM-DDThh:mm:ss.ssss
        # the time following T is optional, so we get T from UTMIDDLE
        if (date[4] == '-'):
            year=int(date.split('-')[0])
            month=int(date.split('-')[1])
            day=int(date.split('-')[2].split('T')[0])
            printdate='{}{:02d}{:02d}'.format(year,month,day)
        else:
            # Old date format DD/MM/YY
            year=int(date[6:8])
            month=int(date[3:5])
            day=int(date[0:2])
            # try to catch old format written after 2000
            # good until 2050, by which time no one should
            # be doing this
            if (year > 50):
                year=year+1900
            else:
                year=year+2000
            printdate='{}{:02d}{:02d}'.format(year,month,day)
        try:
            ut=mshead['UTMIDDLE'].strip()
        except KeyError:
            ut=mshead['UT'].strip()
        if ('T' in ut):
            ut=ut.split('T')[1]
        # some UT/UTMIDDLE values are ISO format
        hour=int(ut.split(':')[0])
        minute=int(ut.split(':')[1])
        second=float(ut.split(':')[2])
        julian=jdcnv(year,month,day,(hour+(minute/60.)+(second/3600.)))
        vh, vb=baryvel(julian)
        # Earth's velocity toward object
        v=vb[0]*np.cos(dec)*np.cos(ra) + vb[1]*np.cos(dec)*np.sin(ra) \
           + vb[2]*np.sin(dec)
        # Correct for earth's rotation.
        # Note that this isn't strictly correct because ha is set to
        # the beginning of the observation, while it really should be
        # the middle.  But this is a small difference and doesn't affect
        # the results in any way that we would notice...
        v =  v - ( 0.4651 * np.cos(latitude) * np.sin(ha) * np.cos(dec) )
        print('\nThe velocity of the Earth toward the target is {} km/s'.format(v))
        print('\nThere are {} apertures in the spectrum.\n'.format(num_apertures))
        # clean up header
        """        delkeylist=['WAT0_001', 'WAT1_001', 'WAT2_001', 'WAT0_002', \
                'WAT1_002', 'WAT2_002', 'WAT3_001', 'WAT2_003', \
                'CTYPE1', 'CTYPE2', 'CTYPE3', 'CD1_1', 'CD2_2', \
                'CD3_3', 'LTM1_1', 'LTM2_2', 'LTM3_3', 'WCSDIM']
        for k in delkeylist:
            try:
                mshead.remove(k)
            except KeyError:
                pass """
        for i in range(0,num_apertures):
            print('\nAperture {}:'.format(i+1))
            wave=getmswave(mshead,i)
            wdelt=wave[1]-wave[0]
            mean=np.mean(multispec[0,i,:])
            ymin,ymax=finalscaler(multispec[0,i,:])
            plt.clf()
            gs=gridspec.GridSpec(2,1,height_ratios=[4,1])
            ax0=plt.subplot(gs[0])
            ax1=plt.subplot(gs[1])
            fig.subplots_adjust(hspace=0)
            ax0.plot(wave,multispec[1,i,:],drawstyle='steps-mid',color='r')
            ax0.plot(wave,multispec[0,i,:],drawstyle='steps-mid',color='k')
            plt.pause(0.01)
            ax0.set_ylim((ymin,ymax))
            ax1.semilogy(wave,np.abs(multispec[1,i,:]-multispec[0,i,:])/mean, \
                         drawstyle='steps-mid',color='k')
            ax0.set_ylabel('Flux')
            ax1.set_ylabel('Log of fractional residuals')
            ax1.set_xlabel('Wavelength')
            plt.pause(0.01)
            print('\nPlotting optimal as black, normal as red\n')
            extract=inputter_single('Do you want to use the (n)ormal or the (o)ptimal extraction? ','no')
            if (extract == 'o'):
                object=multispec[0,i,:]
                extractcode='optimal'
                if (secondord):
                    object2=multispec2[0,i,:]
            else:
                object=multispec[1,i,:]
                extractcode='normal'
                if (secondord):
                    object2=multispec2[1,i,:]
            skyband=2
            if (num_bands == 2):
                skyband=1
            if (num_bands > 3):
                sigma=multispec[3,i,:]
            else:
                sigma=np.sqrt(multispec[skyband,i,:])
            if (np.abs(np.mean(object)) < 1e-7):
                print('\nScaling data up by 10^15\n')
                object=object*1.e15
                sigma=sigma*1.e15
            if (secondord):
                if (num_bands > 3):
                    sigma2=multispec2[3,i,:]
                else:
                    sigma2=np.sqrt(multispec2[skyband,i,:])
                    if (np.abs(np.mean(object)) < 1e-7):
                        object2=object2*1.e15
                        sigma2=sigma2*1.e15
            nandslist=['gemini-n','gemini-s','lco-imacs']
            if (observat in nandslist):
                skyfile=inputfile.replace('tot','sky')
                skyfits=fits.open(skyfile)
                sky=skyfits[0].data[2,0,:]
            else:
                sky=multispec[skyband,i,:]
            envelope_size=25
            mx,mn=envelope(sky,envelope_size)
            skycont=congrid(mn,(len(sky),),minusone=True)
            sky=sky-skycont
            if (np.abs(np.mean(sky)) < 1.e-7):
                sky=sky*1.e15
            print(len(mskywave),len(mskydata))
            msky=scipyrebinsky(mskywave,mskydata,wave)
            xfactor=10
            maxlag=200
            shift=xcor(msky,sky,xfactor,maxlag)
            angshift=shift*wdelt
            print('wdeltf',wdelt)
            print('The x-cor shift in Angstroms is {}'.format(angshift))
            wave=wave+angshift
            skyshiftdone=False
            npixsky2=len(sky)//2
            while (not skyshiftdone):
                plt.clf()
                axarr=fig.subplots(2)
                fig.subplots_adjust(hspace=0)
                waveplus=wave-angshift
                axarr[0].plot(waveplus[0:npixsky2],msky[0:npixsky2], \
                              drawstyle='steps-mid',color='k')
                axarr[0].plot(wave[0:npixsky2],sky[0:npixsky2], \
                              drawstyle='steps-mid',color='r')
                axarr[1].plot(waveplus[npixsky2-1:],msky[npixsky2-1:], \
                              drawstyle='steps-mid',color='k')
                axarr[1].plot(wave[npixsky2-1:],sky[npixsky2-1:], \
                              drawstyle='steps-mid',color='r')
                plt.pause(0.01)
                print('\nBlack spectrum = master sky')
                print('Red spectrum   = object sky shifted to match master sky')
                print('Is this ok?')
                answer=yesno('y')
                if (answer == 'n'):
                    wave=wave- angshift
                    angshift=inputter('Enter desired shift in Angstroms: ','float',False)
                    wave=wave+angshift
                else:
                    skyshiftdone = True
            # B star removal
            bstarpass=bstarstar
            bobj, bsig, bangshift=telluric_remove(bstarwave,bstarpass, bstarairmass, wave, \
                                       object, airmass, sigma)
            if (secondord):
                bobj2, bsig2, bangshift2 =telluric_remove(bstarwave2,bstarpass2, bstarairmass2, \
                                             wave, object2, airmass, sigma2)

            print('\nRemoving redshift due to motion of Earth...')
            z=-1*v/2.99792458e5
            wave=wave/(1+z)
            # rebin
            plt.clf()
            axarr=fig.subplots(1,2)
            ymin,ymax=finalscaler(bobj[0:100])
            axarr[0].plot(wave[0:100],bobj[0:100],drawstyle='steps-mid')
            axarr[0].set_xlabel('Wavelength')
            axarr[0].set_ylabel('Flux')
            #axarr[0].set_ylim((ymin,ymax))
            if (secondtime):
                axarr[0].plot([wavesave0,wavesave0],[ymin,ymax],color='r')
            plt.pause(0.01)
            ymin,ymax=finalscaler(bobj[-100:])
            axarr[1].plot(wave[-100:],bobj[-100:],drawstyle='steps-mid')
            axarr[1].set_xlabel('Wavelength')
            axarr[1].set_ylabel('Flux')
            try:
                axarr[1].set_ylim((ymin,ymax))
            except ValueError:
                bobj_median = np.nanmedian(bobj)
                axarr[1].set_ylim((bobj_median/10.,bobj_median*10.))
            if (secondtime):
                axarr[1].plot([wavesaven,wavesaven],[ymin,ymax],color='r')

                newdelt=0.0
            plt.pause(0.01)
            print('\nCurrent A/pix is {}'.format(wave[1]-wave[0]))
            if (secondtime):
                print('\nPrevious resolution choice: {}'.format(deltsave))
            else:
                deltsave = 0.0
            newdelt=0.0
            while (newdelt <= 0) or (newdelt > wave[-1]):
                print('Rebin to how many Angstroms per pixel? ')
                newdelt=inputter('         <CR> selects previous choice: ','float',True,deltsave)
                if (newdelt <= 0) or (newdelt > wave[-1]):
                    print('Need positive resoution and smaller than the')
                    print('entire spectrum.  Try again')
            print('\nCurrent range: {} {}'.format(wave[0],wave[-1]))
            if (secondtime):
                print('\nPrevious selection was {} {} (marked in red on plot)'.format(wavesave0,wavesaven))
                print('\n<CR> selects previous choice\n')
            else:
                wavesave0 = 0.0
                wavesaven = 0.0
                print('Enter the new wavelength range desired: ')
            
            waveb,waver=waveparse(wave,wavesave0,wavesaven)
            newbin=(waver-waveb)/newdelt +1.0
            frac,whole=np.modf(newbin)
            if (frac > 0.000001):
                print('NON-INTEGER number of bins')
                testwave=newdelt*whole+waveb
                print('Closest match is: {} to {}'.format(waveb,testwave))
                print('Would you like this wavelength range?')
                answer=yesno('y')
                if (answer == 'y'):
                    waver=testwave
            deltsave=newdelt
            wavesave0=waveb
            wavesaven=waver
            waverange=str(waveb) + ' ' + str(waver)
            secondtime = True
            nwave=np.arange(0,whole)*newdelt + waveb
            vartmp=bsig*bsig
            olddeltvec=wave-np.roll(wave,1)
            olddelt=np.mean(olddeltvec)
            finalobj=womashrebin(wave,bobj,nwave)
            finalvar=womashrebin(wave,vartmp,nwave)
            finalsig=np.sqrt(finalvar)*np.sqrt(olddelt/newdelt)
            if (secondord):
                vartmp2=bsig2*bsig2
                finalobj2=womashrebin(wave,bobj2,nwave)
                finalvar2=womashrebin(wave,vartmp2,nwave)
                finalsig2=np.sqrt(finalvar2)*np.sqrt(olddelt/newdelt)
            if (secondord and not seconddone):
                finalobj, finalsig, wavebeg, waveend, brscale=secondcat(nwave,finalobj,finalobj2, finalsig,finalsig2,secondtime, wavebeg, waveend, brscale)
            seconddone = True
            ymin,ymax=finalscaler(finalobj)
            plt.clf()
            plt.plot(nwave,finalobj,drawstyle='steps-mid')
            plt.xlabel('Wavelength')
            plt.ylabel('Flux')
            plt.title(objectname)
            outputdone = False
            while (not outputdone):
                print('\nThe file is: {}'.format(inputfile))
                print('The object is: {}'.format(objectname))
                print('The DATE-OBS is: {}'.format(date))
                print('The aperture is: {}'.format(i+1))
                print('The previous name was: {}'.format(objname))
                print('\nEnter the object name for the final fits file: ')
                objname=inputter('(UT date and .fits will be added): ','string',False)
                fname=objname+'-'+printdate+'.fits'
                sname=objname+'-'+printdate+'-sigma.fits'
                if (os.path.isfile(fname)):
                    print('{} already exists!!!!'.format(fname))
                    print('Do you wish to overwrite it? ')
                    answer=yesno('y')
                    if (answer == 'y'):
                        outputdone = True
                else:
                    outputdone = True
            # add to header
            mshead.set('CRPIX1', 1)
            mshead.set('CRVAL1',  nwave[0])
            mshead.set('CDELT1', nwave[1] - nwave[0])
            mshead.set('CTYPE1', 'LINEAR')
            mshead.set('W_RANGE', waverange)
            mshead.set('BSTAR_Z', bstarairmass)
            mshead.set('BSTARNUM', bstarnum)
            mshead.set('BSTAROBJ', bstarname)
            mshead.set('BARYVEL', v)
            mshead.set('SKYSHIFT', angshift)
            mshead.set('ATMSHIFT', bangshift)
            mshead.set('EXTRACT', extractcode)
            mshead.set('REDUCER', user)
            mshead.set('RED_DATE', datetime.now().strftime("%Y-%M-%d %I:%M%p"), 'EPOCH OF REDUCTION')
            mshead.set('OBJECT', objectname)
            if (secondord):
                mshead.set('SECOND', 'yes',  'Second order correction attempted')
                mshead.set('COMBRANGE',combrange)
                mshead.set('BSTAR_Z2',bstarairmass2)
                mshead.set('BSTARNU2', bstarnum2)
                mshead.set('BSTAROB2', bstarname2)
                fluxairmass2=mshead2['FLX2_Z']
                fluxnum2=mshead2['FLX2_NUM']
                fluxname2=mshead2['FLX2_OBJ']
                mshead.set('FLX2_Z',fluxairmass2)
                mshead.set('FLX2_NUM', fluxnum2)
                mshead.set('FLX2_OBJ', fluxname2)
            outdata=np.zeros((len(finalobj),2))
            outdata[:,0]=finalobj.copy()
            outdata[:,1]=finalsig.copy()
            outhdu=fits.PrimaryHDU(outdata)
            hdul=fits.HDUList([outhdu])
            mshead.set('NAXIS2',2)
            hdul[0].header=mshead.copy()
            hdul.writeto(fname,overwrite=True)
            hdul.close()
                
    print('final')        
    print(objectlist,gratcode,secondord,gratcode2)
    return
