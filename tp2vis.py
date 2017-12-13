#
# A collection of TP2VIS functions to aid in combining ALMA
# Total Power and Visibilities in a Joint Deconvolution
#
# Public functions:
#    tp2vis_version()
#    tp2vis(imagename, msname, ptg, maxuv=10.0, rms=None, nvgrp=4, deconv=True)
#    tp2viswt(msTP, int_ms, tp_beam, int_beam, rms, factor)
#
# Helper functions:
#    getptg()
#    axinorder()
#    arangeax()
#
# PJT functions:    see pjt.py
# PLOT functions:   see tp2visplot.py
# These and regressions tests are available from the full distribution in
# https://github.com/kodajn/tp2vis
#

import os, sys, shutil, time
import numpy as np
import matplotlib.pyplot as plt

## ========================
## Define support functions
## ========================

    
def tp2vis_version():
    print "version 0.6 (12-dec-2017)"

   
def axinorder(image):
    """
        Ensure we have the image in RA-DEC-POL-FREQ axes order.
        Helper function for tp2vis()
    """
    ia.open(image)
    h0 = ia.summary()
    ia.done()
    aname = h0['axisnames']
    print "AXIS NAMES:",aname
    print "AXIS SHAPES:",list(h0['shape'])


    order = ''
    for ax in ['Right Ascension','Declination','Stokes','Frequency']:
        if not ax in aname:
            print "ERROR: No %s axis in %s" % (ax,image)
            raise
        else:
            iax   = list(aname).index(ax)
            order = order + '%1d' % (iax)

    if order == '0123':
        return True
    else:
        return False

def arangeax(image):
    """
        Re-arrange axes and make a RA-DEC-POL-FREQ cube. assume
        axinorder() is already run and 4 axes exist in order.
        Helper function for tp2vis()    
    """
    imageout = 'tmp_trans_' + image        # danger: using fixed name
    os.system('rm -rf %s' % imageout)
    ia.open(image)
    h0 = ia.summary()
    aname = h0['axisnames']

    order = ''
    for ax in ['Right Ascension','Declination','Stokes','Frequency']:
        if ax in aname:
            iax   = list(aname).index(ax)
            order = order + '%1d' % (iax)
    if len(order) ==4:                         # all axes exist
        # on older CASA before 5.0 you will loose beam and
        # object name (bugs.txt #017)
        ia2 = ia.transpose(outfile=imageout,order=order)
        ia2.done()
        ia.done()
        print "Written transposed ",imageout
    else:
        return

    return imageout


def getptg(pfile):
    """ get the ptg's (CASA pointings) from a ptg file into a list                                                            
    'J2000 19h00m00.00000 -030d00m00.000000',...
    required for the simulation tools
            Helper function for tp2vis()
    """
    fp = open(pfile)
    lines = fp.readlines()
    fp.close()
    ptg = []
    for line in lines:
        if line[0] == '#': continue
        w = line.strip().split()
        ptg.append("%s %s %s" % (w[0],w[1],w[2]))
    return ptg

## ==========================================================
## TP2VIS: main function to convert TP cube into visibilities
## ==========================================================

def tp2vis(infile, outfile, ptg, maxuv=10.0, rms=None, nvgrp=4, deconv=True):
    """
    Required:
    ---------
    infile    Input IM filename,  in Jy/beam.  Must exist
    outfile   Output MS filename. Must not exist
    ptg       this can be one of:
              None:     NOT ALLOWED ANYMORE UNTIL RE-IMPLEMENTED TO AUTO_FILL
              string:   ptg file  (or list of strings?) - see also pjt_ms_ptg()
              list:     list of strings (from e.g. pjt_ms_ptg())
              ms:       MS containing the pointings
    Optional:
    ---------
    maxuv     maximum uv distance of TP vis distribution (in m)
              default=10m for 12m ALMA dish

    rms       set the RMS (by default this will be ignored) in the infile
              in order to compute an initial guess for the weights
              Should be Jy/beam
    nvgrp     Number of visibility group (nvis = 1035*nvgrp)
              The number of antenna is hardcoded as 46
    deconv    Use deconvolution as input model (True) -
              you almost never want to change this

    """

    # Bug fixes
    # =========

    use_vp       = False
    bug001_Fixed = False

    # Useful constants
    # ================
    cms = qa.constants('c')['value']            # Speed of light in m/s
    apr = 180.0 * 3600.0 / np.pi                # arcsec per radian

    # Parameters
    # ==========

    seed            = 123                       # for random number

    # Total Power (TP) dish parameters
    tp_beamFWHM0   = 56.6                       # arcsec, TP beam at ref freq.
    tp_beamFreq0   = 115.2e9                    # Hz, ref freq of TP beam

    # Virtual Interferometer (VI) parameters
    if use_vp:              # use vp.setpbgauss(telescope='VIRTUAL',....)
        vi_antname      = 'VIRTUAL'             # should adjust these params.
        vi_dish         = 10.0
        vi_beamFWHM0    = 56.6
        vi_beamFreq0    = 115.2e9
    else: 
        vi_antname      = 'ALMA'                # use ALMA
        vi_dish         = 12.0                  # meters, VI dish size
        vi_beamFWHM0    = 56.6                  # arcsec, VI beam at ref freq.
        vi_beamFreq0    = 115.2e9               # Hz, ref freq of VI beam

    # Query the input image
    # =====================

    # Check if exists
    if os.path.isfile(outfile):
        print "Cannot overwrite",outfile
        return

    # Ensure RA-DEC-POL-FREQ axis order (CASA simulator needs it)
    if axinorder(infile):                       # if 4 axes in order
        imagename = infile                      # use original file
    else:                                       # if not, rearrange
        imagename = arangeax(infile)            # and use re-aranged data

    fax   = 3                                   # freq axis=3

    # Parameters from header
    h0        = imhead(imagename,mode='list')
    shape     = h0['shape']
    nx        = h0['shape'][0]                  # RA
    ny        = h0['shape'][1]                  # DEC
    dx        = np.abs(h0['cdelt1'])            # in radians!
    dy        = np.abs(h0['cdelt2'])
    objname   = h0['object']                    # object name
    nchan     = h0['shape'][3]                  # # of channels
    fstart    = h0['crval4']-h0['crpix4']*h0['cdelt4']  # Hz, channel edge
    fwidth    = h0['cdelt4']                    # Hz, freq width
    reffreq   = fstart + 0.5*fwidth             # Hz, central freq
    reflambda = cms / reffreq                   # wavelength [m]
    refcode   = h0['reffreqtype']               # e.g. 'LSRK'
    bunit     = h0['bunit'].upper()             # JY/BEAM or JY/PIXEL

    # Parameters for TP beam
    # ======================

    # Calculate deconvolution beam (TP beam)
    tp_beamFWHM    = tp_beamFWHM0*(tp_beamFreq0/reffreq) # arcsec at obs freq
    tp_beamSigma   = tp_beamFWHM/apr /(2*np.sqrt(2*np.log(2.))) # radian 
    tp_beamSigmaFT = 1.0/(2.0*np.pi*tp_beamSigma)  # sigma in fourier land
    print "tp_sigma, tp_sigmaFT: ",tp_beamSigma,tp_beamSigmaFT

    # Number of pixels per TP beam
    apixel   = np.abs((dx*apr)*(dy*apr))        # asec2, area in pixel
    abeam    = (np.pi/(4*np.log(2.)))*tp_beamFWHM**2 # asec2, TP bm area
    nppb     = abeam/apixel                     # To convert Jy/bm to Jy/pix
    print "Number of pixels per beam:",nppb

    # Cutoff of deconvolution beam
    eps      = 0.001                            # cutoff amp of gauss tail
    uvMax    = np.sqrt(-2.0*tp_beamSigmaFT**2*np.log(eps)) # & uvmax there
    uvMax    = np.minimum(maxuv/reflambda,uvMax) # compare with maxuv
    print "UVMAX:", uvMax/1000.0,"kLambda"

    # Parameters for VI beam
    # ======================

    vi_beamFWHM    = vi_beamFWHM0*(tp_beamFreq0/reffreq)        # arcsec
    vi_beamSigma   = vi_beamFWHM/apr /(2*np.sqrt(2*np.log(2.))) # radian
    vi_beamSigmaFT = 1.0/(2.0*np.pi*vi_beamSigma)  # sigma in fourier land
    print "vi_sigma, vi_sigmaFT: ",vi_beamSigma,vi_beamSigmaFT

    # Generate fake primary beam for virtual interferometric observations
    #    this beam does not need to be the beam of TP antennas
    # ===================================================================
    
    # make a vpmanager beam for us
    # BUG?  Despite that we use vp.save and vp.load in pjt_tp() before tclean()
    #       and even vp.summary says the 12M is present, tclean claims
    #       VPskyjones says "PB used ALMA"
    #       Kumar says to use tclean's vptable= argument
    #       Remy says to use VIRTUAL for 
    #
    if use_vp:
        print "Using vpmanager VIRTUAL beam"
        vp.reset()
        vp.setpbgauss(telescope='VIRTUAL',
                      halfwidth='1.5arcsec',     # PJT test 
                      maxrad='2.5arcmin',
                      reffreq='115.2GHz',
                      dopb=True)
        vp.summarizevps()

    # Obtain pointing coordinates
    # ===========================

    print "Using ptg = ",ptg
    pointings = getptg(ptg)


    # Deconvolution of TP cube (images) by TP beam
    #   (if deconv=False the input image will be used instead)
    # ========================================================

    # Open TP cube
    ia.open(imagename)                         

    # Generate uvdist^2 image
    du        = 1.0/(nx*dx)
    dv        = 1.0/(ny*dy)
    ugrd,vgrd = np.meshgrid(np.arange(ny),np.arange(nx))
    ugrd      = (ugrd-0.5*(nx-1.0))*du
    vgrd      = (vgrd-0.5*(ny-1.0))*dv
    uvgrd2    = ugrd**2+vgrd**2
    uvgrd2    = np.fft.fftshift(uvgrd2)                   # peak at corner

    del ugrd,vgrd

    # Output deconvolved cube
    if not deconv: 
        print "No deconvolution done"
        imagedecname = ''
    else:
        imagedecname = 'tmp_imagedec.im'
        ia2 = ia.newimagefromimage(imagename,imagedecname,overwrite=True)

        # Loop over channels
        print "Deconvolution loop starts"
        for iz in range(nchan):

            # Beam in Fourier domain
            freq        = fstart+fwidth*(0.5+iz)          # Hz, centfreq of chan
            beamSigmaFT = tp_beamSigmaFT * freq/reffreq
            beamFT      = np.exp(-uvgrd2/(2.0*beamSigmaFT**2))

            # Channel image to be deconvolved
            image       = ia.getchunk([-1,-1,-1,iz],[-1,-1,-1,iz])
                                                          # image[ix][iy][0][0]
            image       = image / nppb                    # scale to Jy/pixel
            imageFT     = np.fft.fft2(image,axes=(0,1))
            nnx         = imageFT.shape[0]
            nny         = imageFT.shape[1]
            imageFT     = imageFT.reshape((nnx,nny))      # remove 3rd & 4th ax

            # Deconvolution 
            imageFTdec       = imageFT.copy()
            idx0             = (uvgrd2   > (uvMax**2))    # idx of outer uv
            idx1             = np.logical_not(idx0)       # idx of inner uv
            imageFTdec[idx1] = imageFT[idx1]/beamFT[idx1] # just for inner uv
            imageFTdec[idx0] = 0.0                        # set outer uv zero
            imagedec         = np.fft.ifft2(imageFTdec)
            ia2.putchunk(abs(imagedec), blc=[0,0,0,iz])

        print "Done. nnx,nny=",nnx,nny
        ia2.close()            

    ia.close()

    # List parameters for virtual interferometric obs
    # ===============================================

    # Due to CASA construction, we cannot set some params directly
    # and have to define many indirect params. E.g., Nvis cannot be set,
    # but is calculated as Nvis=npair*(ttot/tint), where npair=num of ant
    # pairs, ttot=total integ time, and tint=integ time per vis.

    nant            = 46                        # # of fake antennas
    npair           = nant*(nant-1)/2           # # of baselines
    nvis            = npair * nvgrp             # # of vis per point
    source          = objname                   # object name
    npnt            = len(pointings)

    # Spectral windows
    spw_nchan       = nchan                     # # of channels
    spw_fstart      = fstart/1e9                # GHz, start freq
    spw_fwidth      = fwidth/1e9                # GHz, freq width
    spw_fresolution = spw_fwidth
   
    spw_fband       = 'bandtp'                  # just a name
    spw_stokes      = 'I'                       # 1 pol axis (or e.g. 'XX YY')
    spw_refcode     = refcode                   # e.g. 'LSRK'

    # Feed
    fed_mode        = 'perfect X Y'    
    fed_pol         = ['']
  
    # Fields
    fld_calcode     = 'OBJ'
    fld_distance    = '0m'                      # infinite distance

    # Observatory (note the peculiar order where 'ALMA' is needed to get
    # obs_obspos, but then we may to need to reset it to 'VIRTUAL'
    # Is that abuse/overload of the ALMA name?
    # MS->OBSERVATION->TELESCOPE_NAME
    obs_obsname     = 'ALMA'
    obs_obspos      = me.observatory(obs_obsname) # observatory coordinate
    #if use_vp:
    #    obs_obsname = 'VIRTUAL'                   # vp needs this, otherwise ALMA is the pb
        
    # Telescopes
    tel_pbFWHM      = vi_beamFWHM0*(vi_beamFreq0/(spw_fstart*1e9)) # arcsec
    tel_mounttype   = 'alt-az'
    tel_coordsystem = 'local'                   # coordinate of antpos
    tel_antname     = vi_antname
    tel_dish        = vi_dish

    # Fake antenna parms
    tel_antposx=tel_antposy=tel_antposz=np.arange(nant)*1000.0
                                                # fake ant positions
    tel_antdiam     = [tel_dish] * nant         # all dish sizes the same

    # Numbers of vis per pointing and for all pointings
    tvis            = 1.0                       # integ time per vis
    tpnt            = tvis * nvgrp              # integ time per ptgs
    ttot            = tpnt * npnt               # tot int time for all pnt
    tstart          = -ttot/2                   # set start/end time of obs,
    tend            = +ttot/2                   # so that (tend-tstart)/tvis
                                                # = num of vis per point

    # Print parameters & set them in CASA
    # ===================================

    print "TP2VIS Parameters"
    print "   input image name:                    %s" % (imagename)
    print "   image shape:                         %s" % (repr(shape))
    print "   output measurement set name:         %s" % (outfile)
    print "   number of pointings:                 %d" % (npnt)
    print "   number of visibilities per pointing: %d" % (tpnt*npair)
    print "   start frequency [GHz]:               %f" % (spw_fstart)
    print "   frequency width [GHz]:               %f" % (spw_fwidth)
    print "   frequency resolution [GHz]:          %f" % (spw_fresolution)
    print "   freq axis [0=first]:                 %d" % (fax)
    print "   freq channels:                       %d" % (spw_nchan)
    print "   polarizations:                       %s" % (spw_stokes)
    print "   antenna name:                        %s" % (tel_antname)
    print "   VI primary beam fwhm [arcsec]:       %f" % (tel_pbFWHM)
    print "   VI primary beam sigmaFT [m]:         %f" % (vi_beamSigmaFT*reflambda)
    print "   ttot                                 %f" % (ttot)
    print "   frame:                               %s" % (spw_refcode)
    print "   seed:                                %d" % (seed)

    spw_fstart        = str(spw_fstart)      + 'GHz'
    spw_fwidth        = str(spw_fwidth)      + 'GHz'
    spw_fresolution   = str(spw_fresolution) + 'GHz'

    if seed >= 0:
        np.random.seed(seed)
  
    sm.open(outfile)

    if use_vp:                                  # turn on once CASA settled
        vptable = outfile + '/TP2VISVP'
        vp.saveastable(vptable)
    else:
        vptable = None

    sm.setconfig(telescopename=obs_obsname,      # VIRTUAL (ALMA would cause it to use ALMA pb)
                referencelocation=obs_obspos,    # ALMA
                antname=tel_antname,             # VIRTUAL
                mount=tel_mounttype,
                coordsystem=tel_coordsystem,
                x=tel_antposx,y=tel_antposy,z=tel_antposz,
                dishdiameter=tel_antdiam)

    sm.setspwindow(spwname=spw_fband,
                freq=spw_fstart,
                deltafreq=spw_fwidth,
                freqresolution=spw_fresolution,
                nchannels=spw_nchan,
                refcode=spw_refcode,
                stokes=spw_stokes)

    sm.setfeed(mode=fed_mode,
                pol=fed_pol)

    for k in xrange(0,npnt):
        this_pointing = pointings[k]
        src = source + '_%d' % (k)
        sm.setfield(sourcename=src,
                sourcedirection=this_pointing,
                calcode=fld_calcode,
                distance=fld_distance)

    sm.setlimits(shadowlimit=0.001,
                elevationlimit='10deg')

    sm.setauto(autocorrwt=0.0)

    sm.settimes(integrationtime=str(tvis)+'s',
                usehourangle=True,
                referencetime=me.epoch('utc', 'today'))


    # Generate (empty) visibilities
    # =============================

    # This step generates (u,v,w), based on target coord and ant pos
    # following current CASA implementation, but (u,v,w) will be
    # replaced in the next step.
    
    print "Running sm.observemany"
    sources    = []
    starttimes = []
    stoptimes  = []

    tstart_src = tstart
    tend_src   = tstart_src + tpnt
    
    for k in xrange(npnt):
        src = source + '_%d' % (k)
        sources.append(src)
        starttimes.append(str(tstart_src)+'s')
        stoptimes.append(str(tend_src)+'s')
        tstart_src = tstart_src + tpnt
        tend_src   = tstart_src + tpnt

    sm.observemany(sourcenames=sources,
            spwname=spw_fband,
            starttimes=starttimes,
            stoptimes=stoptimes)

    # Genarate (replace) (u,v,w) to follow Gaussian
    # =============================================

    pbsigmaFT       = vi_beamSigmaFT*reflambda  # sigmaF=D/lambda -> D [m]

    # Include (u,v) = (0,0)
    uu = np.array([0.0])
    vv = np.array([0.0])

    # Rest follows Gaussian distribution with < uvMax^2
    nuv = 1                                     # (0,0) exists already
    uvMax2 = (uvMax*reflambda)**2               # 1/lambda -> meter
    while (nuv<nvis):                           # loop until enough
        nrest = nvis-nuv
        utmp,vtmp = np.random.normal(scale=pbsigmaFT,size=(2,nrest))
        ok    = utmp**2+vtmp**2 < uvMax2        # generate gauss and
        uu    = np.append(uu,utmp[ok])          # ok for uvdist<uvMax
        vv    = np.append(vv,vtmp[ok])
        nuv   = uu.size

    uu = uu[:nvis]
    vv = vv[:nvis]
    ww = np.zeros(nvis)

    del utmp,vtmp, ok

    # Replicate the same uv set for all pointings
    if (npnt > 1):
        uu = np.ravel([uu,]*npnt)
        vv = np.ravel([vv,]*npnt)
        ww = np.ravel([ww,]*npnt)

    nuvw = uu.shape[0]
    tb.open(outfile,nomodify=False)
    uvw  = tb.getcol('UVW')
    print "UVW shape",uvw.shape,nuvw,uvw[:,1]
    if len(np.ravel(uvw)) > 0:
        nrow = uvw.shape[1]
        if nrow == nuvw:
            uvw = np.array([uu,vv,ww])
            tb.putcol('UVW',uvw)
            print "UVW0",uu[1],vv[1],ww[1]
        else:
            print "Bad UVW",nrow,nuvw
    else:
        print "WARNING: no uvw?"


    # Set WEIGHT column temporarily
    # =============================

    # Deal with the WEIGHT column in the MS table
    # Adjust weights
    #   weight of individual vis = sqrt(nvis)*rmsJy
    #   after natural wtg --> sqrt(nvis)*rmsJy / sqrt(nvis) = rmsJy

    if rms != None:
        print "Adjusting the weights using rms = %g  nvis=%d" % (rms,nvis)
        weight = tb.getcol('WEIGHT')
        w = rms * np.sqrt(nvis)
        w = 1.0/(w*w)
        print "WEIGHT: Old=%s New=%g Nvis=%d" % (weight[0,0],w,nvis)
        weight[:,:] = w                         # weight[npol,nvis*npnt]
        tb.putcol('WEIGHT',weight)
    tb.close()

    del uvw

    # Fill vis amp/phase based on deconvolved TP image
    # ================================================

    print "setdata"                             # Get all fields
    sm.setdata(fieldid=range(0,npnt))

    print "setvp",vptable                       # Set primary beam
    if use_vp:                                  # turn on once CASA settled
        # according to Kumar
        sm.setvp(dovp=True,usedefaultvp=False, vptable=vptable)
    else:
        sm.setvp(dovp=True,usedefaultvp=False)
    vp.summarizevps()

    print "predict"                             # Replace amp/pha - key task
    if deconv:
        sm.predict(imagename=imagedecname)      # deconvolved cube
    else:
        sm.predict(imagename=imagename)         # input TP cube

    # Print Summary
    sm.summary()

    # Close measurement set
    # =====================
    sm.done()


    # Corrections on CASA header
    #   Most of following should not be necessary if CASA has no bug
    # ==============================================================


    if not bug001_Fixed:
        print "Correcting CASA header inconsistencies [due to CASA bugs]"
        # REST_FREQUENCY in /SOURCE
        #    Having REST_FREQUENCY in header does not make sense (since
        #    multiple lines, but CASA MS does have it. So, put it in.
        h0       = imhead(imagename,mode='list')

        if 'restfreq' in h0.keys():
            restfreq = h0['restfreq'][0]        # restfreq from image header
        else:
            if h0['cunit4'] == 'Hz':            # set it to ref freq [Hz]
                restfreq = h0['crval4']
            elif h0['cunit4'] == 'MHz':
                restfreq = h0['crval4'] * 1.0e6
            elif h0['cunit4'] == 'GHz':
                restfreq = h0['crval4'] * 1.0e9

        print "SET RESTFREQ:::",restfreq/1e9," GHz"
        print "   Set restfreq option in (t)clean if this restfreq is incorrect"

        tb.open(outfile + '/SOURCE',nomodify=False)
        rf = tb.getcol('REST_FREQUENCY')
        rf = rf * 0 + restfreq
        tb.putcol('REST_FREQUENCY',rf)
        tb.close()

        # REF_FREQUENCY in /SPECTRAL_WINDOW
        #    Not clear what should be in this data column, but all ALMA data
        #    seem to have REF_FREQUENCY = REST_FREQUENCY, so we follow.
        tb.open(outfile + '/SPECTRAL_WINDOW',nomodify=False)
        rf = tb.getcol('REF_FREQUENCY')
        rf = rf * 0 + restfreq
        tb.putcol('REF_FREQUENCY',rf)
        tb.close()

    #-end of tp2vis()


## =======================================================
## TP2VISWT: Explore different weights for TP visibilities
## =======================================================
           
def tp2viswt(msTP=None, ms07=None, ms12=None, tp_beam=None, int_beam=None, rms=None, value=1.0, mode=0):
    """

    Parameters
    -----------
    Need one of msTP, ms07, and ms12
    msTP       MS to report weights, or compute new weights for
    
    ms07       MS for  7m, needed to accumulate weights for msTP weights
    ms12       MS for 12m, needed to accumulate weights for msTP weights

    mode       0 or 'statistics'   (default)
               1 or 'constant'
               2 or 'multiply'
               3 or 'rms'
               4 or 'beammatching'

    value      (mode=1,2,3) value to be applied applied to weight

    tp_beam    either a beamsize in arcsec, or a CASA image
    int_beam   either a beamsize in arcsec, or a CASA image (image option
               not supported yet)

    Example usage:
    --------------

    Report weights for "v1.ms"  (TP or INT allowed)
    > tp2viswt("v1.ms",mode=0)              
                       
    Set weights to 0.01
    > tp2viswt("v1.ms",mode=1,value=0.01)               
  
    Multiply weights by 3
    > tp2viswt("v1.ms",mode=2,value=3.0)

    Report beta for beam size weights
    > tp2viswt(tp_beam='dirtymap.TP.psf',int_beam=2.0)

    Apply beam size weights to "tp.ms"
    > tp2viswt("tp.ms","v1.ms","v2.ms",'dirtymap.TP.psf',2.0,mode=4)
      
    """

    # Parameters
    # ----------

    # Lists of measurement sets
    msINT = []
    if ms12 != None:
        if type(ms12) != type([]): ms12 = [ms12]
        for ims in ms12:
            msINT.append(ims)
    if ms07 != None:
        if type(ms07) != type([]): ms07 = [ms07]
        for ims in ms07:
            msINT.append(ims)

    mslist = list(msINT)
    if msTP != None:
        if type(msTP) != type([]): msTP = [msTP]
        for ims in msTP:
            mslist.append(ims)

    if len(msINT) == 0:
        msINT  = None

    if len(mslist) == 0:
        mslist = None

    # Translate mode if set with chars
    if type(mode) is str:
        if   'sta' in mode.lower(): mode = 0 # statistics
        elif 'con' in mode.lower(): mode = 1 # constant
        elif 'mul' in mode.lower(): mode = 2 # multiply
        elif 'rms' in mode.lower(): mode = 3 # rms-based
        elif 'bea' in mode.lower(): mode = 4 # beam matching
        elif 'plo' in mode.lower(): mode = 6 # plot

    # Calculate WEIGHT max & min (default)
    # --------------------------
    if (mode == 0):      # "statistics"
        print "TP2VISWT: statistics of weights"

        print "%20s %12s %13s %13s %13s %13s" \
            % ('MS','nvis','min','max','mean','std')
        for ims in mslist:
            tb.open(ims)
            weight = tb.getcol('WEIGHT')
            nvis   = weight.size
            imin   = weight.min()
            imax   = weight.max()
            imean  = weight.mean()
            istd   = weight.std()
            print "%20s %12d %13.6f %13.6f %13.6f %13.6f" \
                % (ims,nvis,imin,imax,imean,istd)
            tb.close()

        return

    # Set WEIGHT to const.
    # --------------------
    elif (mode == 1):      # "constant"
        print "TP2VISWT: set the weights = %g" % (value)

        for ims in mslist:
            tb.open(ims,nomodify=False)
            weight = tb.getcol('WEIGHT')
            nvis=weight.size
            w = value
            print "%20s     Old=%s New=%g Nvis=%d" % (ims,weight[0,0],w,nvis)
            weight[:,:] = w           # weight[npol,nvis]
            tb.putcol('WEIGHT',weight)
            tb.close()

        return

    # Multipy constant to current WEIGHT
    # ----------------------------------
    elif (mode == 2):      # "multiply"
        print "TP2VISWT: multiply the weights by %g" % (value)

        for ims in mslist:
            tb.open(ims,nomodify=False)
            weight = tb.getcol('WEIGHT')
            nvis=weight.size
            print "%20s      Old=%s New=%g Nvis=%d" \
                % (ims,weight[0,0],weight[0,0]*value,nvis)
            weight = value * weight
            tb.putcol('WEIGHT',weight)
            tb.close()

        return

    # From RMS in TP cube and Nvis
    # ----------------------------
    elif (mode == 3):      # "rms-based"
        print "TP2VISWT: adjust the weights using rms = %g" % (value)

        if value == None:
            print "WEIGHT mode=3  rms not given. set value=rms"
            return

        # we need a proper nvis per field:
        # in tp2vis we observe every field, nvis from the MS needs to be
        # reduced by # fields (not the POINTINGS table, the SOURCE or FIELD)

        for ims in mslist:
            tb.open(ims + '/FIELD')
            src_id = tb.getcol('SOURCE_ID')
            npnt = len(src_id)
            tb.close()
        
            tb.open(ims,nomodify=False)
            weight = tb.getcol('WEIGHT')
            nvis   = weight.size/npnt

            w = value * np.sqrt(nvis)
            w = 1.0/(w*w)
            print "%20s      Old=%s New=%g Nvis=%d (per npnt=%d)" \
                % (ims,weight[0,0],w,nvis,npnt)
            weight[:,:] = w           # weight[npol,nvis]

            tb.putcol('WEIGHT',weight)
            tb.close()

        return

    # Size matching
    # -------------
    elif (mode == 4):      # "beam size matching"
        print "TP2VISWT: equalize synthesized and convolution beams"
        print "          TP weights are set w.r.t. interferometer weights"

        # calculate Omega_conv
        # --------------------

        if type(int_beam) == float or type(int_beam) == int:
            print "Assuming spherical int_beam = %g arcsec" % int_beam
        else:
            print "fitting not yet implemented for INT_BEAM"
            print "unexpected int_beam",type(int_beam),int_beam
            return

        bmaj_int = int_beam / apr              # in radians
        bmin_int = int_beam / apr
        omega_conv = np.pi/(4*np.log(2.0)) * bmaj_int*bmin_int

        # Calculate Omega_syn
        # -------------------

        # from math, W_TP(0,0) = 2*pi*sigma^2=2*pi*(FWHM/2sqrt(2ln2))^2
        if type(tp_beam) == float or type(tp_beam) == int:
            print "Assuming spherical tp_beam = %g arcsec" % tp_beam
            bsize1 = bsize2 = tp_beam / apr
            omega_syn = np.pi/(4*np.log(2.0)) * bsize1 * bsize2

        # read dirty beam of TP and derive W_TP(0,0)
        else:
            print "Assuming CASA tp_beam = %s" % tp_beam        
            ia.open(tp_beam)
            h0 = ia.summary()
            naxis1 = h0['shape'][0]
            naxis2 = h0['shape'][1]
            naxis3 = h0['shape'][3]   # RA-DEC-POL-FREQ cube !!!
            cdelt1 = h0['incr'][0]    # radians
            cdelt2 = h0['incr'][0]
            # Use middle channel 
            iz=naxis3
            if naxis3>1: iz = naxis3/2
            image_tp  = ia.getchunk([-1,-1,-1,iz],[-1,-1,-1,iz]).squeeze() # image[ix][iy]   
            ia.close()

            # Calculate W_TP(0,0)
            #   forward fft requires tricky normalization since B(0,0) is
            #   set to =1 already. Instead, we use inverse fft
            image_tp_fft = np.fft.ifft2(image_tp)
            wtp00 = np.max(np.abs(image_tp_fft))
        
            # Calculate Omega_syn [including terms for discrete FT]
            bsize1 = np.abs(naxis1*cdelt1)
            bsize2 = np.abs(naxis2*cdelt2)
        
            omega_syn = bsize1 * bsize2 * wtp00
            print "Derived TP beam:",np.sqrt(omega_syn)*apr

        # Derive beta
        #   omega_syn  = beta/(1+beta)*W_TP(0,0)
        # --------------------------------------
        beta = omega_conv / (omega_syn - omega_conv)

        # Print
        print "Omega_conv [strad] = ",omega_conv
        print "Omega_syn  [strad] = ",omega_syn 
        print "Beta      = ",beta
    
        if beta < 0:
            print "Problem: you cannot have a TP beam that's smaller than an INT beam"
            return

        # Adjust weights of TP visibilities
        #   w_TP,k = beta/Nvis_TP * sum_k(w_INT,k)
        # ----------------------------------------

        if msINT == None:
            print "Need a list of MS in order to set the beam size weights for the msTP"
            return

        # sumup the weights of INT visibilities
        sumw = 0.0
        for ims in msINT:
            tb.open(ims)
            weight_int = tb.getcol('WEIGHT')
            sumw = sumw + weight_int.sum()              # both POLs
            tb.close()
        del weight_int
        print "INT sumw = ",sumw

        # count the number of TP visibilities
        nvis = 0
        for ims in msTP:
            tb.open(ims)
            weight = tb.getcol('WEIGHT')    
            nvis   = nvis + weight.size
            tb.close()
            
        # weight to be set
        w = beta * sumw / nvis
        print "TP nvis,weight=",nvis,w

        # set TP weights with respect to INT weights
        for ims in msTP:
            tb.open(ims,nomodify=False)
            weight = tb.getcol('WEIGHT')    
            weight[:,:] = w
            tb.putcol('WEIGHT',weight)
            tb.close()

        return

    else:

        print "WEIGHT: Unknown mode ",mode

    #-end of tp2viswt()

## =================================
## TP2VISPL: Plot visibility weights
## =================================

def tp2vispl(msTP=None, ms07=None, ms12=None, ampPlot=True, show=False):
    """
    Plotting TP, 7m, and 12m MSs

    Parameters:
    -----------
    # need at least msTP
    msTP      TP measurement set
    ms07      07m measurement set
    ms12      12m measurement set

    ampPlot   True     amp-uvdistance plot
              False    weight-uvdistance plot
    show      True     plot in display as well as in file
              False    not plot in display, but in file
    """

    print "TP2VISPL: Plot MSs - Takes time."

    # Parameters
    # ----------

    # Parameters for plot
    bin     = 0.5   # m
    uvMax   = 150.0 # m
    uvZoom  = 50.0  # m
    outfig  = 'plot_tp2viswt.png'

    # Generate MS list [order: 12m -> 7m -> TP]
    mslist = []
    clist  = [] # color for plot
    if ms12 != None:
        if type(ms12) != type([]): ms12 = [ms12]
        for ims in ms12:
            mslist.append(ims)
            clist.append('b')
    if ms07 != None:
        if type(ms07) != type([]): ms07 = [ms07]
        for ims in ms07:
            mslist.append(ims)
            clist.append('g')
    if msTP != None:
        if type(msTP) != type([]): msTP = [msTP]
        for ims in msTP:
            mslist.append(ims)
            clist.append('r')

    # Find the freq of max flux channel [preferably from msTP]
    msfile = mslist[-1]
    print "Pick up max flux channel in %s" % (msfile)
    ms.open(msfile)
    ms.selectinit(datadescid=0)
    axinfo     = ms.getdata(['axis_info'])['axis_info']
    ampsum     = ms.getdata(['amplitude'])['amplitude'][0,:,:].sum(1)
    imaxchan   = ampsum.argmax(0)
    targetfreq = axinfo['freq_axis']['chan_freq'][imaxchan][0]
    ms.close()

    del ampsum

    # Start figure
    plt.ioff()                    # no interactive mode
    fig  = plt.figure()
    axtl = fig.add_subplot(2,2,1) # Top-left
    axtr = fig.add_subplot(2,2,2) # Top-right
    axbr = fig.add_subplot(2,2,4) # Bot-right
    axbl = fig.add_subplot(2,2,3) # Bot-left

    # Loop over MSs
    # -------------
    for ims in range(len(mslist)):

        # Set up MSs
        msfile = mslist[ims]
        color  = clist[ims]
        
        # Read and calculate parameters
        # -----------------------------

        # Open MS and set constraints
        ms.open(msfile,nomodify=True)
        ms.selectinit(datadescid=0)
        ms.select({'uvdist':[0.0,uvMax]})

        # Read params
        uu       = ms.getdata(['u'])['u'] # meter
        vv       = ms.getdata(['v'])['v'] # meter
        uvd      = ms.getdata(['uvdist'])['uvdist'] # meter
        wt       = ms.getdata(['weight'])['weight'][0,:] # 1st pol

        # Channel freq, width, and amplitude
        axinfo   = ms.getdata(['axis_info'])['axis_info']
        chanfreq = axinfo['freq_axis']['chan_freq']
        chanwidt = axinfo['freq_axis']['resolution']
        if ampPlot: # channel closest to target frequency
            ichan = np.argmin(np.abs(chanfreq - targetfreq))
            ms.selectchannel(1,ichan,1,1) # nchan,start,width,inc
            amp   = ms.getdata(['amplitude'])['amplitude'][0,0,:] # 1st pol
        else:       # middle channel
            ichan = axinfo['freq_axis']['chan_freq'].size / 2
        chanfreq = chanfreq[ichan][0]/1.0e9
        chanwidt = chanwidt[ichan][0]/1.0e9
        ms.close()

        # Weight in 1GHz width
        wt       = wt / np.abs(chanwidt)

        # print
        print "%30s: (chan, freq, fwid, wtmax) = \
            (%5d, %10.6f, %10.6f, %10.6f)" % \
            (msfile,ichan,chanfreq,chanwidt,wt.max())

        # Calculate avarages
        uvbins = np.arange(0.0,uvd.max() + bin, bin)
        digit  = np.digitize(uvd,uvbins)
        uvarea = np.pi*np.diff(uvbins*uvbins)
        wtbins = [wt[digit==i].sum() for i in range(1,len(uvbins))]
        wtbins = wtbins/uvarea
        uvbins = uvbins[1:]

        uvMax = np.maximum(uvMax,uvd.max())

        del digit, uvarea

        # Plot
        # ----

        # Top-left:     UU     vs VV
        axtl.scatter( uu, vv,marker='.',s=0.2,c=color,lw=0)
        axtl.scatter(-uu,-vv,marker='.',s=0.2,c=color,lw=0)
        # Top-right:    UVdist vs Amplitude or Weight [Zoom-up]
        if ampPlot:
            axtr.scatter(uvd,amp,marker='.',s=0.2,c=color,lw=0)
        else:
            axtr.scatter(uvd,wt,marker='.',s=0.2,c=color,lw=0)
        # Bottom-right: UVdist vs WTdens [Zoom-up]
        axbr.plot(uvbins,wtbins,c=color,drawstyle='steps-mid')
        # Bottom-left:  UVdist vs WTdens
        axbl.plot(uvbins,wtbins,c=color,drawstyle='steps-mid',label=ms)

        del uu,vv,uvd,wt,uvbins,wtbins

    # Cosmetics
    # ---------
    fontsize=10

    # Top-left
    axtl.set_aspect(1.0)
    axtl.set_xlim(-uvZoom, uvZoom)
    axtl.set_ylim(-uvZoom, uvZoom)
    axtl.set_xlabel("u (meter)",fontsize=fontsize)
    axtl.set_ylabel("v (meter)",fontsize=fontsize)
    axtl.tick_params(axis='both',which='major',labelsize=fontsize)
    # Top-right
    axtr.set_xlim(0.0, uvZoom)
    axtr.set_xlabel("uvdistance (meter)",fontsize=fontsize)
    if ampPlot:
        axtr.set_ylim(0.0, 20.0)
        axtr.set_ylabel("amplitude",fontsize=fontsize)
    else:
        axtr.set_ylabel("weight [per visibility]",fontsize=fontsize)
    axtr.tick_params(axis='both',which='major',labelsize=fontsize)
    # Bottom-right
    axbr.set_xlim(0.0,uvZoom)
    axbr.set_yscale('log')
    axbr.set_xlabel("uvdistance (meter)",fontsize=fontsize)
    axbr.set_ylabel("weight density [BW=1GHz]",fontsize=fontsize)
    axbr.tick_params(axis='both',which='major',labelsize=fontsize)
    # Bottom-left
    axbl.set_xlim(0.0, uvMax)
    axbl.set_yscale('log')
    axbl.set_xlabel("uvdistance (meter)",fontsize=fontsize)
    axbl.set_ylabel("weight density [BW=1GHz]",fontsize=fontsize)
    axbl.tick_params(axis='both',which='major',labelsize=fontsize)

    # Legend [conflict with the MS functions of CASA Toolkit]
#    axtr.legend(loc=1,prop={'size':'x-small'})

    # Draw
    plt.savefig(outfig)
    print "Output fig in %s" % (outfig)
    if (show):
        plt.show()
    else:
        plt.close(fig)

    #-end of tp2vispl()
