# tp2vis
Total Power to Visibilities (TP2VIS): an ALMA Development Study (2016/17)

This project provides tools to create visibilities from a single dish cube using the method
of [Koda et al. 2011](http://adsabs.harvard.edu/abs/2011ApJS..193...19K).
The TP visibilities can then be combined with the interferometric visibilities in
a joint deconvolution using for example CASA's
[**tclean()**](https://casa.nrao.edu/casadocs/casa-5.0.0/global-task-list/task_tclean/about)
method.

Our github repo for distribution : https://github.com/tp2vis/distribute

Efforts behind the scenes will be released : https://github.com/kodajn/tp2vis (currently not shown yet)

## Release Schedule

The current beta release is for experts with experience on interferometer data reduction with CASA. We seek feedback from the experts for further development. Release to non-experts is planned in future.


## Usage

To give you a quick idea how to run TP2VIS, here are the basic commands in CASA, broken into 5 pieces. Examples are listed below.


### 1. Preparations:

#### 1.1: Make a pointing (**ptg**) file

       listobs('calibrated_12m.ms',listfile='calibrated_12m.log')

after which you can cut and paste the appropriate **Fields**
in the following format into a pointing file, **12m.ptg**, which has to be a simple text file:

       cat 12m.ptg
       
       J2000 05h39m45.660s -70d07m57.524s
       J2000 05h39m54.340s -70d07m57.524s
       J2000 05h39m41.320s -70d07m19.175s
       ...

This can be a little cumbersome, so in the examples listed below you can see examples using awk and grep.

#### 1.2: Find a reasonable RMS

We need to know the RMS in the TP cube from some line free channels. For example, you might be able
to use the first 10 channels of your TP cube

       imstat('tp.im',axes=[0,1])['rms'][:10].mean()
       -> 0.67

#### 1.3: Cutting down the MS dataset sizes.

CASA tasks such as
[**split()**](https://casa.nrao.edu/casadocs/casa-5.0.0/global-task-list/task_split/about)
and
[**mstransform()**](https://casa.nrao.edu/casadocs/casa-5.0.0/global-task-list/task_mstransform/about)
can be very useful to cut down the size
of a measurement set before using it in tclean, including a spectral gridding.
There are a number of examples in the workflows linked below.


### 2. Load and run TP2VIS:

       execfile('tp2vis.py')                                          # load tp2vis 

       tp2vis('tp.im','tp.ms','12m.ptg',rms=0.67)                     # make visibilities

### 3. Expert mode Weighting Schemes (optional)

       tp2viswt('tp.ms', ...)                                         # set weights
       tp2vispl(msTP='tp.ms')                                         # (optional) plot weights

A practical weighting scheme is under investigation. For now "tp2viswt" provides functions to manipulate weights. There are a number of modes how you can set weight, described in [example1](example1.md).

### 4. Some CASA workarounds to get files on a common velocity grid

Currently we are suffering from CASA crashing when **tclean()** uses a list of MS files that included
the TPMS, so we need to
[**concat()**](https://casa.nrao.edu/casadocs/casa-5.0.0/global-task-list/task_concat/about)
them.

       concat(vis=['12m.ms','7m.ms','tp.ms'], concatvis='all.ms',copypointing=False)  # also ignore POINTING

The "copypointing=False" option is important.

### 5. finally the joint deconvolution using CASA's ``tclean()``

       tclean(vis='all.ms', imagename='all', ...)                     # run clean the way you like it

Where ... represents the large number of options to control the deconvolution. For example, users may try "robust" and "uvtaper" options.


## Examples

* [example1:](example1.md)  M100 (data from CASA guide on [**Feather**](https://casaguides.nrao.edu/index.php/M100_Band3_Combine_4.3))

## References

* Koda et al. 2011, ApJS, 193, 19 : http://adsabs.harvard.edu/abs/2011ApJS..193...19K

* CASA reference manual and cookbook : http://casa.nrao.edu/docs/cookbook/
  * Measurment Set: https://casa.nrao.edu/casadocs/casa-5.0.0/reference-material/measurement-set
  * MS V2 document: [MS v2 memo](https://casa.nrao.edu/casadocs/casa-5.0.0/reference-material/229-1.ps/@@download/file/229.ps)
* CASA simulations: https://casa.nrao.edu/casadocs/casa-5.0.0/simulation
  * Simulations (in 4.4) https://casaguides.nrao.edu/index.php/Simulating_Observations_in_CASA_4.4
  * See also our workflow4 below
* CASA single dish imaging:  https://casa.nrao.edu/casadocs/casa-5.0.0/single-dish-imaging
* CASA feather: https://casa.nrao.edu/casadocs/casa-5.0.0/image-combination/feather
* Nordic Tools SD2VIS: https://www.oso.nordic-alma.se/software-tools.php
* CASA data weights and combination:  https://casaguides.nrao.edu/index.php/DataWeightsAndCombination
* Kauffman's *Adding Zero-Spacing* workflow: https://sites.google.com/site/jenskauffmann/research-notes/adding-zero-spa
