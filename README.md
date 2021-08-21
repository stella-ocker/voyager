This repository contains data products from our analysis of the Voyager 1 PWS wideband receiver data spanning the years 2012 through early 2020 as published in doi:10.1038/s41550-021-01363-7 (arXiv:2105.04000). The raw PWS data can be downloaded from https://space.physics.uiowa.edu/voyager/data/. 

The data products are as follows:

- maskspecv4.npz contains the 2D dynamic spectrum spanning 2017 through 2020. The dynamic spectrum already has interference lines masked and excludes epochs where the rms noise threshold was greater than 0.15 relative intensity units. The file can be loaded in python as follows:

  data = load('maskspecv4.npz')
  
  freqprofs = data['spec']
  
  times = data['time']
  
  freqs = data['freq']

- fullspecout2.npz contains the 2D dynamic spectrum spanning 2012 through early 2020. The spectrum can be loaded as follows:

  data = load('fullspecout2.npz')
  
  specarray = data['spec']
  
  times = data['time']
  
  freqs = data['freqs']

- The folder Sonification_Spectra contains the frequencies and time-stamps of the weak plasma line. Comments on the format of those files is provided separately in that folder.

A program widespec.py is also provided to read in raw data files from the PWS wideband receiver. 
