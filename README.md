This repository contains data products from our analysis of the Voyager 1 PWS wideband receiver data spanning the years 2012 through early 2020.

The data products are as follows:

- maskspecv4.npz contains the 2D dynamic spectrum spanning 2017 through 2020. The dynamic spectrum already has interference lines masked and excludes epochs where the rms noise threshold was greater than 0.15 relative intensity units. The file can be loaded in python as follows:

data = load('maskspecv4.npz') # maskspecv4 = ignoring obs w/ rms>0.15 and masking lines
freqprofs = data['spec']
times = data['time']
freqs = data['freq']

- fullspecout2.npz contains the 2D dynamic spectrum spanning 2012 through early 2020. The spectrum can be loaded as follows:

data = load('fullspecout2.npz')
specarray = data['spec']
times = data['time']
freqs = data['freqs']

- The folder Sonification_Spectra contains the frequencies and time-stamps of the weak plasma line. Comments on the format of those files is provided separately iin that folder.
