from numpy import *
from matplotlib.pyplot import *
from datetime import *
from scipy.signal import find_peaks

'''
This is a program to read in high resolution V1 PWS data from the wideband receiver and create spectra.

Usage (command line): python widespec.py filename(s)

Input: V1 PWS high resolution data from the wideband receiver (60 ms sampling time)

Output: quicklook plot and relevant arrays

History:
simple program created -- Stella 05/28/20
copied over to ella for testing -- Stella 06/22/20
'''

def readhighres(fname):

# reads in data and computes the spectrum

    # 1600 samples in time at 28,800 samples per second
    time = arange(0.,1600.,1) / 28.8 # in ms

    # 801 frequency components will be produced by the FFT
    # The highest frequency will be the Nyquist frequency (half the sample rate)
    freq = (arange(0.,801.,1) / 800.) * 14.4 # in kHz

    with open(fname,'rb') as f:

        # 1024-byte records per line
        #rec = bytearray(1024)

        # Read the first (file header) record and extract the PWS ASCII label info
        rec = f.read()
        label = rec[248:298]
        label = label.decode("utf-8")
        print(label)
    
        # we'll use an array of the even indices in order to extract first the high 4-bit values
        # and then the low 4-bit values
        even = arange(0,800,1) * 2

        # Loop through all of the records
        records = arange(0,len(rec),1024)
        fullrecord = []
        for r in records:
    
            # There will be 1600 4-bit waveform samples per record from bytes 220-1019
            samples = zeros(1600)
    
            # Extract the SCLK line count from bytes 23 and 24 (MSB 16-bit integer)
            line = int(rec[r+22]) * 256 + int(rec[r+23])
            #print(rec[r+23], int(rec[r+23]))
            #print(line)

            # The time offset relative to the time in the label, in seconds
            offset = 0.06 * (line - 1)
            #print('offset',offset)
    
            # Extract the high-order 4-bit samples
        
            for j in range(800):
                samples[even[j]] = int(rec[r+j+220]/16)
                samples[even[j]+1] = rec[r+j+220] % 16
            
            fullrecord.append(samples)
        fullrecord = array(fullrecord)
    arrshape = shape(fullrecord)
    nrows = len(fullrecord[:,0])
    fulltime = arange(0,60*nrows,60)*10.**(-3.)
    print('length of data set is ', max(fulltime),' s')

    spectrum = zeros(arrshape)
    for i in range(len(fullrecord[:,0])):
        wfrm = (fullrecord[i,:] - 7.5)*hamming(len(fullrecord[i,:]))
        spectrum[i,:] = (abs(fft.fft(wfrm))**2.)/len(wfrm)
    spectrumswap = swapaxes(spectrum,0,1)

    # make array of sample times for entire raw data time series

    nrows = len(fullrecord[:,0])
    recordtime = zeros(shape(fullrecord))
    for i in range(nrows):
        recordtime[i,:] = time
    # add last element of row + 4.44 ms to all values in next row
    for i in range(nrows-1):
        recordtime[i+1,:] = recordtime[i+1,:] + recordtime[i,-1] + 4.44

    print(log10(min(spectrumswap.flatten())),log10(max(spectrumswap.flatten())))
    
    return label,recordtime,freq,fullrecord,fulltime,spectrumswap

def specavg(spec,freq):
    
    mini,maxi = 2,int(shape(spec)[0]/2)
    timeavg = mean(spec[mini:maxi,:],axis=1)
    freqavg = mean(spec[mini:maxi,:],axis=0)
    
    cut1 = where(freq>1.)[0][0]
    cut2 = where(freq>3.)[0][0]
    cut3 = where(freq>5.)[0][0]
    
    freqavg_cut1 = mean(spec[cut1:maxi,:],axis=0)
    freqavg_cut2 = mean(spec[cut2:maxi,:],axis=0)
    freqavg_cut3 = mean(spec[cut3:maxi,:],axis=0)
    
    return timeavg,freqavg,freqavg_cut1,freqavg_cut2,freqavg_cut3

def specstats(spec,freq,label):
    
    # calculate rms and variance of spectrum (excluding spectral lines)
    
    mini,maxi = 2,int(shape(spec)[0]/2)
    spec = spec[mini:maxi,:]
    
    specflat = spec.flatten()
    n = 1/len(specflat)
    sumsqr = sum(specflat**2.)
    m2 = sqrt(n*sumsqr)
    vari = m2**2. - (n*sum(specflat))**2.
    rms = sqrt(vari)
    
    for i in range(20):
    
        sumsqr = sum(specflat**2.)
        
        clipinds = where(spec > 3*m2)
        sumclip = sum(spec[clipinds]**2.)
        n = 1./(len(specflat)-len(spec[clipinds]))
        meanspec = n*(sum(specflat) - sum(spec[clipinds]))
        m2 = sqrt(n*(sumsqr - sumclip))
        
        vari = m2**2. - meanspec**2.
        rmsfull = sqrt(vari)
        clipspec = delete(spec,clipinds)
        #print(rmsfull)
        
    # calculate rms of spectrum in frequency bins
    
    bins = arange(0.5,14.5,1.)#[1.,3.,5.,7.,9.,11.,13.] #center of each freq. bin
    print(bins)
    bin_rms = []
    bin_mean = []
    
    for b in bins:
    
        l = where((freq<(b+0.5))&(freq>(b-0.5)))[0]
        #print(freq[l])
        spec2 = spec[l]
        specflat = spec2.flatten()
        
        n = 1/len(specflat)
        sumsqr = sum(specflat**2.)
        m2 = sqrt(n*sumsqr)
        vari = m2**2. - (n*sum(specflat))**2.
        rms = sqrt(vari)
        
        for i in range(20):
        
            sumsqr = sum(specflat**2.)
            
            clipinds = where(spec2 > 3*m2)
            sumclip = sum(spec2[clipinds]**2.)
            n = 1./(len(specflat)-len(spec2[clipinds]))
            meanbin = n*(sum(specflat) - sum(spec2[clipinds]))
            m2 = sqrt(n*(sumsqr - sumclip))
            
            vari = m2**2. - meanbin**2.
            rmsbin = sqrt(vari)
            #print(b, rms)
            
        #print(b, rms, meanspec)
        bin_rms.append(rmsbin)
        bin_mean.append(meanbin)
    
    bins = array(bins)
    bin_rms = array(bin_rms)
    bin_mean = array(bin_mean)
    
    plot(bins,bin_rms,label='rms',linewidth=1)
    plot(bins,bin_mean,label='mean',linewidth=1)
    xlabel('Center Frequency [kHz]')
    ylabel('Relative Power')
    title(label)
    legend()
    show()
        
    return rmsfull, meanspec, clipspec

def ACFcalc(spec,meanspec,T,F):

    # calculate ACF of spectrum

    freqlen = int((shape(spec)[0])/2)+1
    timelen = shape(spec)[1]
    bigspec = spec[0:freqlen,:] - meanspec

    preacf = absolute(fft.fft2(bigspec))**2.
    postacf = fft.ifft2(preacf)
    acf = real(fft.fftshift(postacf))

    one = ones(shape(spec[0:freqlen,:]))
    prenorm = absolute(fft.fft2(one))**2.
    postnorm = fft.ifft2(prenorm)
    norm = real(fft.fftshift(postnorm))

    acf = divide(acf,norm)

    ind1,ind2=int(shape(acf)[0]/2),int(shape(acf)[1]/2)
    freqslice = acf[ind1,:]
    timeslice = acf[:,ind2]

    acfT = concatenate((-1*T[:int(timelen/2)+1][::-1],T[1:int(timelen/2)+1]))
    acfF = concatenate((-1*F[:int(freqlen/2)+1][::-1],F[1:int(freqlen/2)+1]))

    return acf,acfT,acfF,timeslice,freqslice

def speclines(timevg,meanspec,rms,freq):
    
   # find frequencies and average powers of spectral lines
    mini,maxi = 2,int(shape(timevg)[0])
    freq = freq[mini:maxi]
    freqinds = where(freq > 1.) # cut out noise below 1 kHz
    freq = freq[freqinds]
    freq1d = timevg[freqinds]
    
    #linfitt = polyfit(freq,freq1d,3)
    #p = poly1d(linfitt)
    #noramp = freq1d - p(freq)
    #plot(freq,freq1d)
    #plot(freq,p(freq))
    #yscale('log')
    
    thresh = (freq1d - meanspec)/rms
    #peaks,heights = find_peaks(freq1d,prominence=1)
    peaks=[]
    for i in range(len(freq1d)):
        if thresh[i]>1.5:
            peaks.append(i)
    plot(freq,freq1d)
    plot(freq[peaks],freq1d[peaks],'o')
    #plot(freq[peaks],thresh[peaks],'o')
    xscale('log')
    yscale('log')
    xlabel('Frequency [kHz]')
    ylabel('Peak Height (# rms units away from mean)')
    show()
    return freq[peaks],freq1d[peaks]

def highpass(fullrecord):
    fc = 0.03  # Cutoff frequency as a fraction of the sampling rate (in (0, 0.5)). sampling rate is 28.8 kHz
    b = 0.007  # Transition band, as a fraction of the sampling rate (in (0, 0.5)).
    N = int(ceil((4 / b)))
    if not N % 2: N += 1  # Make sure that N is odd.
    n = arange(N)
     
    # Compute a low-pass filter.
    h = sinc(2 * fc * (n - (N - 1) / 2))
    w = blackman(N)
    h = h * w
    h = h / sum(h)
     
    # Create a high-pass filter from the low-pass filter through spectral inversion.
    h = -h
    h[(N - 1) // 2] += 1

    filt = []
    for i in range(len(fullrecord[:,0])):
        s = fullrecord[i,:]
        filt.append(convolve(s,h))
    filt = array(filt)
    
    return filt[:,285:1885]

def plots(spectrum,label,recordtime,fullrecord,fulltime,freq,timeavg1,freqavg1,acf,acfT,acfF,timeslice,freqslice,rms,meanspec,clipspec,filtered,freqavg_cut1,freqavg_cut2,freqavg_cut3):
    
    # plot the time series, power spectrum, & ACF

    #set indices for plotting the spectrum
    #min,max = 0,-1
    mini,maxi = 2,int(shape(spectrum)[0]/2) #plots half the frequencies

    fig = figure(figsize=(15,7))
    fig.suptitle(label,fontsize=12)

    ax1 = fig.add_axes([0.06,0.77,0.82,0.15])
    ax1.plot(10.**(-3.)*recordtime.flatten(),fullrecord.flatten()-7.5,linewidth=0.5)
    #ax1.plot(10.**(-3.)*recordtime.flatten(),filtered.flatten(),linewidth=0.5)
    ax1.set_xlim(0,recordtime[-1,-1]*10**(-3))
    ax1.set_xlabel('Time [s]',fontsize=8)
    ax1.set_ylabel('Raw Amplitude',fontsize=8)
    ax1.set_title('Time Series of Raw Wideband Data',fontsize=8)
    ax1.tick_params(labelsize='small')
    
    ax9 = fig.add_axes([0.36,0.6,0.29,0.1])
    ax9.plot(10.**(-3.)*recordtime[2,:],fullrecord[2,:]-7.5,linewidth=1)
    #ax9.plot(10.**(-3.)*recordtime[2,:],filtered[2,:],linewidth=1)
    ax9.set_xlim(min(10.**(-3.)*recordtime[2,:]),max(10.**(-3.)*recordtime[2,:]))
    ax9.set_xlabel('Time [s]',fontsize=8)
    ax9.set_ylabel('Raw Amplitude',fontsize=8)
    ax9.tick_params(labelsize='small')
    
    ax8 = fig.add_axes([0.885,0.77,0.09,0.15])
    ax8.hist(fullrecord.flatten()-7.5,orientation='horizontal',bins=15,color='green',alpha=0.6)
    ax8.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    ax8.xaxis.offsetText.set_fontsize(8)
    ax8.set_yticklabels([])
    ax8.set_xlabel('# Samples',fontsize=8)
    ax8.tick_params(labelsize='small')

    ax2 = fig.add_axes([0.31,0.18,0.08,0.38])
    ax2.plot(timeavg1,freq[mini:maxi],linewidth=1)
    ax2.set_xscale('log')
    ax2.set_xticks([0.1,1,10])
    ax2.set_yticks([])
    ax2.set_xlabel('Avg. Power',fontsize=8)
    ax2.set_ylim(freq[mini],freq[maxi])
    ax2.tick_params(labelsize='small')

    ax3 = fig.add_axes([0.06,0.18,0.25,0.38])
    specim = ax3.imshow(log10(spectrum[mini:maxi,:]),aspect='auto',origin='lower',extent=[fulltime[0],fulltime[-1],freq[mini],freq[maxi]],vmin=-6.5,vmax=3.)
    ax3.set_xlabel('Time [s]',fontsize=8)
    ax3.set_ylabel('Frequency [kHz]',fontsize=8)
    ax3.tick_params(labelsize='small')
    cax = fig.add_axes([0.06,0.08,0.25,0.03])
    cbar = colorbar(specim,cax=cax,orientation='horizontal')
    cbar.set_label('log$_{10}$(Relative Power)',fontsize=8)
    cbar.ax.tick_params(labelsize='small')

    ax4 = fig.add_axes([0.06,0.56,0.25,0.12])
    ax4.plot(fulltime,freqavg1,linewidth=1,label='all freqs')
    ax4.plot(fulltime,freqavg_cut1,linewidth=1,label='> 1 kHz',color='green')
    ax4.plot(fulltime,freqavg_cut2,linewidth=1,label='> 3 kHz',color='orange')
    ax4.plot(fulltime,freqavg_cut3,linewidth=1,label='> 5 kHz',color='purple')
    ax4.legend(fontsize=6,ncol=4,loc='lower left',frameon=False,bbox_to_anchor=(0.03, 1.))
    ax4.set_xticks([])
    ax4.set_yscale('log')
    ax4.set_yticks([0.1,1])
    ax4.set_xlim(fulltime[0],fulltime[-1])
    ax4.set_ylabel('Avg. Power',fontsize=8)
    ax4.tick_params(labelsize='small')
    
    ax10 = fig.add_axes([0.44,0.18,0.14,0.31])
    ax10.hist(clipspec,bins=40,color='green',alpha=0.6)
    ax10.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
    ax10.xaxis.offsetText.set_fontsize(8)
    ax10.set_xlabel('Rel. Spectral Power',fontsize=8)
    ax10.set_ylabel('# Samples',fontsize=8)
    ax10.set_yscale('log')
    ax10.tick_params(labelsize='small')

    ax11 = fig.add_axes([0.41,0.08,0.1,0.1])
    ax11.text(0,0,'spectral rms: '+str('%.3f'%(rms))+' | mean: '+str('%.3f'%(meanspec)),fontsize=10)
    ax11.text(0,-0.3,r'time duration: '+ str('%.3f'%(max(fulltime)))+' s',fontsize=10)
    ax11.axis('off')

    ax5 = fig.add_axes([0.69,0.18,0.25,0.38])
    acfim = ax5.imshow(log10(acf),aspect='auto',origin='lower',extent=[acfT[0],acfT[-1],acfF[0],acfF[-1]],vmin=-7.,vmax=3.)
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position('right')
    ax5.set_xlabel('Time Lag [s]',fontsize=8)
    ax5.set_ylabel('Frequency Lag [kHz]',fontsize=8)
    ax5.tick_params(labelsize='small')
    cax = fig.add_axes([0.69,0.08,0.25,0.03])
    cbar = colorbar(acfim,cax=cax,orientation='horizontal')
    cbar.set_label(r'log$_{10}$(ACF)',fontsize=8)
    cbar.ax.tick_params(labelsize='small')
        
    ax6 = fig.add_axes([0.61,0.18,0.08,0.38])
    y = linspace(0,len(timeslice),len(timeslice))
    ax6.plot(timeslice,y,drawstyle='steps',label='0 s',linewidth=1)
    ax6.set_xlim(max(timeslice)+1000,10.**(-4))
    ax6.set_xscale('log')
    ax6.set_xticks([0.01,1,100])
    ax6.set_ylim(min(y),max(y))
    ax6.set_yticks([])
    ax6.tick_params(labelsize='small')
    ax6.legend(fontsize=6)

    ax7 = fig.add_axes([0.69,0.56,0.25,0.12])
    ax7.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    x = linspace(0,len(freqslice),len(freqslice))
    ax7.plot(x,freqslice,label='0 kHz',drawstyle='steps',linewidth=1)
    #ax7.set_yscale('log')
    #ax7.set_yticks([10**3])
    ax7.set_xlim(min(x),max(x))
    #ax7.plot(acf[700,:],label=str(acfF[700])+' kHz')
    #ax7.plot(acf[300,:],label='-8.98 kHz')
    ax7.legend(fontsize=6)
    ax7.set_xticks([])
    ax7.tick_params(labelsize='small')

    show()
    #savefig('/home/ella1/stella/voyager/wideband_out/quicklook/figs'+str(label),format='png')

inpath = '/Users/stellaocker/Research/NASA-OH-GI/highres-calibration/C2640311.DAT'
#inpath = '/Users/stellaocker/Research/NASA-OH-GI/binned_rms/C0982711.DAT'
#inpath = '/Users/stellaocker/Research/NASA-OH-GI/binned_rms/C2121358.DAT'
#inpath = '/Users/stellaocker/Research/NASA-OH-GI/binned_rms/C4642358.DAT'

def main(file):

    # ultimately will need to loop through all files given on command line & save arrays + figure to separate directories

    label,recordtime,freq,fullrecord,fulltime,spectrum = readhighres(file)
    rms,meanspec,clipspec = specstats(spectrum,freq,label)
    acf,acfT,acfF,timeslice,freqslice = ACFcalc(spectrum,meanspec,fulltime,freq)
    timeavg1,freqavg1,freqavg_cut1,freqavg_cut2,freqavg_cut3 = specavg(spectrum,freq)
    #filtered = highpass(fullrecord)
    speclines(timeavg1,meanspec,rms,freq)
    
    # plot ACF of time series -- probably not quite right b/c of uneven sampling
    #fr = fullrecord.flatten() - 7.5
    #ACFT = fftconvolve(fr,fr[::-1],mode='full')
    #t = 10.**(-3)*recordtime.flatten()
    #acftime = concatenate((-1*t[:-1][::-1], t))
    #plot(acftime,ACFT)
    #xlabel('Time Lag (s)')
    #title('ACF of Wideband Time Series')
    #show()
    
    #plots(spectrum,label,recordtime,fullrecord,fulltime,freq,timeavg1,freqavg1,acf,acfT,acfF,timeslice,freqslice,rms,meanspec,clipspec,filtered,freqavg_cut1,freqavg_cut2,freqavg_cut3)
    
    #savez('/home/ella1/stella/voyager/wideband_out/quicklook/arrays'+str(label),spectrum,recordtime,fullrecord,fulltime,freq,acf)


main(inpath)
