import numpy
from pycbc.filter import lowpass_fir, highpass_fir
from pycbc.frame import read_frame
from pycbc.filter import highpass
from pycbc.types import load_timeseries, TimeSeries
from pycbc.filter import resample_to_delta_t
from scipy.signal import filtfilt, butter, iirdesign, zpk2tf, freqz

# The irr_bandstops, get_filter_coefs, and filter_data functions are taken from
# https://www.gw-openscience.org/s/events/GW150914/GW150914_tutorial.html?cm_mc_uid=68374226047614580312966&cm_mc_sid_50200000=1458031296
# with some minr modifications (which are noted in the functions).
# The ``notches`` list is copied from the ``notchesAbsolute`` in the 
# ``get_filter_coefs`` on the GW OSC site.
notches = [14.00, 34.70, 35.30, 35.90, 36.70, 37.30, 40.95, 60.00, 120.00, 179.99, 304.99, 331.49,
510.02, 1009.99]
order = 4096
ifos = ['H1', 'L1']

# Consistent with LIGO reported results we use 7 ms +- 0.5 ms. 
time_diff =  0.007
time_uncertainty = 0.0005

# The event_time is the UTC time of the detection in GPS seconds, rounded down
# to the nearest second. This corresponds to 14 Sept 2015 09:50:45 UTC, which
# comes from Abbott et al. PRL 116, 061102 (2016) (the GW150914 PRL).
event_time = 1126259462
# The corr_time is the start time of the window that Creswell et al. used for
# their cross correlation statistic. Details can be found here:
# http://www.nbi.ku.dk/gravitational-waves/gravitational-waves.html
# Note that in the plots and notebooks linked from that page, the cross
# correlation window starts at 0.39. The time axis in those plots comes from
# the Fig. 1 data of the GW150914 PRL. In the caption for that plot, it states
# that the times are relative to the GW150914 UTC time, which here is
# event_time. Thus, in GPS seconds, the window used by Creswell et al. starts
# at event_time + 0.39.
corr_time = 1126259462.39

def iir_bandstops(fstops, fs, order=4):
    """ellip notch filter
    fstops is a list of entries of the form [frequency (Hz), df, df2]                           
    where df is the pass width and df2 is the stop width (narrower                              
    than the pass width). Use caution if passing more than one freq at a time,                  
    because the filter response might behave in ways you don't expect.
    """
    nyq = 0.5 * fs

    # Zeros zd, poles pd, and gain kd for the digital filter
    zd = numpy.array([])
    pd = numpy.array([])
    kd = 1

    # Notches
    for fstopData in fstops:
        fstop = fstopData[0]
        df = fstopData[1]
        df2 = fstopData[2]
        low = (fstop - df) / nyq
        high = (fstop + df) / nyq
        low2 = (fstop - df2) / nyq
        high2 = (fstop + df2) / nyq
        z, p, k = iirdesign([low,high], [low2,high2], gpass=1, gstop=6,
                            ftype='ellip', output='zpk')
        zd = numpy.append(zd,z)
        pd = numpy.append(pd,p)

    # Set gain to one at 100 Hz...better not notch there                                        
    bPrelim,aPrelim = zpk2tf(zd, pd, 1)
    outFreq, outg0 = freqz(bPrelim, aPrelim, 100/nyq)

    # Return the numerator and denominator of the digital filter                                
    b,a = zpk2tf(zd,pd,k)
    return b, a

def get_filter_coefs(fs, lowcut, highcut):
    """This is the same as the funciton on the GW OSC site, but we make
    lowcut and highcut input arguments (on the GW OSC site they are hard coded
    to 43 and 260); these are the frequencies (in Hz) at which to apply the
    bandpass filter.
    """
    
    # assemble the filter b,a coefficients:
    coefs = []

    # bandpass filter parameters
    order = 4

    # bandpass filter coefficients 
    nyq = 0.5*fs
    low = lowcut / nyq
    high = highcut / nyq
    bb, ab = butter(order, [low, high], btype='band')
    coefs.append((bb,ab))

    # Frequencies of notches at known instrumental spectral line frequencies.
    # You can see these lines in the ASD above, so it is straightforward to make this list.
    notchesAbsolute = numpy.array(notches)

    # notch filter coefficients:
    for notchf in notchesAbsolute:                      
        bn, an = iir_bandstops(numpy.array([[notchf,1,0.1]]), fs, order=4)
        coefs.append((bn,an))

    # Manually do a wider notch filter around 510 Hz etc.          
    bn, an = iir_bandstops(numpy.array([[510,200,20]]), fs, order=4)
    coefs.append((bn, an))

    # also notch out the forest of lines around 331.5 Hz
    bn, an = iir_bandstops(numpy.array([[331.5,10,1]]), fs, order=4)
    coefs.append((bn, an))
    
    return coefs

# and then define the filter function:
def filter_data(data_in, coefs):
    data = data_in.copy()
    for coef in coefs:
        b,a = coef
        # filtfilt applies a linear filter twice, once forward and once backwards.
        # The combined filter has linear phase.
        data = filtfilt(b, a, data)
    return TimeSeries(data, delta_t=data_in.delta_t, epoch=data_in.start_time)

def get_residual_strain(window=512):
    strain = {}
    for ifo in ifos:
        strain[ifo] = load_timeseries('residuals.hdf', group='%s/residual' % ifo)
        strain[ifo] = strain[ifo].time_slice(event_time-window/2, event_time+window/2)
    return strain

def get_raw_strain(window=512):
    """This function loads the entire frame file downloaded from GW OSC,
    applies a high-pass filter at 15 Hz, then keeps the data around
    ``event_time`` +/-``window``. The ``event_time`` is from above; the
    default window is 512 (s).
    """
    strain = {}
    for ifo in ifos:
        fname = '%s-%s_LOSC_4_V2-1126257414-4096.gwf' % (ifo[0], ifo)
        channel_name = '%s:LOSC-STRAIN' % ifo
        strain[ifo] = read_frame(fname, channel_name)
        strain[ifo] = highpass(strain[ifo], 15.0)
        strain[ifo] = strain[ifo].time_slice(event_time-window/2, event_time+window/2)
    return strain

def get_nrsub_strain(strain):
    nr = {}
    for ifo in ifos:
        nr[ifo] = load_timeseries('fig1-waveform-%s.txt' % ifo[0]) * 1e-21
        nr[ifo] = resample_to_delta_t(nr[ifo], strain[ifo].delta_t)
        nr[ifo].start_time += event_time
        
    ts2 = {}
    for ifo in ifos:
        ts2[ifo] = strain[ifo].copy()
        spart = ts2[ifo].time_slice(nr[ifo].start_time, nr[ifo].end_time)    
        nr[ifo].start_time = spart.start_time
        spart -= nr[ifo] 
    return ts2

def get_fig1_observed():
    strain = {}
    for ifo in ifos:
        strain[ifo] = load_timeseries('fig1-observed-%s.txt' % ifo[0]) * 1e-21
    return strain

def get_fig1_res():
    strain = {}
    for ifo in ifos:
        strain[ifo] = load_timeseries('fig1-residual-%s.txt' % ifo[0]) * 1e-21
    return strain

def bandpass(strain, flow=35.0, fhigh=350.0, bw=3, method='losc'):
    ts = {}
    if method=='fir':
        for ifo in ifos:
            ts[ifo] = strain[ifo].lowpass_fir(fhigh, order).highpass_fir(flow, order)

            for f in notches:
                ts[ifo] = ts[ifo].notch_fir(f-bw/2.0, f+bw/2.0, order)

    if method=='losc':
        # get filter coefficients
        coefs = get_filter_coefs(strain['H1'].sample_rate, flow, fhigh)
        for ifo in ifos:
            ts[ifo] = filter_data(strain[ifo], coefs)

    return ts

def whiten(strain, flow=35.0, fhigh=350.0):
    ts = {}
    for ifo in ifos:
        ts[ifo] = strain[ifo].crop(8, 8).whiten(32, 32, trunc_method=None,
                                                        low_frequency_cutoff=15.0)
        ts[ifo] =  ts[ifo].lowpass_fir(fhigh, order).highpass_fir(flow, order)  
    return ts    


def cross_correlation(h1, l1, time, w=0.04):
    """This is the cross correlation function used by Creswell et al. We copied
    it directly from here:

    http://www.nbi.ku.dk/gravitational-waves/correlations.html

    (That notebook is linked from:
    http://www.nbi.ku.dk/gravitational-waves/gravitational-waves.html)

    The only changes we make here are to accomodate our data types: in that
    notebook the strain and times come from 2D numpy arrays that were loaded
    from text files. Here, ``h1`` and ``l1`` are pycbc.TimeSeries types. We
    note these differences below.
    """
    # The time steps can be retrieved from the TimeSeries.sample_times
    # attribute, which is a pycbc.types.Array; the `.numpy()` converts that
    # to a numpy.array.
    h1_res_times = h1.sample_times.numpy()
    h1_res_strain = h1.numpy()

    l1_res_times = l1.sample_times.numpy()
    l1_res_strain = l1.numpy()

    fs = 1./(h1_res_times[1]-h1_res_times[0])

    min_indxt = numpy.where(abs(h1_res_times-time) == abs(h1_res_times-time).min())[0][0]
    max_indxt = numpy.where(abs(h1_res_times-(time+w)) == abs(h1_res_times-(time+w)).min())[0][0]

    deltatau = 0.01

    tauind_min = int(-deltatau*fs); tauind_max = int(+deltatau*fs)
    tauind = numpy.arange(tauind_min, tauind_max)
    tau = tauind / fs

    corr = []

    for i in tauind:
        corr.append(numpy.corrcoef(h1_res_strain[min_indxt+abs(tauind_min)+i:max_indxt-abs(tauind_max)+i],l1_res_strain[min_indxt+abs(tauind_min):max_indxt-abs(tauind_max)])[0][1])

    corr = numpy.array(corr)
    return tau, TimeSeries(corr, epoch=tauind_min*h1.delta_t, delta_t=h1.delta_t)


def indices_within_window(corr, width=time_uncertainty):
    # The start/end points correspond to rounding to the nearest sample in the
    # outwards direction so we inlcude the entire time range.
    times = corr.sample_times
    l = numpy.searchsorted(times, time_diff - width, side='right') - 1
    r = numpy.searchsorted(times, time_diff + width, side='left') + 1
    return l, r, times[l], times[r-1] 


def corr_near_ml(corr, width=time_uncertainty):
    """Get the value of the peak correlation within the expected time offset
    region.
    """
    l, r, start, end = indices_within_window(corr, width=width)
    return abs(corr[l:r]).max()

def get_peak(strain, time, width=time_uncertainty):
    """Given two time series, cross correlate them, and return the peak
    correlation within the standard time offset.

    Parameters
    ----------
    strain : dict
        Dictionary mapping H1, L1 -> TimeSeries.
    time : float
        The time to analyze.
    width : float, optional
        The size of the time window around time to get the peak. Default is
        ``time_uncertainty``.
    """
    ts = {}
    for ifo in ifos:
        ts[ifo] = strain[ifo].time_slice(time-3, time+3)
    tau, corr = cross_correlation(ts['H1'], ts['L1'], time) 
    peak = corr_near_ml(corr, width=width)
    return peak

