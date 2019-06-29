import numpy as np
from scipy.io import wavfile

deg = np.pi / 180
sclip = lambda x: x/(1 + np.abs(x))

# here lies the fudge factory, where the fudge factors are made (with love!)
FF_AMPLITUDE = 1/4
FF_H_PRESCALE = 2
FF_DURATION = 60
FF_DUR_POWER = 1/2
FF_MIDPT_DRIVE = 2
FF_LOW_DRIVE = 4
FF_LOW_BOOST = 1

full_length = 3600                 # [s]
abridgement_factor = 10
wave_length = full_length//abridgement_factor
wave_SR = 48000                    # sampling rate to render [Hz]
wave_smp = wave_SR*wave_length
nep_freq = 1440                     # [Hz]
xf_halflives = 10   # how many halflives of decay of each object's bleep to render. 10 corresponds roughly to -60 dB?


def bleep(f, T, a, sr=48000):
    return a*np.sin(2*np.pi*f*np.linspace(0, T, int(T*sr), endpoint=False))


def pan_mono(y, rightness=0):
    """ pan_mono: Given a mono signal, i.e. 1-d array, and a pan position, return a panned stereo signal.

    :param y:
    :param float rightness: pan position from -1 (hard left) to 1 (hard right). Outside this range clips.
    :return: stereo
    """

    assert y.ndim == 1
#    stereo = np.zeros((2, y.shape[0]))
    if rightness < -1:
        shh = np.zeros_like(y)
        stereo = np.vstack((y*(-rightness), shh))
    elif rightness > 1:
        shh = np.zeros_like(y)
        stereo = np.vstack((shh, y*rightness))
    else:
        L = (1-rightness)/2
        R = (1+rightness)/2
        stereo = np.vstack((L*y, R*y))
    return stereo


def exp_fade(y, midpt=0.5, halflives=10):
    # halflives = 10 => envelope minimum at 2**-10 ≈ -60dB
    # halflives = 16 => envelope minimum at 2**-16 ≈ -96dB
    N = len(y)

    N_in = int(midpt * N)
    N_out = N - N_in

    in_win = np.power(2, np.linspace(-halflives, 0, N_in, endpoint=False))
    y[:N_in] *= in_win
    if N_out == N_in:
        y[-N_out:] *= in_win[::-1]
    elif N_out > 0:
        out_win = np.power(2, np.linspace(0, -halflives, N_out, endpoint=False))
        y[-N_out:] *= out_win

    return


def add_wrap_st(cumulant, summand, pos, offset=0.5):
    """
    add_wrap_st: add a stereo (i.e. shape (2,n)) array elementwise to a longer one, in place, wrapping around if necessary
        Note: wraps at most once. If the remaining head or tail is still longer than cumulant, the function will fail)
        Use multiwrap_st instead.

    :param ndarray cumulant: the longer array
    :param ndarray summand: the shorter array
    :param pos: the index of cumulant with which to align summand
    :param offset: the relative position within summand to line up with pos
            the default, 0.5, places summand's center at pos
            0.0 should place the 0th element of summand at pos
            1.0 should place the last element of summand at pos
    :return: tuple (prewrap, postwrap): the number of samples that were wrapped around from beginning to end or end to beginning, respectively
    """

    pos = int(np.rint(pos))  # for safety's sake

    c_len = cumulant.shape[1]
    s_len = summand.shape[1]

    s_pos = int(np.around(s_len * offset))
    s_pre = s_pos
    s_post = s_len - s_pos - 1

    prewrap = 0
    postwrap = 0

    if pos - s_pre < 0:
        prewrap = s_pre - pos
        #        print('pre{}'.format(prewrap))
        # wrap around to the end
        cumulant[:, -prewrap:] += summand[:, :prewrap]
        s_pre -= prewrap
    if pos + s_post > (c_len - 1):
        # wrap around to the start
        postwrap = pos + s_post - (c_len - 1)
        #        print('post{}'.format(prewrap))
        cumulant[:, :postwrap] += summand[:, -postwrap:]
        s_post -= postwrap

    #    print('s_pre{}s_post{}'.format(s_pre,s_post))
    cumulant[:, (pos - s_pre):(pos + s_post + 1)] += summand[:, (s_pos - s_pre):(s_pos + s_post + 1)]
    return prewrap, postwrap


def multiwrap_st(cumulant, summand, pos, offset=0.5):
    """
    multiwrap_st: add a stereo (i.e. shape (2,n)) array elementwise to a longer one, in place, wrapping around as many times as necessary
        could it be so easy? probably not?? but maybe...some of the time???

    :param cumulant: the longer array
    :param summand: the shorter array
    :param pos: the index of <cumulant> with which to align <summand>
    :param offset: the relative position within <summand> to line up with <pos>
            the default, 0.5, places <summand>'s center at pos
            0.0 should place the 0th element of <summand> at pos
            1.0 should place the last element of <summand> at pos
    """

    pos = int(np.rint(pos))  # for safety's sake

    c_len = cumulant.shape[1]
    s_len = summand.shape[1]

    s_pos = int(np.around(s_len * offset))
    s_pre = s_pos
    s_post = s_len - s_pos - 1

    prewrap = 0
    postwrap = 0

    while pos - s_pre < 0:
        prewrap = min(s_pre - pos, c_len)
        #        print('pre{}'.format(prewrap))
        # wrap around to the end
        cumulant[:, -prewrap:] += summand[:, :prewrap]
        s_pre -= prewrap
    while pos + s_post > (c_len - 1):
        # wrap around to the start
        postwrap = min(pos + s_post - (c_len - 1), c_len)
        #        print('post{}'.format(prewrap))
        cumulant[:, :postwrap] += summand[:, -postwrap:]
        s_post -= postwrap

    #    print('s_pre{}s_post{}'.format(s_pre,s_post))
    cumulant[:, (pos - s_pre):(pos + s_post + 1)] += summand[:, (s_pos - s_pre):(s_pos + s_post + 1)]
#    return (prewrap, postwrap)


def rendont(dist_data, datestamp):
    """
    rendont: like "render" but with all the actual rendering commented out. i.e. only calculate parameters.
    :param dist_data:
    :param datestamp:
    :return:
    """
    neptune = dist_data[0]

#    wave = np.zeros((2, wave_smp))

#    N = len(dist_data)-1
#    last_msg = 0
#    N_msg = 100

    for n, obj in enumerate(dist_data[1:], 1):
        # calculate frequency, duration, & amplitude of bloop
        f = nep_freq * obj['f_wrt_Nep']
        T = min(np.power(obj['T']/neptune['T'], FF_DUR_POWER) * FF_DURATION / abridgement_factor, 2 * wave_length)
        a = np.power(10, -FF_H_PRESCALE * obj['H'] / 20)
#        bloop_mono = bleep(f, T, a, sr=wave_SR)
        N_bloop = int(T * wave_SR)

        if f <= 60:  # for no other reason than to make these audible at all...
            scl_parm = FF_LOW_DRIVE * (1 - np.log2(f / 60))
#           bloop_mono = sclip(scl_parm * bloop_mono / a) * 2 * np.power(a, 1 / FF_LOW_BOOST)

        # calculate & apply window to bloop
        # for high eccentricity (near 1), xf_parm is small (near 0): fast attack, slow release
        # for low eccentricity (near 0), xf_parm is near 0.5: balanced attack/release
        xf_parm = np.power(1 - obj['e'], FF_MIDPT_DRIVE) / 2
#        exp_fade(bloop_mono, xf_parm, xf_halflives)

        # calculate & apply panning, and add to the cumulant
        rightness = np.sin(np.sin(obj['hEcl-lat'] * deg) * np.pi / 2)
#        bloop_stereo = pan_mono(bloop_mono, rightness) * FF_AMPLITUDE
#        del bloop_mono

        bloop_ctr_smp = obj['hEcl-lon'] / 360 * wave_smp
        # the peak of the bloop (midpt of exp_fade) should be centered:
        midpt_offset = int((0.5 - xf_parm) * T / 2 * wave_SR)
        ctr_smp = bloop_ctr_smp + midpt_offset

        # for future animation use:
        # T is length in seconds
        obj['lon_on'] = obj['hEcl-lon'] - (T/2)*360/wave_length
        obj['lon_off'] = obj['hEcl-lon'] + (T/2)*360/wave_length

#        multiwrap_st(wave, bloop_stereo, ctr_smp)
#        del bloop_stereo

#        if (n / N > last_msg / N_msg) or n == N:
#            print('{}/{}'.format(n, N))
#            last_msg += 1

#    wavfile.write('Distants.{}.wav'.format(datestamp), wave_SR, wave.transpose())


def render(dist_data, datestamp):
    neptune = dist_data[0]

    wave = np.zeros((2, wave_smp))

    N = len(dist_data)-1
    last_msg = 0
    N_msg = 100

    for n, obj in enumerate(dist_data[1:], 1):
        # calculate frequency, duration, & amplitude of bloop
        f = nep_freq * obj['f_wrt_Nep']
        T = min(np.power(obj['T']/neptune['T'], FF_DUR_POWER) * FF_DURATION / abridgement_factor, 2 * wave_length)
        a = np.power(10, -FF_H_PRESCALE * obj['H'] / 20)
        bloop_mono = bleep(f, T, a, sr=wave_SR)
        N_bloop = int(T * wave_SR)

        if f <= 60:  # for no other reason than to make these audible at all...
            scl_parm = FF_LOW_DRIVE * (1 - np.log2(f / 60))
            bloop_mono = sclip(scl_parm * bloop_mono / a) * 2 * np.power(a, 1 / FF_LOW_BOOST)

        # calculate & apply window to bloop
        # for high eccentricity (near 1), xf_parm is small (near 0): fast attack, slow release
        # for low eccentricity (near 0), xf_parm is near 0.5: balanced attack/release
        xf_parm = np.power(1 - obj['e'], FF_MIDPT_DRIVE) / 2
        exp_fade(bloop_mono, xf_parm, xf_halflives)

        # calculate & apply panning, and add to the cumulant
        rightness = np.sin(np.sin(obj['hEcl-lat'] * deg) * np.pi / 2)
        bloop_stereo = pan_mono(bloop_mono, rightness) * FF_AMPLITUDE
        del bloop_mono

        bloop_ctr_smp = obj['hEcl-lon'] / 360 * wave_smp
        # the peak of the bloop (midpt of exp_fade) should be centered:
        midpt_offset = int((0.5 - xf_parm) * T / 2 * wave_SR)
        ctr_smp = bloop_ctr_smp + midpt_offset

        # for future animation use:
        obj['lon_on'] = obj['hEcl-lon'] - (T/2)*360/wave_length
        obj['lon_off'] = obj['hEcl-lon'] + (T/2)*360/wave_length

        multiwrap_st(wave, bloop_stereo, ctr_smp)
        del bloop_stereo

        if (n / N > last_msg / N_msg) or n == N:
            print('{}/{}'.format(n, N))
            last_msg += 1

    wavfile.write('Distants.{}.wav'.format(datestamp), wave_SR, wave.transpose())
