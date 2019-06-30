import datetime as dt
import numpy as np
from scipy.optimize import newton

AU = 149597870700  # 1 astronomical unit (AU) in [meters]. (Square brackets henceforth denote units.)
deg = np.pi / 180
mu = 1.3271244001798698e20  # Kepler's constant (GM) for heliocentric system in [m^3/s^2]
                            # src: HORIZONS doc
si_yr = 365.256363004       # sidereal year in SI days
                            # src: ...hmm
UT1_minus_UTC = 0.3554784   # deviation of UT1 from UTC on 2000-01-01 in [s]
                            # src: http://maia.usno.navy.mil/ser7/finals2000A.data
J2000_UT1 = dt.datetime(2000, 1, 1, 11, 58, 56, 171000)
J2000_UTC = J2000_UT1 + dt.timedelta(microseconds=int(-1e6 * UT1_minus_UTC))

dB = lambda x: 20 * np.log10(np.abs(x))

neptune = {  # Osculating elements from HORIZONS for UT 2019-Jun-28 00:00
    # TODO: figure out how to fetch & parse this from HORIZONS automatically
    # from ELEMENTS table for body 8 on HORIZONS, strictly speaking the barycenter of Neptune and its moons
    'Name': 'Neptune',
    'e': 8.745186326850268E-03,  # eccentricity [1] ('EC' on HORIZONS)
    'i': 1.770711776857225E+00,  # inclination [deg] ('IN')
    'a': 3.015250948814777E+01,  # semimajor axis [AU] ('A')
    'T': 6.047452073920306E+04,  # sidereal orbital period [days] ('PR')
    # from OBSERVER table, QUANTITIES=18,19
    'hEcl-lon': 346.5858,  # longitude on the heliocentric celestial sphere w/r/t the ecliptic
    'hEcl-lat': -1.0108,  # latitude on same
    'r': 29.93583185581,  # heliocentric radial distance [AU]
    'rdot': -0.0272290,  # heliocentric radial velocity [km/s]
}
neptune['Ty'] = neptune['T'] / si_yr  # sidereal orbital period in sidereal years...sidereally
neptune['q'] = neptune['a'] * (1 - neptune['e'])  # perihelion distance [AU]
neptune['Q'] = neptune['a'] * (1 + neptune['e'])  # aphelion distance [AU]
nep_freq = 720  # [Hz]

neptune_mean = {
    # "Mean Orbital Elements (J2000)" per https://nssdc.gsfc.nasa.gov/planetary/factsheet/neptunefact.html
    'Name': 'Neptune',
    'a_mean': 30.06896348,
    'e_mean': 0.00858587,
    'i_mean': 1.76917,
    'Node_mean': 131.72169,
    'q_mean': 30.06896348 * (1 - 0.00858587),  # perihelion distance [AU]
    'Q_mean': 30.06896348 * (1 + 0.00858587),  # aphelion distance [AU]
    #    'Longitude_of_perhelion_mean': 44.97135,
    #    'Mean_longitude': 304.88003
    'f_wrt_Nep': 1.0,
    'H': -7.0,  # https://arxiv.org/ftp/arxiv/papers/1604/1604.00518.pdf passband V ¯\_(ツ)_/¯
}
neptune.update(neptune_mean)


def sfloat(string):
    string = string.strip()
    if string is not '':
        fl = float(string)
    else:
        fl = np.nan
    return fl


def sint(string):
    string = string.strip()
    if string is not '':
        if string[0] == '(' and string[-1] == ')':
            string = string[1:-1]
        integer = int(string)
    else:
        integer = np.nan
    return integer


def unpack_epoch(packed_epoch):
    year = 100 * (ord(packed_epoch[0]) - 55) + int(packed_epoch[1:3])

    m_pack = packed_epoch[3]
    m_ord = ord(m_pack)
    if m_ord >= 49 and m_ord <= 57:
        month = m_ord - 48
    elif m_ord >= 65 and m_ord <= 67:
        month = m_ord - 55
    else:
        raise ValueError('Invalid epoch month: {}'.format(m_pack))

    d_pack = packed_epoch[4]
    d_ord = ord(d_pack)
    if d_ord >= 49 and d_ord <= 57:
        day = d_ord - 48
    elif d_ord >= 65 and d_ord <= 86:
        day = d_ord - 55
    else:
        raise ValueError('Invalid epoch day: {}'.format(d_pack))

    dtobj = dt.datetime(year=year, month=month, day=day) - dt.timedelta(minutes=1, seconds=4, milliseconds=180)
    # the above timedelta is the offset between TT and UTC as of 2000-01-01
    # which is already probably pointless, and the correct offset for any epoch I'm ever likely to encounter will be at most a few more seconds off from that
    # so
    # fine
    return dtobj


def rel_err(obs, exp):
    if exp != 0:
        result = abs((obs-exp)/exp)
    else:
        result = np.nan
    return result


def idx(token, db):
    """
    idx: Find the index of an object in dictionary `db` with the name or number `token`. It had better be unique.
    """
    by_name = [record for record in db if 'Name' in record and record['Name'] == token]
    by_number = [record for record in db if 'Number' in record and record['Number'] == token]
    if len(by_name) == 1:
        return db.index(by_name[0])
    elif len(by_number) == 1:
        return db.index(by_number[0])
    else:
        raise ValueError


def dt_from_julian(J_date):
    """"
    dt_from_julian: Construct the datetime object corresponding to a given Julian date using the algorithm from Meeus, "Astronomical Formulae for Calculators", as described in http://mathforum.org/library/drmath/view/51907.html

    :param float J_date: Julian date
    :return: datetime object
    """
    Jp = J_date + 0.5
    Z = int(Jp)
    F = Jp - Z

    if Z < 2299161:
        A = Z
    else:
        alpha = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - int(alpha / 4)

    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((B - D) / 30.6001)

    dayfrac = B - D - int(30.6001 * E) + F
    if E < 13.5:
        month = E - 1
    else:
        month = E - 13

    if month > 2.5:
        year = C - 4716
    else:
        year = C - 4715

    day = int(dayfrac)
    frac = dayfrac - day

    frac *= 24
    hour = int(frac)
    frac -= hour

    frac *= 60
    minute = int(frac)
    frac -= minute

    frac *= 60
    second = int(frac)
    frac -= second

    frac *= 1e6
    microsecond = int(frac)

    return dt.datetime(year=year, month=month, day=day,
                       hour=hour, minute=minute, second=second, microsecond=microsecond)


def hEcl_from_kep(a, ecc, omega, Omega, inc, M_0, t_0=J2000_UTC, t=dt.datetime.utcnow()):
    """
    hEcl_from_kep: calculate heliocentric ecliptic spherical coordinates from Keplerian orbital elements
        - shoutout to Dr. J.B. Tatum for "Computation of an Ephemeris",
        http://astrowww.phys.uvic.ca/~tatum/celmechs/celm10.pdf, pages 4-7
        - epochs for omega, Omega, & inc all assumed, hoped, prayed to be J2000.0,
        because that's what the Minor Planet Center files are...right???

    :param a: semi-major axis [AU]
    :param ecc: eccentricity [1]
    :param omega: argument of periapsis, J2000.0 [deg]
    :param Omega: longitude of ascending node, J2000.0 [deg]
    :param inc: inclination, J2000.0 [deg]
    :param M_0: mean anomaly, epoch t_0 [deg]
    :param t_0: reference epoch [datetime.datetime object, UTC please & thank you]
    :param t: time for which to calculate the vector [datetime.datetime object, UTC please & thank you]
    :returns: (r_dist, hEcl_lon, hEcl_lat)
        - r_dist: radial distance from (center of) Sun [AU]
        - hEcl_lon: heliocentric longitude w/r/t J2000.0 equinox [deg]
        - hEcl_lat: heliocentric latitude w/r/t J2000.0 equinox [deg]
    """

    # astronomers like degrees, but numpy likes radians
    # when in Rome...
    a *= AU
    omega *= deg
    Omega *= deg
    inc *= deg
    M_0 *= deg

    delta_t = (t - t_0).total_seconds()
    M_t = M_0 + delta_t * np.sqrt(mu / a ** 3)  # mean anomaly at epoch t

    # Step 1: solve Kepler's equation for the eccentric anomaly, using, say, Newton's method
    kep_func = lambda E: E - ecc * np.sin(E) - M_t
    kep_prime = lambda E: 1 - ecc * np.cos(E)
    # Note: Providing the second derivative (for parabolic Halley's method) appears to be slower overall.
    # Apparently the faster convergence isn't worth the extra function calls.
    # The next 2 lines have been retained for historical interest.
    #    kep_prime2 = lambda E: ecc*np.sin(E)
    #    E_t = newton(kep_func, M_t, fprime=kep_prime, fprime2=kep_prime2)
    E_t = newton(kep_func, M_t, fprime=kep_prime)

    # Step 2: calculate true anomaly & radial distance from eccentric anomaly
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E_t / 2), np.sqrt(1 - ecc) * np.cos(E_t / 2))
    r_dist = a * (1 - ecc * np.cos(E_t))

    # Step 3: convert from orbital coordinate frame to ecliptic frame
    # refs: pages 15-16 of http://astrowww.phys.uvic.ca/~tatum/celmechs/celm3.pdf
    # and pages 6-7 of http://astrowww.phys.uvic.ca/~tatum/celmechs/celm10.pdf:
    # if X is the object, P is its periapsis, and N is the orthogonal projection of X onto the ecliptic,
    # then you've got yourself a right spherical triangle with corners at the ascending node, X, and N,
    # with hypotenuse |nodePX| = omega+nu
    # longitudinal side |nodeN| = hEcl_lon-Omega
    # and latitudinal side |XN| = hEcl_lat
    # all three sides are arcs of great circles.
    # the cotangent formula for spherical triangles relates two angles, the side (arc) between them, and a second side.
    # with the inclination, the right angle, and the hypotenuse known, you can solve for hEcl_lon
    # then with the inclination, the right angle, and the longitudinal side you can solve for hEcl_lat
    # so: do that
    hEcl_lon = np.arctan2(np.cos(inc) * np.sin(omega + nu), np.cos(omega + nu)) + Omega
    hEcl_lon %= 2 * np.pi  # would you jump off a bridge if HORIZONS jumped? probably
    #    hEcl_lat = np.arctan2(np.sin(hEcl_lon-Omega)*np.sin(inc), np.cos(inc))
    hEcl_lat = np.arctan(np.sin(hEcl_lon - Omega) * np.tan(inc))

    # ...emoR ni nehw
    r_dist /= AU
    hEcl_lon /= deg
    hEcl_lat /= deg

    return r_dist, hEcl_lon, hEcl_lat


def cart_from_kep(a, ecc, omega, Omega, inc, M_0, t_0=J2000_UTC, t=dt.datetime.utcnow()):
    """
    cart_from_kep: calculate heliocentric Cartesian orbital state vector from Keplerian orbital elements
        this began by basically copying https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
        so, you know, shout-out to René Schwarz
        epoch for omega, Omega, & inc assumed, hoped, prayed to be J2000.0 because that's what the Minor Planet Center files have...right???

    :param a: semi-major axis [AU]
    :param ecc: eccentricity [1]
    :param omega: argument of periapsis, J2000.0 [deg]
    :param Omega: longitude of ascending node, J2000.0 [deg]
    :param inc: inclination, J2000.0 [deg]
    :param M_0: mean anomaly, epoch t_0 [deg]
    :param t_0: reference epoch [datetime.datetime object, UTC please & thank you]
    :param t: time for which to calculate the vector [datetime.datetime object, UTC please & thank you]
    :returns: tuple (r, rdot)
        - r: heliocentric Cartesian position vector [AU]
        - rdot: heliocentric Cartesian velocity vector [km/s]
    """

    # astronomers like degrees, but numpy likes radians
    # when in Rome...
    a *= AU
    omega *= deg
    Omega *= deg
    inc *= deg
    M_0 *= deg

    delta_t = (t - t_0).total_seconds()
    M_t = M_0 + delta_t * np.sqrt(
        mu / a ** 3)  # mean anomaly at epoch t, assuming simple Keplerian orbit (unperturbed etc)

    # step 1: solve Kepler's equation for eccentric anomaly
    kep_func = lambda E: E - ecc * np.sin(E) - M_t
    kep_prime = lambda E: 1 - ecc * np.cos(E)
    # Note: Providing the second derivative (for parabolic Halley's method) appears to be slower overall.
    # Apparently the faster convergence isn't worth the extra function calls.
    # The lines below have been retained for historical interest.
    #    kep_prime2 = lambda E: ecc*np.sin(E)
    #    E_t = newton(kep_func, M_t, fprime=kep_prime, fprime2=kep_prime2)
    E_t = newton(kep_func, M_t, fprime=kep_prime)
    cos_E, sin_E = np.cos(E_t), np.sin(E_t)

    # step 2: calculate true anomaly & radial distance from eccentric anomaly
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E_t / 2), np.sqrt(1 - ecc) * np.cos(E_t / 2))
    r_dist = a * (1 - ecc * cos_E)

    # step 3: construct position and velocity vectors in orbital frame
    # (z-axis normal to orbital plane, x-axis oriented towards periapsis of orbit)
    o = r_dist * np.array([np.cos(nu), np.sin(nu), 0])
    odot = np.sqrt(mu * a) / r_dist * np.array([-sin_E, np.sqrt(1 - ecc ** 2) * cos_E, 0])

    # step 4: transform orbital-frame vectors into heliocentric Cartesian coordinates
    # ho boy

    cos_om, sin_om = np.cos(omega), np.sin(omega)
    cos_Om, sin_Om = np.cos(Omega), np.sin(Omega)
    cos_inc, sin_inc = np.cos(inc), np.sin(inc)

    R_mat = np.array([[
        cos_om * cos_Om - sin_om * cos_inc * sin_Om,
        -(sin_om * cos_Om + cos_om * cos_inc * sin_Om),
        0
    ],
        [
            cos_om * sin_Om + sin_om * cos_inc * cos_Om,
            -sin_om * sin_Om + cos_om * cos_inc * cos_Om,
            0
        ],
        [
            sin_om * sin_inc,
            cos_om * sin_inc,
            0
        ]])

    r = R_mat.dot(o)
    rdot = R_mat.dot(odot)

    '''
    # matmul hater version:
    r = [
        o[0]*(cos_om*cos_Om - sin_om*cos_inc*sin_Om) - o[1]*(sin_om*cos_Om + cos_om*cos_inc*sin_Om),
        o[0]*(cos_om*sin_Om + sin_om*cos_inc*cos_Om) + o[1]*(-sin_om*sin_Om + cos_om*cos_inc*cos_Om),
        o[0]*sin_om*sin_inc + o[1]*cos_om*sin_inc
    ]
    rdot = [
        odot[0]*(cos_om*cos_Om - sin_om*cos_inc*sin_Om) - odot[1]*(sin_om*cos_Om + cos_om*cos_inc*sin_Om),
        odot[0]*(cos_om*sin_Om + sin_om*cos_inc*cos_Om) + odot[1]*(-sin_om*sin_Om + cos_om*cos_inc*cos_Om),
        odot[0]*sin_om*sin_inc + odot[1]*cos_om*sin_inc
    ]
    '''

    # ...emoR ni nehw
    r /= AU
    rdot /= 1e3  # from m/s to km/s

    return r, rdot


def nep_res(j, k):
    """
    nep_res: Calculates approximate semimajor axis for specified mean-motion resonances with Neptune.

    :param j: base of (j+k):j period resonances
    :param k: offset, e.g. for 2:1 resonances, j=1, k=1
    :returns: semimajor axis for the specified resonance [AU]
    """

    if 'a_mean' in neptune:
        a_nep = neptune['a_mean']
    elif 'a' in neptune:
        a_nep = neptune['a']
    else:
        raise KeyError("what's Neptune's semimajor axis again??")
    return a_nep * np.power((j+k)/j, 2/3)


def nep_freq_ratio(a):
    """
    nep_freq_ratio: Approximate orbital frequency relative to Neptune
    bypassing T_from_a() because pointless

    :param a: semimajor axis of distant object [AU]
    :return: orbital frequency relative to Neptune
    """

    if 'a_mean' in neptune:
        a_nep = neptune['a_mean']
    elif 'a' in neptune:
        a_nep = neptune['a']
    else:
        raise KeyError("what's Neptune's semimajor axis again??")
    return np.power(a_nep/a, 3/2)


def T_from_a(a):
    """
    T_from_a: Keplerian approximation of orbital period around the Sun of much smaller object, given its semimajor axis

    :param a: semimajor axis [AU]
    :return: T: sidereal orbital period [days]
    """
    a *= AU
    T = 2*np.pi*np.sqrt(np.power(a, 3.0)/mu)/86400
    return T


def time_from_lon(lon, wave_length):
    """
    time_from_lon: Given the piece's total rendered length, map a heliocentric longitude to a time within the piece
    as a min:sec.ms string

    :param lon: heliocentric longitude [deg]
    :param wave_length: total length of the piece [s]
    :return: formatted string
    """
    time_sec = wave_length*lon/360
    time_min = int(time_sec//60)
    time_sec %= 60
    return '{}:{:06.3f}'.format(time_min, time_sec)


# TODO: possibly factor all of the above core functions out from the below scriptier functions


def fetch(datestamp=None):
    """
    fetch: read the MPC's distant.txt file, downloading if not locally present

    :param datestamp: date, as a string in the form YYYYMMDD. It's assumed that either a distant.txt file from that date
        is locally cached, or we should download today's.
    :return: currents_data: lines of the distant.txt file, to be parsed by parse()
    """
    import urllib.request
    import urllib.error

    if datestamp is not None:
        currents_filename = "distant.{}.txt".format(datestamp)
        currents_fd = open(currents_filename, 'r')
        currents_data: list = currents_fd.readlines()
        currents_fd.close()
    else:
        datestamp = dt.date.today().strftime("%Y%m%d")
        currents_filename = "distant.{}.txt".format(datestamp)
        try:
            currents_fd = open(currents_filename, 'r')
            currents_data: list = currents_fd.readlines()
            currents_fd.close()
        except FileNotFoundError:
            print("fetching today's distant.txt from the Minor Planet Center...")
            mpc_http_response = urllib.request.urlopen("https://www.minorplanetcenter.net/iau/ECS/MPCAT/current/distant.txt")
            if mpc_http_response.code != 200:
                raise urllib.error.HTTPError
            currents_data = mpc_http_response.read().decode('utf-8').splitlines()
            with open(currents_filename, 'w') as currents_fd:
                currents_fd.writelines(currents_data)

    return currents_data


def parse(mpc_data):
    """
    parse: Read the fixed-width columns out of the lines of MPC data

    :param mpc_data: list of lines from MPC data file
    :return: dist_data: parsed data structure as list of dicts with named fields
    """
    dist_data = [neptune]
    for txt in mpc_data:
        try:
            txt_tmp = {
                # 'PackedDesig': txt[0:7].strip(),
                'H': sfloat(txt[8:13]),  # absolute magnitude
                'G': sfloat(txt[14:19]),  # slope parameter
                'Epoch': unpack_epoch(txt[20:25]),  # epoch (read in packed form, .0 TT)
                # NB: .0 TT
                # to date, I have not cared enough yet to calculate/correct for the difference between TT and UTC here
                'M': sfloat(txt[26:35]),  # mean anomaly at the epoch [deg]
                'Peri': sfloat(txt[37:46]),  # argument of perihelion, J2000.0 [deg]
                'Node': sfloat(txt[48:57]),  # longitude of the ascending node, J2000.0 [deg]
                'i': sfloat(txt[59:68]),  # inclination to the ecliptic, J2000.0 [deg]
                'e': sfloat(txt[70:79]),  # orbital eccentricity
                'n': sfloat(txt[80:91]),  # mean daily motion [deg/day]
                'a': sfloat(txt[92:103]),  # semimajor axis [AU]
                'Ref': txt[107:116],
                'Num_obs': sint(txt[117:122]),
                'Num_opps': sint(txt[123:126]),
                'Number': sint(txt[166:175]),
                'Name': txt[175:194].strip(),
                'Last_obs': txt[194:202],
            }
        except:
            print(txt)
            raise
        if not np.isnan(txt_tmp['H']):
            dist_data.append(txt_tmp)

    for d in dist_data[1:]:
        r, hEcl_lon, hEcl_lat = hEcl_from_kep(d['a'], d['e'], d['Peri'], d['Node'], d['i'], d['M'], d['Epoch'])
        d['r'] = r
        d['hEcl-lon'] = hEcl_lon
        d['hEcl-lat'] = hEcl_lat
        d['T'] = T_from_a(d['a'])
        d['Ty'] = d['T'] / si_yr
        d['f_wrt_Nep'] = nep_freq_ratio(d['a'])
        d['q'] = d['a'] * (1 - d['e'])
        d['Q'] = d['a'] * (1 + d['e'])
        # Tisserand's parameter w/r/t Neptune. Shoutout to Wikipedia! This could be useful for characterizing objects
        # and changing their sounds. Or not.
        d['Tiss_Nep'] = neptune['a'] / d['a'] + 2 * np.sqrt(d['a'] / neptune['a'] * (1 - d['e'] ** 2)) * np.cos(d['i'])

    return dist_data
