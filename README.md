# Study for distant objects (Distants)
To produce a piece of generative music based on the properties of distant solar system bodies, as defined and conveniently tabulated by the Minor Planet Center. This piece changes, slowly, every time it's rendered--moreso when the MPC's data file is updated, especially as new objects are discovered or significant corrections are made to the orbits of known objects.

The gist for these versions is for each distant object to generate a sound-event. The timing of the sound-event is based on the object's heliocentric longitude, like a radar beam sweeping around the celestial sphere, or a music box. The frequency of each sound is based on the object's orbital frequency relative to Neptune, where Neptune itself is somewhat arbitrarily assigned 1440 Hz (but does not generate a sound). A large number of distant objects (including Pluto) are in 3:2 orbital resonance with Neptune: they complete 2 orbits for every 3 of Neptune's. This leads to a large number of sound-events with frequencies near 960 Hz (2/3rds of Neptune's frequency). Other properties of the object determine other parameters of the sound-event, with varying degrees of arbitrariness & futzing involved:
* Absolute magnitude: maximum amplitude
* Orbital period: duration of sound-event (higher pitches correlate to shorter durations)
* Orbital eccentricity: asymmetry of amplitude envelope (circular orbits fade in & out smoothly; highly eccentric orbits have fast attack and slow release)
* Heliocentric latitude: panning (stereo position)

## Future Work
My intent is for the full piece to be 60 minutes long, but rendering the piece takes an extravagantly long time even at the 10x "abridged" version this code currently generates. Most of this time is spent calculating the (exponential) amplitude envelopes for every sound-event. It might be better to quantize the durations and asymmetries into a limited number of categories and generate a lookup table (or cache them, so as not to calculate unused combinations at all).

The `render` function in `render_video.py` generates an animation showing which objects are sounding as the piece progresses, but this display could be a lot clearer: https://youtu.be/avZkv_wdZDY 
Currently all points and labels look the same. I would like to change either the points or the labels or both to visualize the pitch, amplitude, and/or eccentricity of each object, to make it easier to tell which sound corresponds to which object. Also, the audio has to be manually added to the video because I haven't figured out how to make ffmpeg do it from matplotlib.

Osculating elements for Neptune, from NASA's HORIZONS web service, are hardcoded, because the alternative is Telnet, and I haven't figured out how to automate that.

## Credit Where Credit's Due
* This research has made use of data and/or services provided by the International Astronomical Union's Minor Planet Center.
    * namely: https://www.minorplanetcenter.net/iau/MPCORB/Distant.txt (initially)
    * https://www.minorplanetcenter.net/iau/ECS/MPCAT/current/distant.txt (currently)
    * http://minorplanetcenter.net/Extended_Files/distant_extended.json.gz
    * and the MPC's documentation regarding each of them
 * and by the NASA Jet Propulsion Laboratory HORIZONS ephemeride system: https://ssd.jpl.nasa.gov/horizons.cgi

## References
* Giorgini, J.D., Yeomans, D.K., Chamberlin, A.B., Chodas, P.W., Jacobson, R.A., Keesey, M.S., Lieske, J.H., Ostro, S.J., Standish, E.M., Wimberly, R.N., "JPL's On-Line Solar System Data Service", Bulletin of the American Astronomical Society, Vol 28, No. 3, p. 1158, 1996.
    * HORIZONS documentation. Available at [ftp://ssd.jpl.nasa.gov/pub/ssd/Horizons_doc.pdf](ftp://ssd.jpl.nasa.gov/pub/ssd/Horizons_doc.pdf). Accessed 23 January 2018.
* The Math Forum. 2001. "Julian to Calendar Date Conversion". Available at http://mathforum.org/library/drmath/view/51907.html. Accessed 23 January 2018.
    * Meeus, Jean. 1979. "Astronomical Formulae for Calculators". Cited in ibid.
* Tatum, J.B. 2017. Celestial Mechanics. Available at http://astrowww.phys.uvic.ca/~tatum/celmechs.html. Accessed 06 February 2018.
    * in particular http://astrowww.phys.uvic.ca/~tatum/celmechs/celm3.pdf pp. 15-16 and http://astrowww.phys.uvic.ca/~tatum/celmechs/celm10.pdf pp. 6-7
