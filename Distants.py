if __name__ == "__main__":
    import sys
    import datetime as dt
    import prep
    import render_audio
    import render_video

    if len(sys.argv) == 2:
        datestamp = sys.argv[1]
    else:
        datestamp = dt.date.today().strftime("%Y%m%d")

    mpc_data = prep.fetch(datestamp)
    dist_data = prep.parse(mpc_data)
    try:
        open("Distants.%s.wav" % datestamp).close()
        render_audio.rendont(dist_data, datestamp)
    except FileNotFoundError:
        print("Distants.%s.wav not found; attempting to render..." % datestamp)
        render_audio.render(dist_data, datestamp)
    try:
        open("Distants.%s.mp4" % datestamp).close()
    except FileNotFoundError:
        print("Distants.%s.mp4 not found; attempting to render..." % datestamp)
        render_video.render(dist_data, datestamp)
