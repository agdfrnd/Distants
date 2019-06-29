import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import logging

logging.basicConfig(level=logging.INFO)

# plt.rcParams['animation.ffmpeg_path'] = "D:\\Users\\alex\\Downloads\\lib\\ffmpeg.exe"

# TODO: obviously this needs to be set one place and passed from there
full_length = 3600  # [s]
abridgement_factor = 10
wave_length = full_length // abridgement_factor

fps = 5

H_cutoff = 6
anim_width = 1920
anim_height = 1080


def generate_annos(ax, dist_data):
    vid_data = []
    for dist in dist_data[1:]:
        if dist['H'] < H_cutoff:
            vid_data.append({
                'anno': ax.annotate(dist['Name'],
                                    xy=(dist['hEcl-lon'], dist['hEcl-lat']),
                                    xycoords='data',
                                    alpha=0.8
                                    ),
                'lon_on': dist['lon_on'],
                'lon_off': dist['lon_off'],
                })
    print("generated annotations for %d out of %d objects..." % (len(vid_data), len(dist_data)-1))
    return vid_data


def render(dist_data, datestamp):
    fig, ax = plt.subplots(figsize=(anim_width / 120, anim_height / 120), dpi=120)

    x = [d['hEcl-lon'] for d in dist_data[1:]]
    y = [d['hEcl-lat'] for d in dist_data[1:]]
    x_nep = dist_data[0]['hEcl-lon']
    y_nep = dist_data[0]['hEcl-lat']

    # Lagrange points of Neptune-Sun system...not that they're that informative
    L4 = (x_nep + 60) % 360
    L5 = (x_nep - 60) % 360
    L3 = (x_nep + 180) % 360

    scat = ax.scatter(x, y, s=9, alpha=0.1)
    ax.plot(x_nep, y_nep, 'ro', [L4, L5], [y_nep, y_nep], 'r.', L3, y_nep, 'r+')
    ax.set_xlim(0, 360)
    ax.set_ylim(-90, 90)
    ax.set_xlabel('heliocentric longitude')
    ax.set_ylabel('heliocentric latitude')
    ax.set_title(None)

    fig.subplots_adjust(left=0.04, right=0.96, bottom=0.05, top=1, hspace=0, wspace=0)

    ax.annotate('Nep', xy=(x_nep, y_nep), xycoords='data')
    ax.annotate('L3', xy=(L3, y_nep), xycoords='data')
    ax.annotate('L4', xy=(L4, y_nep), xycoords='data')
    ax.annotate('L5', xy=(L5, y_nep), xycoords='data')

    vid_data = generate_annos(ax, dist_data)

    lon_list = np.linspace(0, 360, wave_length * fps + 1)

    cels = []
    vis_anno_count = 0

    for lon in lon_list:
        line, = ax.plot([lon, lon], [-90, 90], 'r', alpha=0.4)
        annos = [dist['anno'] for dist in vid_data
                 if dist['lon_on'] <= lon < dist['lon_off']
                 ]
        vis_anno_count += len(annos)
        cels.append([line] + annos)
    print("avg %0.3f visible annotations per frame" % (vis_anno_count / len(lon_list)))

    ani = animation.ArtistAnimation(fig, cels, interval=1000//fps, blit=True)
    print("generated animation; attempting to write to video...")

    # ffmpeg_extra_args = ['-loglevel', 'verbose']
    # if include_audio:
    #     ffmpeg_extra_args = ffmpeg_extra_args + ['-i', "Distants.%s.wav" % datestamp,
    #                                              '-c:a', 'libvo_aacenc', '-b:a', '320k']

    aniwriter = animation.FFMpegFileWriter(fps=fps, extra_args=[])

    ani.save("Distants.%s.mp4" % datestamp, writer=aniwriter)
    return ani
