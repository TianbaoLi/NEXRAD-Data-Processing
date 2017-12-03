import datetime
import matplotlib
matplotlib.use('agg', warn = False, force = True)
import matplotlib.pyplot as plt
import multiprocessing
import numpy.ma as ma
import numpy as np
import pyart
import tempfile
import boto
import math
import os
from xml.dom import minidom
from urllib import urlopen
from pyart.graph import common
from PIL import Image
import matplotlib.cm as cm
import matplotlib.colors as colors
from pyart import io


plotsAll = [
    # variable-name in pyart, display-name that we want, sweep-number of radar (0=lowest ref, 1=lowest velocity)
    ['reflectivity', 'Reflectivity (dBZ)', 0],
    ['differential_reflectivity', 'Zdr (dB)', 0],
    ['differential_phase', 'Phi_DP (deg)', 0],
    ['cross_correlation_ratio', 'Rho_HV', 0],
    ['velocity', 'Velocity (m/s)', 1],
    ['spectrum_width', 'Spectrum Width', 1]
]
reflectivityPlot = ['reflectivity', 'Reflectivity (dBZ)', 0]
reflectivityQCedPlot = ['reflectivityqc', 'QCed Reflectivity (dBZ)', 0]


def getListText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)


def getDayFiles(date, site):
    bucketURL = "http://noaa-nexrad-level2.s3.amazonaws.com"
    dirListURL = bucketURL + "/?prefix=" + date + "/" + site

    print "listing files from %s" % dirListURL

    xmldoc = minidom.parse(urlopen(dirListURL))
    itemlist = xmldoc.getElementsByTagName('Key')
    print date, len(itemlist), "keys found"

    # For this test, WCT is downloaded and unzipped directly in the working directory
    # The output files are going in 'output'
    # http://www.ncdc.noaa.gov/wct/install.php
    files = []
    for x in itemlist:
        file = getListText(x.childNodes)
        files.append(file)

    return files


def cutInnerCircle(filename):
    r2 = 300
    fig_origin = Image.open(filename + '.png').convert("RGBA")
    fig_origin = fig_origin.resize((r2, r2))
    circle_limit = Image.new('RGBA', (r2, r2), (255, 255, 255, 0))
    pim_fig_origin = fig_origin.load()
    pim_circle_limit = circle_limit.load()

    r = float(r2 / 2)
    for i in range(r2):
        for j in range(r2):
            lx = abs(i - r + 0.5)
            ly = abs(j - r + 0.5)
            l = pow(lx, 2) + pow(ly, 2)
            if l <= pow(r, 2):
                pim_circle_limit[i, j] = pim_fig_origin[i, j]
    circle_limit.convert('L')
    circle_limit.save(filename + '.png')


def cutInnerDiamond(filename):
    r2 = 300
    fig_origin = Image.open(filename + '.png').convert("RGBA")
    fig_origin = fig_origin.resize((r2, r2))
    circle_limit = Image.new('RGBA', (r2, r2), (255, 255, 255, 0))
    pim_fig_origin = fig_origin.load()
    pim_circle_limit = circle_limit.load()

    r = float(r2 / 2)
    for i in range(r2):
        for j in range(r2):
            lx = abs(i - r + 0.5)
            ly = abs(j - r + 0.5)
            if lx + ly <= r:
                pim_circle_limit[i, j] = pim_fig_origin[i, j]
    circle_limit.convert('L')
    circle_limit.save(filename + '.png')


def cutInnerSquare(filename):
    r2 = int(300.0 / (2 ** 0.5))
    fig_origin = Image.open(filename + '.png').convert("RGBA")
    fig_origin = fig_origin.resize((r2, r2))
    circle_limit = Image.new('RGBA', (r2, r2), (255, 255, 255, 0))
    pim_fig_origin = fig_origin.load()
    pim_circle_limit = circle_limit.load()

    r = float(r2 / 2)
    for i in range(r2):
        for j in range(r2):
            lx = abs(i - r + 0.5)
            ly = abs(j - r + 0.5)
            if lx <= r and ly <= r:
                pim_circle_limit[i, j] = pim_fig_origin[i, j]
    circle_limit.convert('L')
    #circle_limit = circle_limit.resize((50, 50))
    circle_limit.save(filename + '.png')


def download(file, date, site, index, targetDir):
    date = date.replace('/', '')
    suffix ='@' + site + '@' + date + '@' + index
    radar = setAWSConnection(file)
    fig = plt.figure(dpi = 100, frameon = False)

    #plt.show()
    #save_Fig(fig, plt, targetDir, suffix)
    #plt.show()

    fig_name = targetDir + '/' + 'Qced' + suffix
    saveFig(fig, plt, radar, fig_name)
    plt.close(fig)

    #cutInnerCircle(fig_name)
    #cutInnerDiamond(fig_name)
    cutInnerSquare(fig_name)


def saveFig(fig, plt, radar, figname):
    fig.clf()
    plot_radar_images(radar, fig, reflectivityQCedPlot)
    fig.savefig(figname, transparent = True, bbox_inches = 'tight', pad_inches = 0.0)


def setAWSConnection(file):
    # read a volume scan file on S3. I happen to know this file exists.
    s3conn = boto.connect_s3()
    bucket = s3conn.get_bucket('noaa-nexrad-level2')
    s3key = bucket.get_key(file)
    print s3key

    # download to a local file, and read it
    localfile = tempfile.NamedTemporaryFile()
    s3key.get_contents_to_filename(localfile.name)
    radar = pyart.io.read_nexrad_archive(localfile.name)

    return radar


def gen_radar_images(display, fig, radar, plots):
    # display the lowest elevation scan data
    ncols = 2
    nrows = len(plots) / 2
    for plotno, plot in enumerate(plots, start=1):
        ax = fig.add_subplot(nrows, ncols, plotno)
        display.plot(plot[0], plot[2], ax=ax, title=plot[1],
                     colorbar_label='',
                     axislabels=('East-West distance from radar (km)' if plotno == 6 else '',
                                 'North-South distance from radar (km)' if plotno == 1 else ''))
        display.set_limits((-300, 300), (-300, 300), ax=ax)
        display.set_aspect_ratio('equal', ax=ax)
        display.plot_range_rings(range(100, 350, 100), lw=0.5, col='black', ax=ax)
        plt.show()


def _mask_outside(flag, data, v1, v2):
    """ Return the data masked outside of v1 and v2 when flag is True.  """
    if flag:
        data = np.ma.masked_invalid(data)
        data = np.ma.masked_outside(data, v1, v2)
    return data

def plot_ppi_mask_fixed(
    display, field, sweep=0, mask_tuple=None,
    vmin=None, vmax=None, norm=None, cmap=None, mask_outside=False,
    title=None, title_flag=True,
    axislabels=(None, None), axislabels_flag=True,
    colorbar_flag=True, colorbar_label=None,
    colorbar_orient='vertical', edges=True, gatefilter=None,
    filter_transitions=True, ax=None, fig=None,
    ticks=None, ticklabs=None, raster=None, **kwargs):

    # parse parameters
    ax, fig = common.parse_ax_fig(ax, fig)
    vmin, vmax = common.parse_vmin_vmax(display._radar, field, vmin, vmax)
    cmap = common.parse_cmap(cmap, field)

    # get data for the plot
    data = display._get_data(
        field, sweep, mask_tuple, filter_transitions, gatefilter)
    x, y = display._get_x_y(sweep, edges, filter_transitions)

    data.fill_value = -100.0
    data = data.filled()
    data = np.add(100.0, data)
    #print data

    # mask the data where outside the limits
    data = _mask_outside(mask_outside, data, vmin, vmax)

    # plot the data
    if norm is not None:  # if norm is set do not override with vmin/vmax
        vmin = vmax = None
    vmin = 0
    vmax = 255
    pm = ax.pcolormesh(
        x, y, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, **kwargs)

    if raster is not None:
        pm.set_rasterized(True)

    if title_flag:
        display._set_title(field, sweep, title, ax)

    if axislabels_flag:
        display._label_axes_ppi(axislabels, ax)

    # add plot and field to lists
    display.plots.append(pm)
    display.plot_vars.append(field)

    if colorbar_flag:
        display.plot_colorbar(
            mappable=pm, label=colorbar_label, orient=colorbar_orient,
            field=field, ax=ax, fig=fig, ticks=ticks, ticklabs=ticklabs)


def gen_single_radar_image(display, fig, radar, plot):
    # display the lowest elevation scan data
    plot_ppi_mask_fixed(display, plot[0], plot[2], title = '', colorbar_label = '', axislabels = ('', ''), colorbar_flag = False, cmap = 'gray')
    limits = 300.0 / (2 ** 0.5)
    display.set_limits((-300, 300), (-300, 300))
    display.set_aspect_ratio('equal')
    ax = fig.add_subplot(1, 1, 1)
    #display.plot_range_rings([300], lw=0.5, col='white', ax=ax)
    plt.xticks([])
    plt.yticks([])
    plt.gca().set_frame_on(False)
    #display.plot_range_rings(range(100, 350, 100), lw=0.5, col='black')


def plot_radar_images(radar, fig, plots):
    refl_grid = radar.get_field(0, 'reflectivity')
    #print refl_grid[0]
    rhohv_grid = radar.get_field(0, 'cross_correlation_ratio')
    zdr_grid = radar.get_field(0, 'differential_reflectivity')

    # apply rudimentary quality control
    reflow = np.less(refl_grid, 20)
    zdrhigh = np.greater(np.abs(zdr_grid), 2.3)
    rhohvlow = np.less(rhohv_grid, 0.95)
    notweather = np.logical_or(reflow, np.logical_or(zdrhigh, rhohvlow))
    #print notweather[0]

    qcrefl_grid = ma.masked_where(notweather, refl_grid)
    #print qcrefl_grid[0]

    # let's create a new object containing only sweep=0 so we can add the QC'ed ref to it for plotting
    qced = radar.extract_sweeps([0])
    qced.add_field_like('reflectivity', 'reflectivityqc', qcrefl_grid)
    display = pyart.graph.RadarDisplay(qced)

    #fig = plt.figure(figsize=(11, 5))

    gen_single_radar_image(display, fig, radar, plots)


def run_simple_process(args):
    file, date_str, site, index, targetDir = args.split('@')
    download(file, date_str, site, index, targetDir)


def run_in_range(date_start, date_end, site, pool_size):
    cwd = os.getcwd()
    targetDir = cwd + '/' + site
    try:
        os.mkdir(targetDir, 0755)
    except OSError, e:
        if e.errno != os.errno.EEXIST:
            raise
        pass

    iter = date_start
    delta = datetime.timedelta(days = 1)
    while iter <= date_end:
        date_str = iter.strftime("%Y/%m/%d")
        print "Start to deal with", date_str, site
        files = getDayFiles(date_str, site)
        date = date_str.replace('/', '')
        indexWidth = int(math.log10(len(files))) + 1
        cwd = os.getcwd()
        targetDir = cwd + '/' + site + '/' + date

        try:
            os.mkdir(targetDir, 0755)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise
            pass

        arg_list = []
        for file in files:
            index = str(files.index(file)).zfill(indexWidth)
            #arg_list.append(file + '@' + date_str + '@' + site  + '@' + index + '@' + targetDir)
            download(file, date_str, site, index, targetDir)

        #pool = multiprocessing.Pool(processes = pool_size)
        #pool.map(run_simple_process, arg_list)
        #pool.close()
        #pool.join()

        print "Done dealing with", date_str, site
        iter += delta


def main():
    date_start = datetime.date(2017, 01, 01)
    date_end = datetime.date(2017, 01, 31)
    site = "KATX"
    pool_size = 2
    run_in_range(date_start, date_end, site, pool_size)


if __name__ == "__main__":
    main()
