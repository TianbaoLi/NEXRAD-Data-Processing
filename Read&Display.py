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


def startDownloading(files, date, site):
    date = date.replace('/', '')
    indexWidth = int(math.log10(len(files))) + 1
    cwd = os.getcwd()
    targetDir = cwd + '/' + site + '/' + date
    try:
        os.mkdir(targetDir, 0755)
    except OSError, e:
        if e.errno != os.errno.EEXIST:
            raise
        pass
    fig = plt.figure(frameon = False)
    for file in files:
        index = str(files.index(file)).zfill(indexWidth)
        suffix ='@' + site + '@' + date + '@' + index
        radar = setAWSConnection(file)

        #plt.show()
        #save_Fig(fig, plt, targetDir, suffix)

        #plt.show()
        saveFig(fig, plt, radar, targetDir, suffix)


def saveFig(fig, plt, radar, targetDir, suffix):
    fig.clf()
    plot_radar_images(radar, fig, reflectivityQCedPlot)
    fig.savefig(targetDir + '/' + 'Qced' + suffix, bbox_inches = 'tight', pad_inches = 0)


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


def gen_single_radar_image(display, fig, radar, plot):
    # display the lowest elevation scan data
    display.plot(plot[0], plot[2], title = '', colorbar_label = '', axislabels = ('', ''), colorbar_flag = False)
    display.set_limits((-300, 300), (-300, 300))
    display.set_aspect_ratio('equal')
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
    iter_str, site = args.split('@')
    print "Start to deal with", iter_str, site
    files = getDayFiles(iter_str, site)
    startDownloading(files, iter_str, site)
    print "Done dealing with", iter_str, site


def run_in_range(date_start, date_end, site, pool_size):
    iter = date_start
    delta = datetime.timedelta(days = 1)
    date_list = []
    while iter <= date_end:
        date_list.append(iter.strftime("%Y/%m/%d") + '@' + site)
        iter += delta

    cwd = os.getcwd()
    targetDir = cwd + '/' + site
    try:
        os.mkdir(targetDir, 0755)
    except OSError, e:
        if e.errno != os.errno.EEXIST:
            raise
        pass

    pool = multiprocessing.Pool(processes = pool_size)
    #print pool.get()
    pool.map(run_simple_process, date_list)

    pool.close()
    pool.join()



def main():
    date_start = datetime.date(2017, 03, 21)
    date_end = datetime.date(2017, 03, 23)
    site = "KATX"
    pool_size = 5
    run_in_range(date_start, date_end, site, pool_size)


if __name__ == "__main__":
    main()
