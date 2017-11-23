import matplotlib.pyplot as plt
import numpy.ma as ma
import numpy as np
import pyart.graph
import tempfile
import pyart.io
import boto



# read a volume scan file on S3. I happen to know this file exists.
s3conn = boto.connect_s3()
bucket = s3conn.get_bucket('noaa-nexrad-level2')
s3key = bucket.get_key('2015/05/15/KVWX/KVWX20150515_080737_V06.gz')
print s3key


# download to a local file, and read it
localfile = tempfile.NamedTemporaryFile()
s3key.get_contents_to_filename(localfile.name)
radar = pyart.io.read_nexrad_archive(localfile.name)



# display the lowest elevation scan data
display = pyart.graph.RadarDisplay(radar)
fig = plt.figure(figsize=(9, 12))

plots = [
    # variable-name in pyart, display-name that we want, sweep-number of radar (0=lowest ref, 1=lowest velocity)
    ['reflectivity', 'Reflectivity (dBZ)', 0],
    ['differential_reflectivity', 'Zdr (dB)', 0],
    ['differential_phase', 'Phi_DP (deg)', 0],
    ['cross_correlation_ratio', 'Rho_HV', 0],
    ['velocity', 'Velocity (m/s)', 1],
    ['spectrum_width', 'Spectrum Width', 1]
]

def plot_radar_images(plots):
    ncols = 2
    nrows = len(plots)/2
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

plot_radar_images(plots)

refl_grid = radar.get_field(0, 'reflectivity')
print refl_grid[0]
rhohv_grid = radar.get_field(0, 'cross_correlation_ratio')
zdr_grid = radar.get_field(0, 'differential_reflectivity')

# apply rudimentary quality control
reflow = np.less(refl_grid, 20)
zdrhigh = np.greater(np.abs(zdr_grid), 2.3)
rhohvlow = np.less(rhohv_grid, 0.95)
notweather = np.logical_or(reflow, np.logical_or(zdrhigh, rhohvlow))
print notweather[0]

qcrefl_grid = ma.masked_where(notweather, refl_grid)
print qcrefl_grid[0]

# let's create a new object containing only sweep=0 so we can add the QC'ed ref to it for plotting
qced = radar.extract_sweeps([0])
qced.add_field_like('reflectivity', 'reflectivityqc', qcrefl_grid)
display = pyart.graph.RadarDisplay(qced)
fig = plt.figure(figsize=(11, 5))

plots = [
    # variable-name in pyart, display-name that we want, sweep-number of radar (0=lowest ref, 1=lowest velocity)
    ['reflectivity', 'Reflectivity (dBZ)', 0],
    ['reflectivityqc', 'QCed Reflectivity (dBZ)', 0],
]
plot_radar_images(plots)
