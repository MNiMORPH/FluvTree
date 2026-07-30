# drainageNetworkGRASS.py
# Based off of:
# ToroChannelNetwork.py from work for AGU 2017
# A. Wickert, 06 DEC 2025

##################
# IMPORT MODULES #
##################
# PYTHON
import numpy as np
# GRASS
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.gis import region
from grass.pygrass import vector # Change to "v"?
from grass.script import vector_db_select
from grass.pygrass.vector import Vector, VectorTopo
from grass.pygrass.raster import RasterRow
from grass.pygrass.raster.abstract import RasterAbstractBase
from grass.pygrass import utils
from grass import script as gscript
from grass.pygrass.vector.geometry import Point
import grass

DEM_original_import='DEM_original_import'
DEM='dem'
accumulation='accumulation'
#accumulation_onmap='accumulation_onmap'
streams_all='streams_all'
draindir='draindir'
cellArea_meters2='cellArea_meters2'

# Wickert g.extension. Still necessary?
r.cell_area(output=cellArea_meters2, units='m2', overwrite=True)

# Fails with non-square cells!
r.watershed(elevation=DEM, flow=cellArea_meters2, accumulation=accumulation, flags='s', overwrite=True)
#r.mapcalc(accumulation_onmap+' = '+accumulation+' * ('+accumulation+' > 0)', overwrite=True)
#r.mapcalc('tmp'+' = if(isnull('+accumulation_onmap+'),null(),'+DEM_original_import+')', overwrite=True)
r.mapcalc(accumulation+' = if(isnull('+DEM+'),null(),'+accumulation+')', overwrite=True)
#g.rename(raster=('tmp',DEM), overwrite=True)
r.stream_extract(elevation=DEM, accumulation=accumulation, stream_raster=streams_all, stream_vector=streams_all, threshold=3000, direction=draindir, d8cut=0, overwrite=True)

v.stream_network(map=streams_all)

# Run in shell or rewrite
# 985: NW
# 1258: Center
v.stream.profiler cat=1258 streams=streams_all outstream=stream direction=upstream elevation=dem accumulation=accumulation units=m2 --o

# Topo and database table are in the same order
colNames = np.array(list(gscript.vector_db_select('stream', layer=1)['columns']))
colValues = np.array(list(gscript.vector_db_select('stream', layer=1)['values'].values()))
number_of_segments = colValues.shape[0]
cats = colValues[:,colNames == 'cat'].astype(int).squeeze() # stream cats: just use cat instead of new col stream
tocats = colValues[:,colNames == 'tostream'].astype(int).squeeze() # tostream cats

streamTopo = VectorTopo('stream')
streamTopo.open('r')
lines = []
for item in streamTopo:
    if type(item) is grass.pygrass.vector.geometry.Line:
        lines.append(item)
streamTopo.rewind()

# Get elevations
zRast = RasterRow(DEM)
zRast.open('r')
ARast = RasterRow(accumulation)
ARast.open('r')

x = []
y = []
z = []
A = []
for line in lines:
    for point in line:
        val = np.abs(ARast.get_value((point.x, point.y)))
        if val >= 10000000:
            A.append(val)
        else:
            A.append(np.nan)
        x.append(point.x)
        y.append(point.y)
        z.append(zRast.get_value((point.x, point.y)))
x = np.array(x)
y = np.array(y)
z = np.array(z)
A = np.array(A)
dist = np.hstack((0, np.cumsum((np.diff(x)**2 + np.diff(y)**2)**.5)))

def running_mean(x, y, dx):
    _x = x[0]
    xout = []
    yout = []
    while _x < x[-1]:
        _y = np.nanmean(y[(x >= _x) * (x < _x+dx)])
        _x = np.nanmean(x[(x >= _x) * (x < _x+dx)])
        xout.append(_x)
        yout.append(_y)
        _x += dx
    return xout, yout
        
S = -np.diff(z) / np.diff(dist)
AS = (A[1:] + A[:-1])/2.
distS = (dist[1:] + dist[:-1])/2.

xm, Sm = running_mean(distS, S, 300)
xm, ASm = running_mean(distS, AS, 300)

plt.plot(dist, z)

# chi = np.cumsum(A**.6 * np.diff(dist)) # 2025 note: ASm?

"""
for line in lines:
    for point in line:
        point.z = zRast.get_value((point.x, point.y))
        
for _li in range(len(lines)):
    for _pi in range(len(lines[_li])):
        lines[_li][_pi].z = zRast.get_value((lines[_li][_pi].x, lines[_li][_pi].y))

# if downstream only, then one after the other:
local_downstream_distance_lines = []
for line in lines:
    tmpList = []
    for point in line:
        tmpList.append(
"""


