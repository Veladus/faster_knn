import geopandas as gp
import sys
from shapely.geometry import MultiPolygon, Point, Polygon

# Taken from https://wiki.openstreetmap.org/wiki/Osmosis/Polygon_Filter_File_Python_Parsing
def parse_poly(lines):
    """ Parse an Osmosis polygon filter file.

        Accept a sequence of lines from a polygon file, return a shapely.geometry.MultiPolygon object.

        http://wiki.openstreetmap.org/wiki/Osmosis/Polygon_Filter_File_Format
    """
    in_ring = False
    coords = []
    
    for (index, line) in enumerate(lines):
        if index == 0:
            # first line is junk.
            continue
        
        elif index == 1:
            # second line is the first polygon ring.
            coords.append([[], []])
            ring = coords[-1][0]
            in_ring = True
        
        elif in_ring and line.strip() == 'END':
            # we are at the end of a ring, perhaps with more to come.
            in_ring = False
    
        elif in_ring:
            # we are in a ring and picking up new coordinates.
            data = list(map(float, line.split()))
            point = Point(data[0], data[1])
            frame = gp.GeoDataFrame({'geometry': [point]}, crs='EPSG:4326')
            frame = frame.to_crs('EPSG:31467')

            x = frame.iloc[0]['geometry'].x
            y = frame.iloc[0]['geometry'].y

            ring.append([x, y])
    
        elif not in_ring and line.strip() == 'END':
            # we are at the end of the whole polygon.
            break
    
        elif not in_ring and line.startswith('!'):
            # we are at the start of a polygon part hole.
            coords[-1][1].append([])
            ring = coords[-1][1][-1]
            in_ring = True
    
        elif not in_ring:
            # we are at the start of a polygon part.
            coords.append([[], []])
            ring = coords[-1][0]
            in_ring = True
    
    return MultiPolygon(coords)
