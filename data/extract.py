#! /usr/bin/env python3

import osmium
import sys
import math
from shapely.geometry import Point
from parse_poly import parse_poly
import geopandas as gp

objects = []

class OsmObject:
    def __init__(self):
        self.lat = None
        self.lon = None
        self.x = None
        self.y = None
        self.obj_id = None
        self.tags = {}

    def from_node(node):
        r = OsmObject()

        r.lat = node.location.lat
        r.lon = node.location.lon
        r.obj_id = node.id
        for k, v in node.tags:
            r.tags[k] = v

        return r

    def set_x_and_y(selves):
        raw = {'id': [], 'geometry': []}
        for o in selves:
            raw['id'].append(o.obj_id)
            raw['geometry'].append(Point(o.lon, o.lat))
        frame = gp.GeoDataFrame(raw, crs='EPSG:4326')
        frame = frame.to_crs('EPSG:31467')

        for o, p in zip(selves, frame.iloc):
            o.x = p['geometry'].x
            o.y = p['geometry'].y

class NamesHandler(osmium.SimpleHandler):
    def node(self, n):
        objects.append(OsmObject.from_node(n))

def main(osmfile, datafile, queryfile, query_type, *query_args):
    if query_type in ["fix", "ratio", "adaptive"]:
        if query_type == "fix" and len(query_args) != 2:
            print('arguments for "fix": <x_cnt> <y_cnt>')
            return 1
        if query_type == "ratio" and len(query_args) != 1:
            print('arguments for "ratio": <cnt>')
            return 1
        if query_type == "adaptive" and len(query_args) != 2:
            print('arguments for "adaptive": <shapefile> <point-count>')
            return 1
    else:
        print('query_type has to be either "fix" or "ratio", or "adaptive"')
        return 1

    print("=== Reading file", file=sys.stderr)
    NamesHandler().apply_file(osmfile)
    print("=== DONE Reading file", file=sys.stderr)

    chosen = []
    print("=== Choosing objects", file=sys.stderr)
    for o in objects:
        if "amenity" in o.tags:
            chosen.append(o)
    print("=== DONE Choosing objects", file=sys.stderr)
    print(len(objects), len(chosen))

    print("=== Converting coordinate system", file=sys.stderr)
    OsmObject.set_x_and_y(chosen)
    print("=== DONE Converting coordinate system", file=sys.stderr)

    with open(datafile, "w") as f:
        print("x,y,data", file=f)
        for o in chosen:
            print(f"{o.x},{o.y},{o.tags['amenity']}", file=f)

    x_min = min(o.x for o in chosen)
    x_max = max(o.x for o in chosen)

    y_min = min(o.y for o in chosen)
    y_max = max(o.y for o in chosen)

    y_len = y_max - y_min
    x_len = x_max - x_min
    ratio = x_len / y_len

    print(f"x-len={x_len}", file=sys.stderr)
    print(f"y-len={y_len}", file=sys.stderr)

    def gen_points(x_cnt, y_cnt):
        x_offset = x_len / x_cnt
        y_offset = y_len / y_cnt
        return [
                (x_min + (x_ind+0.5) * x_offset, y_min + (y_ind+0.5) * y_offset)
                for x_ind in range(x_cnt)
                for y_ind in range(y_cnt)
                ]
                

    if query_type == "fix":
        x_cnt, y_cnt = map(int, query_args)
        print(f"ratio={x_len/y_len} ({x_cnt/y_cnt})", file=sys.stderr)
        points = gen_points(x_cnt, y_cnt)
    elif query_type == "ratio":
        cnt = int(query_args[0])
        x_cnt = math.sqrt(cnt * ratio)
        y_cnt = cnt / x_cnt

        x_cnt = int(x_cnt)
        y_cnt = int(y_cnt)
        print(f"cnt={x_cnt*y_cnt}", file=sys.stderr)
        points = gen_points(x_cnt, y_cnt)
    elif query_type == "adaptive":
        polylines = open(query_args[0], "r").readlines()
        polygon = parse_poly(polylines)

        def gen(x_cnt):
            y_cnt = int(x_cnt * ratio)
            x_offset = x_len / x_cnt
            y_offset = y_len / y_cnt
            diff = math.sqrt(x_offset**2 + y_offset**2) + 1e-4

            points = gen_points(x_cnt, y_cnt)
            return [p for p in points if polygon.contains(Point(p[0], p[1])) or polygon.distance(Point(p[0], p[1])) <= diff]
        
        cnt = int(query_args[1])
        lx = int(math.sqrt(cnt * ratio))

        clx = len(gen(lx))
        # clx has to increase to cnt, so lx has to increase to lx * math.sqrt(cnt/clx) if we assume, that it scales 'normally'
        # we multiply by 4 to give some wiggle roome
        rx = int(4.0 * lx * math.sqrt(cnt / clx))
        while rx - lx > 1:
            mx = (rx + lx) // 2
            if len(gen(mx)) > cnt:
                rx = mx
            else:
                lx = mx
        points = gen(lx)
        print(f"x_cnt={lx} y_cnt={int(lx*ratio)}")
        print(f"cnt={len(points)}")

    with open(queryfile, "w") as f:
        print("x,y", file=f)
        for point in points:
            print(f"{point[0]},{point[1]}", file=f)
    return 0

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python %s <osmfile> <data-file> <query-file> <type> <*parameters>" % sys.argv[0], file=sys.stderr)
        sys.exit(-1)
    exit(main(*sys.argv[1:]))
