#!/bin/bash

curl https://download.geofabrik.de/europe/germany/berlin-latest.osm.pbf --output raw.osm.pbf
curl https://download.geofabrik.de/europe/germany/berlin.poly --output shape.poly
python extract.py raw.osm.pbf data.csv grid.csv adaptive shape.poly 60000
