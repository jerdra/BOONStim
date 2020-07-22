#!/bin/bash

BOONSTIM_DIR=$1
SUBJECT=$2
WORKDIR=$3
VIEW_FILE=$4
OUTPUT_FILE=$5

# Grab files
fields_file=$BOONSTIM_DIR/$SUBJECT/results/${SUBJECT}_fields.msh
coil_geo=$BOONSTIM_DIR/$SUBJECT/results/${SUBJECT}_optimized_coil.geo

tmpdir=$(mktemp -d "$WORKDIR/boonstim_viz.XXXX")

(
	cd $tmpdir

	# Link in input files
	ln -s "$fields_file" "./sub.msh"
	ln -s "$coil_geo" "./coil.geo"
	ln -s "$VIEW_FILE" "./view.geo"

	# Run visualization command
	gmsh -0 view.geo

	# Move output PNG file to destination
	mv "sub.png" $OUTPUT_FILE
)
