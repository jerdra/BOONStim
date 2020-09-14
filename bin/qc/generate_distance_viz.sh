#!/bin/bash

BOONSTIM_DIR=$1
SUBJECT=$2
WORKDIR=$3

# Grab files
fields_file=$BOONSTIM_DIR/$SUBJECT/results/${SUBJECT}_fields.msh
qc_geo=$BOONSTIM_DIR/$SUBJECT/results/${SUBJECT}_qcdist.geo
output_dir=$BOONSTIM_DIR/$SUBJECT/results/qc
mkdir -p $output_dir

tmpdir=$(mktemp -d "$WORKDIR/boonstim_distqc_viz.XXXX")

(
	cd $tmpdir

	# Link in input files
	ln -s "$fields_file" "."
	ln -s "$qc_geo" "./qc.geo"

	# Run visualization command
	gmsh -0 -save -o "${SUBJECT}_qcdist.msh" qc.geo

	# Move output MSH file to destination
	mv "${SUBJECT}_qcdist.msh" "${output_dir}/"
)

rm -rf $tmpdir
