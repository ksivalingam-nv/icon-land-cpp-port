#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2026, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

set -e
##############################################################################################
### This program collects state variables of the jsbach carbon process from an ICON restart
### file and generates a corresponding initial conditions (ic_) file for a follow-up
### experiment. This helps to overcome long spinup times of the carbon cycle.
###
### The program currently only works with single file restart files. In case you have
### multi-file restart files, concatinate the patches using "cdo collgrid".
###
### Note: Tested with nco version 5.0.1 and cdo version 2.4.2
###
# Todo: Include further variables needed in simulations with natural LCC and/or disturbances
##############################################################################################
which ncatted || {
. ${MODULESHOME}/init/bash
  module load nco
}

### User Interface
EXP=picontrol       # Experiment name, just used in the file names
YEAR=2920           # Restart file year, just used in the file names
INFILE=/work/mh1496/m220053/equilibrateC/${EXP}/${EXP}_restart_atm_${YEAR}0101T000000Z_eqC.nc
OUTFILE=/work/mh1496/m220053/equilibrateC/${EXP}/ic_land_carbon_${EXP}_${YEAR}.nc
GRIDFILE=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc

PROCESS="carbon"
VARNAME_LIST="c_green    c_woods     c_reserve
              c_acid_ag1 c_water_ag1 c_ethanol_ag1 c_nonsoluble_ag1
              c_acid_bg1 c_water_bg1 c_ethanol_bg1 c_nonsoluble_bg1 c_humus_1
              c_acid_ag2 c_water_ag2 c_ethanol_ag2 c_nonsoluble_ag2
              c_acid_bg2 c_water_bg2 c_ethanol_bg2 c_nonsoluble_bg2 c_humus_2"
TILE_LIST="box land veg pft01 pft02 pft03 pft04 pft05 pft06 pft07 pft08 pft09 pft10 pft11 pft12"

### Delete existing output files ?
if [ -f  ${OUTFILE} ];then
   echo "The output file ${OUTFILE} already exists. Shall I delete it?"
   echo "RETURN =File will be deleted."
   echo "XXX    =Every other input stops this script here."
   read y
   if [ ! ${y} ] ; then
      rm ${OUTFILE}
      echo "The existing file ${OUTFILE} was deleted."
      echo ""
   else
      exit
   fi
fi

### Generate variable list
varlist=""
for VARNAME in ${VARNAME_LIST} ; do
   for TILE in ${TILE_LIST} ; do
      varlist="${varlist},${PROCESS}_${VARNAME}_${TILE}"
   done
done

### Extract the variables from restart file
cdo selvar${varlist}    ${INFILE}   ${OUTFILE}.tmp1

ncks -C -x -v time  ${OUTFILE}.tmp1 ${OUTFILE}.tmp2  # Delete variable time
ncwa -O -a time     ${OUTFILE}.tmp2 ${OUTFILE}.tmp3  # Delete dimension time
cdo setgrid,${GRIDFILE} ${OUTFILE}.tmp3 ${OUTFILE}

### Clean up
rm -f ${OUTFILE}.tmp?

### Finish
echo ""
echo "$0: Generation of ${OUTFILE} done"
echo ""
exit 0

