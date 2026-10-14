#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------



###############################################################################
# Script to process the ERA5 hourly data into a format that could be used by
# ICON-Land in stand alone mode.
#
# Name of the input files must follow the pattern ${yr}_${mon}_hourly.nc, with
# each file containing all required variables (see below). The year and month
# must provided with the formats YYYY, and MM, respectively.
#
# The script can be submitted independently to process individual years. A job
# for each year can be submitted.
###############################################################################

### Batch Queuing System is SLURM
#SBATCH --job-name=generate_forcing_era5_data
#SBATCH --output=generate_forcing_era5_data.o%j
#SBATCH --error=generate_forcing_era5_data.o%j
#SBATCH --partition=compute
#SBATCH --mem=200GB
#SBATCH --time=08:00:00
#SBATCH --account mj0143


module load cdo nco

# Directory with the input data
path_in=${1}
# Directory to save the data
path_out=${1}
# Pattern for the produced files
fout_name=$2 ### e.g. ERA5_R2B5
# Year to be processed
yr=$3
# File with the desired ICON grid
grid=$4 ### eg., icon_grid_0030_R02B05_G.nc

mkdir -p $path_out

#### Variables --> do not change order!!!
declare -a variables=(total_precipitation                       #0
            2m_temperature                                      #1
            10m_u_component_of_wind                             #2
            10m_v_component_of_wind                             #3
            surface_pressure                                    #4
            2m_dewpoint_temperature                             #5
            surface_solar_radiation_downwards                   #6
            surface_thermal_radiation_downwards                 #7
            )

data_list=""
for mon in $(seq -f "%02g" 1 12)
do
  echo $mon

  fin_name=${yr}_${mon}_hourly.nc
  fout=${path_out}${fout_name}_${yr}_${mon}.nc
  fouttog=${path_out}${fout_name}_${yr}.nc


  ### Temperature
  echo 'Min Temperature'
  cdo -L -subc,273.15 -daymin -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[1]}_${fin_name} ${fout}
  ncatted -hO -a history,global,d,, ${fout} ## Remove history
  ncrename -v t2m,tmin ${fout}
  ncatted -hO -a units,tmin,o,c,"degC" ${fout}
  ncatted -hO -a title,tmin,c,c,"min temperature" ${fout}

  echo 'Max Temperature'
  ttmp=${path_out}temp_sum_NH_${yr}_${mon}.nc
  cdo -L -subc,273.15 -daymax -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[1]}_${fin_name} ${ttmp}
  ncatted -hO -a history,global,d,, ${ttmp} ## Remove history
  ncrename -v t2m,tmax ${ttmp}
  ncatted -hO -a units,tmax,o,c,"degC" ${ttmp}
  ncatted -hO -a title,tmax,c,c,"max temperature" ${ttmp}
  ncks -hA -v tmax ${ttmp} ${fout}
  rm ${ttmp}

  #### Precipitation
  echo 'Total Precip'
  prectmp=${path_out}prec_sum_NH_${yr}_${mon}.nc
  cdo -L -mulc,0.2778 -dayavg -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[0]}_${fin_name} ${prectmp}
  ncrename -hO -v tp,precip ${prectmp}
  ncatted -hO -a units,precip,o,c,"kg m-2 s-1" ${prectmp}
  ncatted -hO -a long_name,precip,o,c,"Total precipitation" ${prectmp}
  ncwa -O -h -a expver ${prectmp} ${prectmp}
  ncks -hA -v precip ${prectmp} ${fout}
  rm -f ${prectmp} ${raintmp} ${raintmp2}

  echo 'Wind speed'
  utmp=${path_out}u10_NH_${yr}_${mon}.nc
  vtmp=${path_out}v10_NH_${yr}_${mon}.nc
  cdo -L -dayavg -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[2]}_${fin_name} ${utmp}
  cdo -L  -dayavg -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[3]}_${fin_name} ${vtmp}
  ncrename -hO -v u10,Wind_E ${utmp}
  ncrename -hO -v v10,Wind_N ${vtmp}
  ncks -hA ${utmp} ${vtmp}

  ncap2 -hO -s 'wspeed=sqrt(Wind_N^2+Wind_E^2)'  ${vtmp} ${vtmp}
  ncatted -hO -a units,wspeed,o,c,"m s-1" ${vtmp}
  ncatted -hO -a long_name,wspeed,o,c,"10m Wind-speed" ${vtmp}
  ncks -hO -x -v 'Wind_N' ${vtmp}  ${vtmp} # Delete var to reduce memory requirements
  ncks -hO -x -v 'Wind_E' ${vtmp}  ${vtmp} # Delete var to reduce memory requirements
  ncks -hA ${vtmp} ${fout}

  echo 'Sea level pressure'
  slptmp=${path_out}slp_tmp_NH_${yr}_${mon}.nc
  cdo -L -dayavg -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[4]}_${fin_name} ${slptmp}
  ncks -hA ${slptmp}  ${fout}
  ncks -hO -x -v sp ${fout} ${fout}
  ncrename -hO -v sp,PSurf ${slptmp}
  ncatted -hO -a long_name,PSurf,o,c,"Pressure" ${slptmp}

  echo 'Specific Humidity'
  qtmp=${path_out}qair_tmp_NH_${yr}_${mon}.nc
  cdo -L -dayavg -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[5]}_${fin_name} ${qtmp}
  ncap2 -hO -s 'Td=d2m-273.15' ${qtmp}  ${qtmp}
  ncks -hO -x -v 'd2m' ${qtmp}  ${qtmp} # Delete var to reduce memory requirements
  ncks -hA -v PSurf ${slptmp} ${qtmp}
  ncap2 -hO -s 'e=6.112*exp((17.67*Td)/(Td+243.5))'  ${qtmp} ${qtmp}
  ncks -hO -x -v 'Td' ${qtmp}  ${qtmp} # Delete var to reduce memory requirements
  ncap2 -hO -s 'qair=(0.622*e) /((PSurf/100.)-(0.378*e))'   ${qtmp} ${qtmp}
  ncks -hO -x -v 'PSurf' ${qtmp}  ${qtmp} # Delete var to reduce memory requirements
  ncatted -hO -a long_name,qair,o,c,"Air_Specific_Humidity" ${qtmp}
  ncatted -hO -a units,qair,o,c,"kg kg-1" ${qtmp}
  ncks -hA -v qair ${qtmp}  ${fout}

  echo "SW Radiation"
  frad=${path_out}rad_tmp_NH_${yr}_${mon}.nc
  cdo -L  -dayavg -divc,3600. -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[6]}_${fin_name}  ${frad}
  ncwa -O -h -a expver ${frad} ${frad}
  ncrename -hO -v ssrd,shortwave ${frad}
  ncatted -hO -a long_name,shortwave,o,c,"Incoming_Short_Wave_Radiation" ${frad}
  ncatted -hO -a units,shortwave,o,c,"W m-2" ${frad}
  ncatted -hO -a standard_name,,d,, ${frad}
  ncks -hA -v shortwave ${frad} ${fout}

  echo "LW radiation"
  flw=${path_out}lwr_tmp_NH_${yr}_${mon}.nc
  cdo -L -b F64 -dayavg -divc,3600. -sellonlatbox,-180,180,-90,90 -remapcon,${grid} ${path_in}/${variables[7]}_${fin_name}  ${flw}
  ncwa -O -h -a expver ${flw} ${flw}
  ncrename -hO -v strd,longwave ${flw}
  ncatted -hO -a long_name,longwave,o,c,"Incoming_Long_Wave_Radiation" ${flw}
  ncatted -hO -a units,longwave,o,c,"W m-2" ${flw}
  ncks -hA -v longwave ${flw} ${fout}


  cdo -L -b F64 -O -setmissval,1.e+20 ${fout} ${fout}.tmp #

  ## Fix after ERA5 update in 2025 (time axis broken)
  nccopy -k classic  ${fout}.tmp ${fout}
  ncrename -hO -d valid_time,time -v valid_time,time ${fout} ${fout}.tmp
  cdo -L -a setcalendar,proleptic_gregorian -settaxis,"${yr}-${mon}-01,00:00:00,1day" ${fout}.tmp ${fout}
  rm ${fout}.tmp

  date_today=$(date +'%Y-%m-%d')
  echo "Global attributes"
  ncatted -hO -a history,global,d,, ${fout}
  ncatted -hO -a delta_tstep_sec,global,c,f,86400.0  ${fout}
  ncatted -hO -a day0,global,c,f,1.0  ${fout}
  ncatted -hO -a month0,global,c,f,1.0  ${fout}
  ncatted -hO -a year0,global,c,f,${yr}.0  ${fout}
  ncatted -hO -a history,global,c,c,"Prepared by F.J. Cuesta-Valero ${date_today} (francisco_jose.cuesta_valero@uni-leipzig.de)" ${fout}


  rm -f ${flw}
  rm -f ${frad}
  rm -f ${qtmp}
  rm -f ${path_out}snow*${yr}_${mon}.nc
  rm -f ${slptmp}
  rm -f ${utmp}
  rm -f ${vtmp}
  rm -f ${path_out}rain*${yr}_${mon}.nc

  data_list="${data_list} ${fout}"

done

rm -f ${fouttog}
cdo -L -copy ${data_list} ${fouttog}


echo "*** Script run until the end ***"