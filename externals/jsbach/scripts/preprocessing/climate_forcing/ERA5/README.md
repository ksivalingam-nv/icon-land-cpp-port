ERA5 data processing script

Script to process ERA5 data so it can be correctly read by ICON-Land in standalone mode.

When downloading ERA5 data, all variables must be saved in the same file. The Name of the file must follow the pattern `${year}_${month}_hourly.nc`, with the year and month following the formats YYYY, and MM, respectively.

Each year can be analyzed independently, thus a job per year can be submitted.