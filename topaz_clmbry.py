#

#Nick: run example with:
#python topaz_clmbry.py --incremental barents 20191201 20191203 /lustre/storeB/project/fou/hi/arktis2030/metroms_run/barents-2.5km/barents_clm.nc fimex_barents-2.5km_topaz.cfg TOPAZ2ROMS_varnames.ncml ./test/

# requires "module load romstools/1.1" on the PPI

# The overall goal of this script is to fetch climatology and boundary condition data for the
# Barents-2.5km ROMS+CICE coupled model. The script fetches climatology data from TOPAZ model
# output at thredds.met.no and interpolates it onto the Barents-2.5km grid, then outputs the
# data in clim-files and bry-files (boundary conditions exctracted from the climatology. The
# data gathering is done for a user specified range of dates and outputed in netcdf files in
# a user-defined directory. See fimex_barents-2.5km_topaz.cfg for what variables are retrieved.

# ===================================== IMPORTS ======================================
import os
import glob
import shutil
import argparse
import subprocess
import re
from copy import copy
from datetime import datetime, timedelta
from subprocess import call
from netCDF4 import Dataset

# =============================== FUNCTION DEFINITIONS ===============================
def cmd(command):
    """Function runs provided command in the system shell.

    Args:
        command (string) : Command to be executed in shell
    Returns:
        result (integer) : Shell returned status of command
    """
    print("> " + command)
    result = subprocess.call(command, shell = True)

    if result != 0:
        print("Command failed: %d" % result)
    else:
        return result

def replace_line(filename, pattern, replacer, precursor="", replace_all=False):
    """
    Function that searches a text file line-by-line and replaces a line with user speified string
    if a user provided regex pattern is matched on that line. The search only starts after an
    optional precursor string is found in the file.

    Args:
        filename (str)     : Filename of text file
        pattern (str)      : Regular expression pattern to search for
        replacer (str)     : String to replace on line where regex was matched
        precursor (str)    : Optional string to limit the regex search start from
                             the line where the precursor was found. If not provided,
                             search starts from the top of file.
        replace_all (bool) : If True, replace all lines where regex matched (after
                             optional precursor), otherwise replace only first occurence
    """
    tmp_file = filename + "_tmp"
    begin_search = False
    done_replace = False

    with open(tmp_file, "w") as new_file:
        with open(filename, "r") as old_file:
            for line in old_file:
                if precursor in line:
                    begin_search=True  # only look for line after precursor

                # if precursor found, no
                if begin_search and not done_replace and re.match(pattern, line):
                    new_file.write(replacer)

                    if not replace_all:
                        done_replace = True
                else:
                    new_file.write(line)  # otherwise just copy existing line

    shutil.move(tmp_file, filename)  # replace original with newly written

def valid_date(date_string):
    """
    Function that validates if a date string is of a certain format
    and returns a datetime object corresponding to that string.

    Args:
        date_string (str) : String of a date to be checked
    Returns:
        datetime (datetime) : Datetime object created from input string
    """
    try:
        return datetime.strptime(date_string, '%Y%m%d')
    except ValueError:
        message = 'Not a valid date: {0}.'.format(date_string)
        raise argparse.ArgumentTypeError(message)

def valid_path(path):
    """
    Function that validates if path exists.

    Args:
        path (str) : Suggested file path
    Returns:
        path (str) : If path exists
    """
    if os.path.exists(path):
        return path
    else:
        raise argparse.ArgumentTypeError("Invalid path {}!".format(path))

# ==================== PARSING COMMAND LINE ARGS FROM USER =======================
parser = argparse.ArgumentParser(description='Download TOPAZ data from thredds.met.no')
parser.add_argument('model_name', type=str, help='name of model (used for output filenames)')
parser.add_argument('start_date', type=valid_date, help='first date YYYYMMDD of data (inclusive)')
parser.add_argument('end_date', type=valid_date, help='final date YYYYMMDD of data (inclusive)')
parser.add_argument('grid_file', type=valid_path, help='filename of ROMS grid file')
#parser.add_argument('cfg_file', type=valid_path, help='filename of fimex config file for interpolation')
parser.add_argument('cfg_file', type=str, help='filename of fimex config file for interpolation')
parser.add_argument('ncml_file', type=valid_path, help='filename of fimex config file for variable names')
parser.add_argument('work_dir', type=valid_path, help='directory to store data')
parser.add_argument('--incremental', dest='incremental', action='store_true', help='fetch data in daily netcdf files')
parser.set_defaults(incremental=False)
args = parser.parse_args()

# =============================== STATIC VARIABLES ===============================
# define filepaths for output files
clmfile = os.path.join(args.work_dir, '{}_clm.nc'.format(args.model_name))
bryfile = os.path.join(args.work_dir, '{}_bry.nc'.format(args.model_name))
tmp_topaz2 = os.path.join(args.work_dir, 'tmp_topaz_interpolated.nc')

# vertical configuration parameters
theta_s, theta_b = 6., 0.3
Tcline = 100.
vtrans, vstretch = 2, 4
roms_nlevs = 42

# =============================== DATA GATHERING ===============================
# replace whatever output file is in cfg file with a tmp filename
#replace_line(args.cfg_file, "file ?= ?.+", "file={}\n".format(tmp_topaz2), precursor="[output]")

print('\nFetching TOPAZ data from {} to {} and storing data in {}...\n'.format(
    args.start_date, args.end_date, args.work_dir))

files_to_be_produced = list()
d = copy(args.start_date)

# main part of script where TOPAZ data is fetched from thredds
while d <= args.end_date:
    if args.incremental:  # only one day at the time
        start = d.strftime('%Y-%m-%d')
        stop = start
    else:                 # the entire period
        start = args.start_date.strftime('%Y-%m-%d')
        stop = args.end_date.strftime('%Y-%m-%d')

    cmd('cp fimex_barents-2.5km_topaz.cfg.base '+args.cfg_file)
    s_tStamp = d.strftime('%Y%m%d')
    #newF = 'https://thredds.met.no/thredds/dodsC/topaz/dataset-topaz4-arc-myoceanv2-{0}'.format(s_tStamp)
    #newF = 'https://thredds.met.no/thredds/dodsC/topaz/dataset-topaz4-arc-myoceanv2-be'
    newF = '/lustre/storeB/project/copernicus/sea/metnotopaz4_2/arctic/mersea-class1/{0}_dm-metno-MODEL-topaz4-ARC-b{0}-fv02.0.nc'.format(s_tStamp)
    replace_line(args.cfg_file, "file ?= ?.+", "file={}\n".format(newF), precursor="[input]")
    newF = args.work_dir+'topaz2barents_'+s_tStamp+'.nc'
    replace_line(args.cfg_file, "file ?= ?.+", "file={}\n".format(newF), precursor="[output]")
    
    cmd('fimex -c {} --extract.reduceTime.start {} --extract.reduceTime.end {} --input.config {}'.format(args.cfg_file, start, stop, args.ncml_file))

    tmp_topaz2 = newF
    if os.path.exists(tmp_topaz2):
        print('Successfull horizontal interpolation!')
    else:
        raise RuntimeError('Failed horizontal interpolation! Aborting...')

    tpz_nlevs = Dataset(tmp_topaz2).dimensions['depth'].size
    tmpnc_filepath = os.path.join(args.work_dir, 'tmp.nc')
    cmd('cp {} {}'.format(tmp_topaz2, tmpnc_filepath))
    fid = Dataset(tmp_topaz2)

    for key in fid.variables.keys():
        if '_FillValue' in fid.variables[key].ncattrs():
            cmd('ncatted -O -h -a _FillValue,{},m,f,0.0 {}'.format(key, tmpnc_filepath))
        else:
            print('Skipping removal of _FillValue for '+key)
    fid.close()

    cmd('mv {} {}'.format(tmpnc_filepath, tmp_topaz2))
    
    if args.incremental:
      d = d + timedelta(days=1)  # advance to the next day
    
    '''
    cmd('roms2roms<< EOF\n%s\n%i 0.0 0.0 0.0\n%s\n%s\n%s\n0 1\n%s\n%i %1.1f %1.1f %1.1f\n%s\n%s\n%i %i\nEOF' % (tmp_topaz2, tpz_nlevs, tmp_topaz2, tmp_topaz2, tmp_topaz2, args.grid_file, roms_nlevs, theta_s, theta_b, Tcline, args.work_dir+'/', args.model_name, vtrans, vstretch))
    fid = Dataset(clmfile, 'r+')
    ctime = fid.variables['clim_time']
    ctime[:] = ctime[:]/24. - 7305
    fid.close()

    cmd('calc_uvbar << EOF\n%s\n%s\n%1.1f %1.1f %1.1f %i %i\nEOF' % (clmfile, args.grid_file, Tcline, theta_b, theta_s, vtrans, vstretch))
    cmd('bry_from_clim<< EOF\n'+clmfile+'\n'+bryfile+'\nEOF')

    if args.incremental:
        clm_wdatetail = clmfile.replace('.nc', '_{}.nc'.format(d.strftime('%Y%m%d')))
        bry_wdatetail = bryfile.replace('.nc', '_{}.nc'.format(d.strftime('%Y%m%d')))
        files_to_be_produced.append(clm_wdatetail)
        files_to_be_produced.append(bry_wdatetail)
        cmd('mv {} {}'.format(clmfile, clm_wdatetail))  # add datestamp to filename
        cmd('mv {} {}'.format(bryfile, bry_wdatetail))  # add datestamp to filename
    else:
        files_to_be_produced.append(clmfile)
        files_to_be_produced.append(bryfile)
        break

# =============================== LOGISTICS AND CLEAN-UP ===============================
cmd('rm {}'.format(tmp_topaz2))  # remove intermediate interpolation file

# message to user if files were generated successfully
print('-------------------------------------------------------------------------')
print('Generated output data files:')
failed_files = list()

for filename in files_to_be_produced:
    if os.path.exists(filename):
        print(filename)

    else:
        failed_files.append(filename)

if len(failed_files) != 0:
    print('Something went wrong along the way and some files were not produced:')

    for filename in failed_files:
        print(filename)

'''

