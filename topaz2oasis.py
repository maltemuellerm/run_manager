import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import glob
import datetime as dt

'''
First put TOPAZ files on AA grid:
/lustre/storeB/users/nicholass/get_topaz/on_ppi$ python topaz_clmbry.py --incremental barents 20191209 20191220 barents_clm.nc fimex_barents-2.5km_topaz.cfg TOPAZ2ROMS_varnames.ncml ./201912/

Convert to OASIS format:

'''

def make_oasisFields(iTime, data, nx=739,ny=949, iLev=0):
  '''
  Given a fimex(topaz->barents), make fields for oasis input of:
  ['sst_ext', 'sss_ext','uocn_ext','vocn_ext','hmix_ext','qdp_ext','aice_bry','vice_bry','vsno_bry']
  '''
  vals_out = {}
  
  vals = data.variables['temperature'][iTime,iLev, :,:]
  vals = np.ravel(vals)
  vals_out['sst_ext'] = vals
  
  vals = data.variables['salinity'][iTime,iLev, :,:]
  vals = np.ravel(vals)
  vals_out['sss_ext'] = vals
  
  vals = data.variables['u'][iTime,iLev, :,:]
  vals = np.ravel(vals)
  vals_out['uocn_ext'] = vals
  
  vals = data.variables['v'][iTime,iLev, :,:]
  vals = np.ravel(vals)
  vals_out['vocn_ext'] = vals
  
  vals = data.variables['mlp'][iTime, :,:]
  vals = np.ravel(vals)
  vals_out['hmix_ext'] = vals
  
  vals_out['qdp_ext'] = vals*0 #TODO: heat flux @ bottom of mixed layer
  
  vals = data.variables['aice'][iTime, :,:]
  vals = np.ravel(vals)
  vals_out['aice_bry'] = vals
  aice = vals
  
  vals = data.variables['hice'][iTime, :,:]
  vals = np.ravel(vals)
  vals = vals*aice
  vals_out['vice_bry'] = vals
  
  vals = data.variables['hsnow'][iTime, :,:]
  vals = np.ravel(vals)
  vals = vals*aice
  vals_out['vsno_bry'] = vals
  
  return vals_out

def demo():
  fDir = '/scratch/topaz/'
  fList = sorted(glob.glob(fDir+'topaz2barents*nc')); print fList
  
  for f in fList:
    with netCDF4.Dataset(f,'r') as data:
      vals_iTime = make_oasisFields(0, data)
      print vals_iTime

def demo_writeFile(ymd='20191206'):
  #User inputs --------------
  nx=739; ny=949; nt = 3; #for output file
  delt=24*3600; ntEachFile = 1 #for input files 
  
  #input files
  fDir = './'
  fList = sorted(glob.glob( fDir+'topaz2barents_{0}.nc'.format(ymd) )); print fList
  
  #output files
  #fNameOut = fDir+'TOPAZ_forcing.test.nc'
  fNameOut = fDir+'TOPAZ2CICE_{0}_1D.nc'.format(ymd)
  #end user inputs ------------
  
  #create and dimension output file
  dataOut = netCDF4.Dataset(fNameOut,'w')
  dataOut.createDimension('x', nx*ny); dataOut.createDimension('y', 1); dataOut.createDimension('time', nt);
  
  #from AA output values from file -> oasis input file
  nf = len(fList); iTimeGlobal=0
  for iFile in xrange(nf):
    data = netCDF4.Dataset(fList[iFile],'r')
    for iTime in xrange(ntEachFile):
      vals_iTime = make_oasisFields(iTime, data)
      for key in vals_iTime:
        print key, vals_iTime[key]
      
      if iTimeGlobal==0: #create variables in netcdf file
        dataOut.createVariable('time','i8',('time'))
        for key in vals_iTime:
          dataOut.createVariable(key,'f8',('time','y','x'), fill_value=0.0)
      #store values in file
      dataOut.variables['time'][iTimeGlobal] = iTimeGlobal*delt
      for key in vals_iTime:
        dataOut.variables[key][iTimeGlobal,0,:] = vals_iTime[key]
      iTimeGlobal = iTimeGlobal+1
    data.close()
  
  dataOut.close()

def make_TOPAZ2CICE_daily():
  t0 = dt.datetime(2019,12,6)
  delt = dt.timedelta(days=1)
  tf = dt.datetime(2019,12,10) #tf = dt.datetime(2019,12,31)
  while (t0<=tf):
    s = t0.strftime('%Y%m%d')
    demo_writeFile(ymd=s)
    t0 = t0+delt
  
if __name__=='__main__':
  #demo()
  #demo_writeFile()
  make_TOPAZ2CICE_daily()
