# -*- coding: utf-8 -*-

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import glob
import datetime as dt

'''
For CICE, I use:
strau:strav:wind
frain:fsnow:Tair:qair:Pair:rhoa
swdvdr:swdvdf:swdidr:swdidf:lwd
From your AA files (say Dec2019_CoupVerification_AA_2019121009.nc for concreteness), I see:
For strau:strav:wind (zonal/meridional wind stresses and wind speed)...x_wind_10m and y_wind_10m or u/v (in WW3) winds. I can calculate an approximate wind stress from the (average?) winds given a drag coefficient. Wind speed is the magnitude of the components.
For frain:snow (instantaneous rain/snow)...rainfall_amount and snowfall_amount.These look like they accumulate over the forecast (despite the metadata). At the initial time, they aren't zero and the spatial field has a physical structure. 

For Tair: air_temperature_2mFor qair: specific_humidity_2mFor Pair: surface_air_pressureFor rhoa...I can do ideal gas with Tair and Pair
For swdvdr:swdvdf:swdidr:swdidf (the downwelling shortwave at the surface)...integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time.Values at initial time are 1.e37 type garbage. I think I just replace them with 0.
For lwd (downwelling longwave at the surface)...integral_of_surface_downwelling_longwave_flux_in_air_wrt_timeValues at initial time are nonzero and both positive and negative. The spatial field has a physical structure
'''

def calc_angle_global2grid(lat2d, lon2d, r=1.,plotAngle=False):
  #Input double longitude(y, x) in radians.
  #Return the static angle used in rotating grid-relative winds to zonal/meridional so
  #angleGlobal = angleGrid+atan2(yvel/xvel)
  
  #Following from Eivind (SUV), we just calculate the orientation from
  #the discretely sampled grid
  #...I'll just look at the neighbor x+1 and copy values for right boundary
  
  latp1 = lat2d.copy(); lonp1 = lon2d.copy()
  latp1[:,0:-1] = lat2d[:,1:]; lonp1[:,0:-1] = lon2d[:,1:];
  
  #tilt of grid
  distLat = r*(latp1-lat2d)
  distLon = r*(lonp1-lon2d)*np.cos(lat2d)
  ang = np.arctan2(distLat,distLon) #Element-wise arc tangent of x1/x2 choosing the quadrant correctly
  #take care of right boundary
  ang[:,-1] = ang[:,-2]
  
  if plotAngle:
    plt.figure()
    #plt.pcolormesh(lon2d,lat2d,ang*180./np.pi)
    plt.pcolormesh(ang*180./np.pi)
    plt.colorbar()
    plt.show()
  
  return ang #angle zonal --> x-dir
  

def make_oasisFields(iTime, data, nx=739,ny=949, iLev=0):
  #Input data is from atmosphere.
  #Input wind is from separate file (waves)
  
  vals_out = {}
  
  #instantaneous values ---------------
  vals = data.variables['air_temperature_2m'][iTime,iLev, :,:]; Tair = vals
  vals = np.ravel(vals)
  vals_out['Tair'] = vals #in K
  
  vals = data.variables['surface_air_pressure'][iTime,iLev, :,:]; Pair =  vals
  vals = np.ravel(vals)
  vals_out['Pair'] = vals #in Pa
  
  vals = data.variables['specific_humidity_2m'][iTime,iLev, :,:]; qair =  vals
  vals = np.ravel(vals)
  vals_out['qair'] = vals #in kg/kg
  
  #"derived" instantaneous values
  vals = Pair / (Tair * 287.058); rhoa = vals
  vals = np.ravel(vals)
  vals_out['rhoa'] = vals #in kg
  
  Cd = 1.e-3
  vals = data.variables['x_wind_10m'][iTime,iLev, :,:]; xvel =  vals
  vals = data.variables['y_wind_10m'][iTime,iLev, :,:]; yvel =  vals
  wind = np.sqrt(xvel*xvel+yvel*yvel)
  vals = np.ravel(wind)
  vals_out['wind'] = vals #in m/s
  
  lat2d = data.variables['latitude'][:,:]*np.pi/180.
  lon2d = data.variables['longitude'][:,:]*np.pi/180.
  angle = calc_angle_global2grid(lat2d, lon2d, r=1.,plotAngle=False)
  
  angleWind_xy = np.arctan2(yvel,xvel)
  uvel = wind*np.cos(angle+angleWind_xy)
  vvel = wind*np.sin(angle+angleWind_xy)
  vals = np.ravel(uvel)
  vals_out['u'] = vals #in m/s
  vals = np.ravel(vvel)
  vals_out['v'] = vals #in m/s
  
  vals = .5*rhoa*Cd*uvel*uvel
  vals = np.ravel(vals)
  vals_out['strau'] = vals #in TODO
  vals = .5*rhoa*Cd*vvel*vvel
  vals = np.ravel(vals)
  vals_out['strav'] = vals #in TODO
  
  #accumulated values ------------------
  #...frain:fsnow:swd:lwd
  #setting values at time 0 accumulated over 0->1 in file at time 1-time 0
  delt = data.variables['time'][iTime+1]-data.variables['time'][iTime]
  
  key = 'rainfall_amount'
  vals1 = data.variables[key][iTime+1,iLev, :,:]; vals0 = data.variables[key][iTime,iLev, :,:]
  if iTime == 0: #unclear what value at initial time represents
    vals0 = 0
  vals = (vals1-vals0)/delt
  vals = np.ravel(vals)
  vals_out['frain'] = vals #in kg/m^2 s
  
  key = 'snowfall_amount'
  vals1 = data.variables[key][iTime+1,iLev, :,:]; vals0 = data.variables[key][iTime,iLev, :,:]
  if iTime == 0: #unclear what value at initial time represents
    vals0 = 0
  vals = (vals1-vals0)/delt
  vals = np.ravel(vals)
  vals_out['fsnow'] = vals #in kg/m^2 s
  
  key = 'integral_of_surface_downwelling_longwave_flux_in_air_wrt_time'
  vals1 = data.variables[key][iTime+1,iLev, :,:]; vals0 = data.variables[key][iTime,iLev, :,:]
  if iTime == 0: #unclear what value at initial time represents
    vals0 = 0
  vals = (vals1-vals0)/delt
  vals = np.ravel(vals)
  vals_out['lwd'] = vals #in W/m2
  
  #single swdown is partitioned into 4 streams based on cice ice_forcing.F90 defaults
  key = 'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time'
  vals1 = data.variables[key][iTime+1,iLev, :,:]; vals0 = data.variables[key][iTime,iLev, :,:]
  if iTime == 0: #unclear what value at initial time represents
    vals0 = 0
  swd = (vals1-vals0)/delt
  vals = np.ravel(swd)
  vals_out['swdvdr'] = vals*.28 #in W/m2
  vals_out['swdvdf'] = vals*.24 #in W/m2
  vals_out['swdidr'] = vals*.31 #in W/m2
  vals_out['swdidf'] = vals*.17 #in W/m2
  
  #
  return vals_out

def demo():
  fDir = '/scratch/topaz/'
  fList = sorted(glob.glob(fDir+'Dec2019*_AA_*nc')); print fList
  
  nf = len(fList)
  for iFile in xrange(nf):
    data = netCDF4.Dataset(fList[iFile],'r')
    for iTime in xrange(3):
      vals_iTime = make_oasisFields(iTime, data)
      for key in vals_iTime:
        print key, vals_iTime[key]
      if True:
        nx=739; ny=949
        lat2d = data.variables['latitude'][:,:]
        lon2d = data.variables['longitude'][:,:]
        plt.figure()
        vals = np.reshape(vals_iTime['u'],(ny,nx),order='F')
        #plt.pcolormesh(lon2d,lat2d,u2d,vmin=-10,vmax=10)
        plt.pcolormesh(vals,vmin=-10,vmax=10)
        plt.colorbar()
        plt.title('u')
        
        plt.figure()
        vals = np.reshape(vals_iTime['v'],(ny,nx),order='F')
        #plt.pcolormesh(lon2d,lat2d,u2d,vmin=-10,vmax=10)
        plt.pcolormesh(vals,vmin=-10,vmax=10)
        plt.colorbar()
        plt.title('v')
        
        plt.show()
      
    data.close()
    
def demo_writeFile(ymd='20191206'):
  #User inputs --------------
  nx=739; ny=949; nt = 24; #for output file
  delt=3600; ntEachFile = 3 #for input files 
  
  #input files
  fDir = '/lustre/storeB/users/nicholass/get_topaz/on_ppi/201912/'
  #fList = sorted(glob.glob(fDir+'Dec2019*_AA_*nc')); print fList
  fList = sorted(glob.glob( fDir+'Dec2019*_AA_{0}*nc'.format(ymd) )); print fList
  
  #output files
  #fNameOut = fDir+'AA_forcing.test.nc'
  fNameOut = fDir+'AA2CICE_{0}_1Hr.nc'.format(ymd)
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
          dataOut.createVariable(key,'f8',('time','y','x'))
      #store values in file
      dataOut.variables['time'][iTimeGlobal] = iTimeGlobal*delt
      for key in vals_iTime:
        dataOut.variables[key][iTimeGlobal,0,:] = vals_iTime[key]
      iTimeGlobal = iTimeGlobal+1
    data.close()
  
  dataOut.close()

def make_oasisFields_waves(iTime, data, nx=737,ny=947, iLev=0):
  #Input data is from atmosphere.
  
  vals_out = {}
  
  vals = data.variables['x_wind_10m'][iTime,iLev, 1:-1,1:-1]; xvel =  vals
  vals = data.variables['y_wind_10m'][iTime,iLev, 1:-1,1:-1]; yvel =  vals
  wind = np.sqrt(xvel*xvel+yvel*yvel)
  vals = np.ravel(wind)
  vals_out['wind'] = vals #in m/s
  
  lat2d = data.variables['latitude'][1:-1,1:-1]*np.pi/180.
  lon2d = data.variables['longitude'][1:-1,1:-1]*np.pi/180.
  angle = calc_angle_global2grid(lat2d, lon2d, r=1.,plotAngle=False)
  
  angleWind_xy = np.arctan2(yvel,xvel)
  uvel = wind*np.cos(angle+angleWind_xy)
  vvel = wind*np.sin(angle+angleWind_xy)
  vals = np.ravel(uvel)
  vals_out['WW3__U10'] = vals #in m/s
  vals = np.ravel(vvel)
  vals_out['WW3__V10'] = vals #in m/s
  
  return vals_out

def demo_writeFile_wave():
  #User inputs --------------
  nx=737; ny=947; nt = 24; #for output file
  delt=3600; ntEachFile = 3 #for input files 
  
  #input files
  fDir = '/lustre/storeB/users/nicholass/get_topaz/on_ppi/201912/'
  fList = sorted(glob.glob(fDir+'Dec2019*_AA_*nc')); print fList
  
  #output files
  fNameOut = fDir+'AA_forcing_WW3.test.nc'
  #end user inputs ------------
  
  #create and dimension output file
  dataOut = netCDF4.Dataset(fNameOut,'w')
  dataOut.createDimension('x', nx*ny); dataOut.createDimension('y', 1); dataOut.createDimension('time', nt);
  
  #from AA output values from file -> oasis input file
  nf = len(fList); iTimeGlobal=0
  for iFile in xrange(nf):
    data = netCDF4.Dataset(fList[iFile],'r')
    for iTime in xrange(ntEachFile):
      vals_iTime = make_oasisFields_waves(iTime, data)
      for key in vals_iTime:
        print key, vals_iTime[key]
      
      if iTimeGlobal==0: #create variables in netcdf file
        dataOut.createVariable('time','i8',('time'))
        for key in vals_iTime:
          dataOut.createVariable(key,'f8',('time','y','x'))
      #store values in file
      dataOut.variables['time'][iTimeGlobal] = iTimeGlobal*delt
      for key in vals_iTime:
        dataOut.variables[key][iTimeGlobal,0,:] = vals_iTime[key]
      iTimeGlobal = iTimeGlobal+1
    data.close()
  
  dataOut.close()

def make_AA2CICE_daily():
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
  #demo_writeFile_wave()
  make_AA2CICE_daily()

