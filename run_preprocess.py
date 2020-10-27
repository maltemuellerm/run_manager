import os
import datetime as dt

import helper_nml

def setup_case_ww3(ICE_RUNDIR, t0, tf):
  ymd0 = t0.strftime('%Y%m%d');
  
  os.system('ln -sf /home/sm_erith/WW3/work_cice/*ww3 {0}'.format(ICE_RUNDIR) )
  os.system('ln -sf /home/sm_erith/WW3/work_cice/ww3_shel {0}'.format(ICE_RUNDIR))
  #os.system('cp -f /home/sm_erith/WW3/work_cice/ww3_shel.inp {0}'.format(ICE_RUNDIR))
  os.system('cp -f /home/sm_erith/WW3/work_cice/ww3_strt.inp {0}'.format(ICE_RUNDIR))
  #ln -sf /nobackup/forsk/sm_nicsz/for_barents/ww3WindForce.nc ${ICE_RUNDIR}
  os.system('ln -sf /home/sm_erith/WW3/wind_forcing_Dec2019/10mWind_Forcing_WW3_{1}.nc {0}/ww3WindForce.nc'.format(ICE_RUNDIR,ymd0) )
  os.system('ln -sf /home/sm_erith/WW3/boundary_files_Dec2019/nest.{1}.ww3 {0}/nest.ww3'.format(ICE_RUNDIR,ymd0))
  os.system('cp -f /nobackup/forsk/sm_nicsz/for_barents/ww3_shel.inp {0}'.format(ICE_RUNDIR))
  os.system('rm -f {0}/out_grd.ww3'.format(ICE_RUNDIR))
  if (os.path.exists('{0}/restart001.ww3'.format(ICE_RUNDIR)) ):
    os.system('ln -sf {0}/restart001.ww3 {0}/restart.ww3'.format(ICE_RUNDIR))

  #
  s = helper_nml.form_ww3_shel(t0, tf)
  with open('ww3_shel.inp', 'w') as f:
    f.write(s)
  os.system('mv ww3_shel.inp {0}'.format(ICE_RUNDIR))

def setup_case_cice(ICE_RUNDIR, t0, deltat = 120.0, npt = 720, restart_fsd = '.false.'):
  
  ymd0 = t0.strftime('%Y%m%d');
  
  os.system('ln -sf /nobackup/forsk/sm_nicsz/for_barents/TOPAZ2CICE_{1}_1D.nc {0}/oceanForce.nc'.format(ICE_RUNDIR,ymd0) )
  os.system('ln -sf /nobackup/forsk/sm_nicsz/for_barents/AA2CICE_{1}_1Hr.nc {0}/atmoForce.nc'.format(ICE_RUNDIR,ymd0) )
  os.system('ln -sf /nobackup/forsk/sm_nicsz/for_barents/cice*.nc {0}'.format(ICE_RUNDIR) )
  
  tStamp = t0.strftime('%Y-%m-%d-00000'); fRestart_cice = 'restart/iced.{tStamp}.nc'.format(tStamp=tStamp)
  if (not os.path.exists('{0}/{1}'.format(ICE_RUNDIR,fRestart_cice)) ):
    os.system('ln -sf /nobackup/forsk/sm_nicsz/for_barents/iced.2020-01-23-00000.nc {0}/{1}'.format(ICE_RUNDIR, fRestart_cice) )

  s = helper_nml.form_ice_in(t0, deltat, npt, restart_fsd)
  with open('ice_in', 'w') as f:
    f.write(s)
  os.system('mv ice_in {0}'.format(ICE_RUNDIR))

def setup_case_oasis(ICE_RUNDIR, t0, tf):
  os.system('ln -sf /nobackup/forsk/sm_nicsz/for_barents/rmp*.nc {0}'.format(ICE_RUNDIR) ) #remapping and grid files for oasis
  if (not os.path.exists('{0}/restart_ice2wave.nc'.format(ICE_RUNDIR)) ):
    os.system('ln -sf /nobackup/forsk/sm_nicsz/for_barents/restart*.nc {0}'.format(ICE_RUNDIR) ) #coupling restarts when w/lag
  #t0 = dt.datetime(2019,12,6); tf = dt.datetime(2019,12,7)
  tRun = (tf-t0).total_seconds(); tRun=int(tRun)
  s = helper_nml.form_oasis_namcouple(tRun)
  with open('namcouple', 'w') as f:
    f.write(s)
  os.system('mv namcouple {0}'.format(ICE_RUNDIR))
#

def setup_case(ICE_RUNDIR, t0, tf, deltat = 120.0, npt = 720, restart_fsd = '.false.'):
  
  setup_case_ww3(ICE_RUNDIR, t0, tf)
  
  setup_case_cice(ICE_RUNDIR, t0, deltat = deltat, npt = npt, restart_fsd = restart_fsd)
  
  setup_case_oasis(ICE_RUNDIR, t0, tf)

def demo_setup_case():
  ICE_RUNDIR='/nobackup/forsk/sm_nicsz/CICE/CICE_RUNS/ciceFSD_WW3.cycle.t00/'
  t0 = dt.datetime(2019,12,6); tf = dt.datetime(2019,12,7)
  setup_case(ICE_RUNDIR, t0, tf)

if __name__=='__main__':
  demo_setup_case()



