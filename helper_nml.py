import datetime as dt

def form_ice_in(t0, deltat, npt, restart_fsd):
  s = '''&setup_nml
    days_per_year  = 365
    use_leap_years = .false.
    year_init      = {year_init}
    istep0         = {istep0}
    dt             = {deltat}
    npt            = {npt}
    ndtd           = 1
    runtype        = 'initial'
    ice_ic         = 'restart/iced.{tStamp}.nc'
    restart        = .true.
    restart_ext    = .true.
    use_restart_time = .false.
    restart_format = 'nc'
    lcdf64         = .true.
    numin          = 21
    numax          = 89
    restart_dir    = './restart/'
    restart_file   = 'iced'
    pointer_file   = './ice.restart_file'
    dumpfreq       = 'd'
    dumpfreq_n     = 1
    dump_last      = .false.
    bfbflag        = 'off'
    diagfreq       = 1
    diag_type      = 'stdout'
    diag_file      = 'ice_diag.d'
    print_global   = .true.
    print_points   = .true.
    latpnt(1)      =  90.
    lonpnt(1)      =  0.
    latpnt(2)      = -65.
    lonpnt(2)      = -45.
    dbug           = .false.
    histfreq       = 'm','h','x','x','x'
    histfreq_n     =  1 , 1 , 1 , 1 , 1
    hist_avg       = .true.
    history_dir    = './history/'
    history_file   = 'iceh'
    write_ic       = .true.
    incond_dir     = './history/'
    incond_file    = 'iceh_ic'
    version_name   = 'CICE_6.1.0'
/

&grid_nml
    grid_format  = 'nc'
    grid_type    = 'rectangular'
    grid_file    = '/nobackup/forsk/sm_nicsz/for_barents/cice.grid.nc'
    kmt_file     = '/nobackup/forsk/sm_nicsz/for_barents/cice.kmt.nc'
    bathymetry_file = '/nobackup/forsk/sm_nicsz/CICE//CICE_data/grid/barents/barents_bath.nc'
    use_bathymetry = .false.
    gridcpl_file = 'unknown_gridcpl_file'
    kcatbound    = 0
    dxrect         = 2.5e5
    dyrect         = 2.5e5
    close_boundaries = .false.
    ncat         = 5
    nfsd         = 12
    nilyr        = 7
    nslyr        = 1
    nblyr        = 7
/

&tracer_nml
    n_aero       = 1
    n_zaero      = 0
    n_algae      = 0
    n_doc        = 0
    n_dic        = 0
    n_don        = 0
    n_fed        = 0
    n_fep        = 0
    tr_iage      = .true.
    restart_age  = .false.
    tr_FY        = .true.
    restart_FY   = .false.
    tr_lvl       = .true.
    restart_lvl  = .false.
    tr_pond_cesm = .false.
    restart_pond_cesm = .false.
    tr_pond_topo = .false.
    restart_pond_topo = .false.
    tr_pond_lvl  = .true.
    restart_pond_lvl  = .false.
    tr_aero      = .false.
    restart_aero = .false.
    tr_fsd       = .true.
    restart_fsd  = {restart_fsd}
/

&thermo_nml
    kitd              = 1
    ktherm            = 2
    conduct           = 'bubbly'
    a_rapid_mode      =  0.5e-3
    Rac_rapid_mode    =    10.0
    aspect_rapid_mode =     1.0
    dSdt_slow_mode    = -5.0e-8
    phi_c_slow_mode   =    0.05
    phi_i_mushy       =    0.85
/

&dynamics_nml
    kdyn            = 1
    ndte            = 240
    revised_evp     = .false.
    kevp_kernel     = 0
    brlx            = 300.0
    arlx            = 300.0
    advection       = 'remap'
    kstrength       = 1
    krdg_partic     = 1
    krdg_redist     = 1
    mu_rdg          = 3
    Cf              = 17.
    Ktens           = 0.
    e_ratio         = 2.
    basalstress     = .false.
    k1              = 8.
    coriolis        = 'latitude'
    kridge          = 1
    ktransport      = 1
/

&shortwave_nml
    shortwave       = 'dEdd'
    albedo_type     = 'ccsm3'
    albicev         = 0.78
    albicei         = 0.36
    albsnowv        = 0.98
    albsnowi        = 0.70 
    ahmax           = 0.3
    R_ice           = 0.
    R_pnd           = 0.
    R_snw           = 1.5
    dT_mlt          = 1.5
    rsnw_mlt        = 1500.
    kalg            = 0.6
/

&ponds_nml
    hp1             = 0.01
    hs0             = 0.
    hs1             = 0.03
    dpscale         = 1.e-3
    frzpnd          = 'hlid'
    rfracmin        = 0.15
    rfracmax        = 1.
    pndaspect       = 0.8
/

&forcing_nml
    formdrag        = .false.
    atmbndy         = 'default'
    calc_strair     = .false.
    calc_Tsfc       = .true.
    highfreq        = .true.
    natmiter        = 5
    ustar_min       = 0.0005
    emissivity      = 0.95
    fbot_xfer_type  = 'constant'
    update_ocn_f    = .false.
    l_mpond_fresh   = .false.
    tfrz_option     = 'mushy'
    oceanmixed_ice  = .true.
    wave_spec_type  = 'profile'
    wave_spec_file  = 'unknown_wave_spec_file'
    nfreq           = 25
    restore_ice     = .true.
    restore_ocn     = .false.
    trestore        =  0
    precip_units    = 'mks'
    default_season  = 'winter'
    atm_data_type   = 'arome'
    ocn_data_type   = 'default'
    bgc_data_type   = 'default'
    fe_data_type    = 'default'
    ice_data_type   = 'default'
    fyear_init      = 2019
    ycycle          = 1
    atm_data_format = 'bin'
    atm_data_dir    = '/glade/u/home/tcraig/cice_data/'
    bgc_data_dir    = 'unknown_bgc_data_dir'
    ocn_data_format = 'bin'
    ocn_data_dir    = '/unknown_ocn_data_dir'
    oceanmixed_file = 'unknown_oceanmixed_file'
/

&domain_nml
    nprocs = 32
    nx_global         = 739
    ny_global         = 949
    block_size_x      = 20
    block_size_y      = 20
    max_blocks        = 56
    processor_shape   = 'slenderX2'
    distribution_type = 'roundrobin'
    distribution_wght = 'latitude'
    ew_boundary_type  = 'open'
    ns_boundary_type  = 'open'
    maskhalo_dyn      = .true.
    maskhalo_remap    = .true.
    maskhalo_bound    = .true.
/

&zbgc_nml
    tr_brine        = .true.
    restart_hbrine  = .false.
    tr_zaero        = .false.
    modal_aero      = .false.
    skl_bgc         = .false.
    z_tracers       = .false.
    dEdd_algae      = .false.
    solve_zbgc      = .false.
    bgc_flux_type   = 'Jin2006'
    restore_bgc     = .false.
    restart_bgc     = .false.
    scale_bgc       = .false.
    solve_zsal      = .false.
    restart_zsal    = .false.
    tr_bgc_Nit      = .true.
    tr_bgc_C        = .true.
    tr_bgc_chl      = .false.
    tr_bgc_Am       = .true.
    tr_bgc_Sil      = .true.
    tr_bgc_DMS      = .false.
    tr_bgc_PON      = .true.
    tr_bgc_hum      = .true.
    tr_bgc_DON      = .false.
    tr_bgc_Fe       = .true. 
    grid_o          = 0.006
    grid_o_t        = 0.006
    l_sk            = 0.024
    grid_oS         = 0.0
    l_skS           = 0.028
    phi_snow        = -0.3
    initbio_frac    = 0.8
    frazil_scav     = 0.8  
    ratio_Si2N_diatoms = 1.8                         
    ratio_Si2N_sp      = 0.0
    ratio_Si2N_phaeo   = 0.0
    ratio_S2N_diatoms  = 0.03  
    ratio_S2N_sp       = 0.03 
    ratio_S2N_phaeo    = 0.03
    ratio_Fe2C_diatoms = 0.0033
    ratio_Fe2C_sp      = 0.0033
    ratio_Fe2C_phaeo   = 0.1
    ratio_Fe2N_diatoms = 0.023 
    ratio_Fe2N_sp      = 0.023
    ratio_Fe2N_phaeo   = 0.7
    ratio_Fe2DON       = 0.023
    ratio_Fe2DOC_s     = 0.1
    ratio_Fe2DOC_l     = 0.033
    fr_resp            = 0.05
    tau_min            = 5200.0
    tau_max            = 173000.0
    algal_vel          = 0.0000000111
    R_dFe2dust         = 0.035
    dustFe_sol         = 0.005
    chlabs_diatoms     = 0.03
    chlabs_sp          = 0.01
    chlabs_phaeo       = 0.05
    alpha2max_low_diatoms = 0.8
    alpha2max_low_sp      = 0.67
    alpha2max_low_phaeo   = 0.67
    beta2max_diatoms   = 0.018
    beta2max_sp        = 0.0025
    beta2max_phaeo     = 0.01
    mu_max_diatoms     = 1.44
    mu_max_sp          = 0.851
    mu_max_phaeo       = 0.851
    grow_Tdep_diatoms  = 0.06
    grow_Tdep_sp       = 0.06
    grow_Tdep_phaeo    = 0.06
    fr_graze_diatoms   = 0.0
    fr_graze_sp        = 0.1
    fr_graze_phaeo     = 0.1
    mort_pre_diatoms   = 0.007
    mort_pre_sp        = 0.007
    mort_pre_phaeo     = 0.007
    mort_Tdep_diatoms  = 0.03
    mort_Tdep_sp       = 0.03
    mort_Tdep_phaeo    = 0.03
    k_exude_diatoms    = 0.0
    k_exude_sp         = 0.0
    k_exude_phaeo      = 0.0
    K_Nit_diatoms      = 1.0
    K_Nit_sp           = 1.0
    K_Nit_phaeo        = 1.0
    K_Am_diatoms       = 0.3
    K_Am_sp            = 0.3
    K_Am_phaeo         = 0.3
    K_Sil_diatoms      = 4.0
    K_Sil_sp           = 0.0
    K_Sil_phaeo        = 0.0
    K_Fe_diatoms       = 1.0
    K_Fe_sp            = 0.2
    K_Fe_phaeo         = 0.1
    f_don_protein      = 0.6
    kn_bac_protein     = 0.03
    f_don_Am_protein   = 0.25
    f_doc_s            = 0.4
    f_doc_l            = 0.4
    f_exude_s          = 1.0
    f_exude_l          = 1.0
    k_bac_s            = 0.03
    k_bac_l            = 0.03
    T_max              = 0.0
    fsal               = 1.0
    op_dep_min         = 0.1
    fr_graze_s         = 0.5
    fr_graze_e         = 0.5
    fr_mort2min        = 0.5
    fr_dFe             = 0.3
    k_nitrif           = 0.0
    t_iron_conv        = 3065.0
    max_loss           = 0.9
    max_dfe_doc1       = 0.2
    fr_resp_s          = 0.75
    y_sk_DMS           = 0.5
    t_sk_conv          = 3.0
    t_sk_ox            = 10.0
    algaltype_diatoms  = 0.0
    algaltype_sp       = 0.5
    algaltype_phaeo    = 0.5
    nitratetype        = -1.0
    ammoniumtype       = 1.0
    silicatetype       = -1.0
    dmspptype          = 0.5
    dmspdtype          = -1.0
    humtype            = 1.0
    doctype_s          = 0.5
    doctype_l          = 0.5
    dontype_protein    = 0.5
    fedtype_1          = 0.5
    feptype_1          = 0.5
    zaerotype_bc1      = 1.0
    zaerotype_bc2      = 1.0
    zaerotype_dust1    = 1.0
    zaerotype_dust2    = 1.0
    zaerotype_dust3    = 1.0
    zaerotype_dust4    = 1.0
    ratio_C2N_diatoms  = 7.0
    ratio_C2N_sp       = 7.0
    ratio_C2N_phaeo    = 7.0
    ratio_chl2N_diatoms= 2.1
    ratio_chl2N_sp     = 1.1
    ratio_chl2N_phaeo  = 0.84
    F_abs_chl_diatoms  = 2.0
    F_abs_chl_sp       = 4.0
    F_abs_chl_phaeo    = 5.0
    ratio_C2N_proteins = 7.0
/

&icefields_nml
    f_tmask        = .true.
    f_blkmask      = .true.
    f_tarea        = .true.
    f_uarea        = .true.
    f_dxt          = .false.
    f_dyt          = .false.
    f_dxu          = .false.
    f_dyu          = .false.
    f_HTN          = .false.
    f_HTE          = .false.
    f_ANGLE        = .true.
    f_ANGLET       = .true.
    f_NCAT         = .true.
    f_VGRDi        = .false.
    f_VGRDs        = .false.
    f_VGRDb        = .false.
    f_VGRDa        = .true.
    f_bounds       = .false.
    f_aice         = 'mh' 
    f_hi           = 'mh'
    f_hs           = 'mh' 
    f_Tsfc         = 'mh' 
    f_sice         = 'm' 
    f_uvel         = 'mh' 
    f_vvel         = 'mh' 
    f_uatm         = 'mh' 
    f_vatm         = 'mh' 
    f_fswdn        = 'mh' 
    f_flwdn        = 'mh'
    f_snowfrac     = 'x'
    f_snow         = 'h' 
    f_snow_ai      = 'm' 
    f_rain         = 'h' 
    f_rain_ai      = 'm' 
    f_sst          = 'mh' 
    f_sss          = 'm' 
    f_uocn         = 'm' 
    f_vocn         = 'm' 
    f_frzmlt       = 'm'
    f_fswfac       = 'm'
    f_fswint_ai    = 'm'
    f_fswabs       = 'x' 
    f_fswabs_ai    = 'm' 
    f_albsni       = 'm' 
    f_alvdr        = 'x'
    f_alidr        = 'x'
    f_alvdf        = 'x'
    f_alidf        = 'x'
    f_alvdr_ai     = 'x'
    f_alidr_ai     = 'x'
    f_alvdf_ai     = 'x'
    f_alidf_ai     = 'x'
    f_albice       = 'x'
    f_albsno       = 'x'
    f_albpnd       = 'x'
    f_coszen       = 'x'
    f_flat         = 'x' 
    f_flat_ai      = 'm' 
    f_fsens        = 'x' 
    f_fsens_ai     = 'm' 
    f_fswup        = 'h' 
    f_flwup        = 'h' 
    f_flwup_ai     = 'm' 
    f_evap         = 'x' 
    f_evap_ai      = 'm' 
    f_Tair         = 'mh' 
    f_Tref         = 'x' 
    f_Qref         = 'x'
    f_congel       = 'm' 
    f_frazil       = 'm' 
    f_snoice       = 'm' 
    f_dsnow        = 'x' 
    f_melts        = 'm'
    f_meltt        = 'm'
    f_meltb        = 'm'
    f_meltl        = 'm'
    f_fresh        = 'x'
    f_fresh_ai     = 'm'
    f_fsalt        = 'x'
    f_fsalt_ai     = 'm'
    f_fbot         = 'm'
    f_fhocn        = 'x' 
    f_fhocn_ai     = 'm' 
    f_fswthru      = 'x' 
    f_fswthru_ai   = 'm' 
    f_fsurf_ai     = 'x'
    f_fcondtop_ai  = 'x'
    f_fmeltt_ai    = 'x' 
    f_strairx      = 'mh' 
    f_strairy      = 'mh' 
    f_strtltx      = 'x' 
    f_strtlty      = 'x' 
    f_strcorx      = 'x' 
    f_strcory      = 'x' 
    f_strocnx      = 'x' 
    f_strocny      = 'x' 
    f_strintx      = 'x' 
    f_strinty      = 'x'
    f_taubx        = 'x'
    f_tauby        = 'x'
    f_strength     = 'm'
    f_divu         = 'm'
    f_shear        = 'm'
    f_sig1         = 'm' 
    f_sig2         = 'm' 
    f_sigP         = 'm' 
    f_dvidtt       = 'm' 
    f_dvidtd       = 'm' 
    f_daidtt       = 'mh'
    f_daidtd       = 'mh' 
    f_dagedtt      = 'm'
    f_dagedtd      = 'm' 
    f_mlt_onset    = 'm'
    f_frz_onset    = 'm'
    f_hisnap       = 'x'
    f_aisnap       = 'x'
    f_trsig        = 'm'
    f_icepresent   = 'm'
    f_iage         = 'm'
    f_FY           = 'x'
    f_aicen        = 'x'
    f_vicen        = 'x'
    f_vsnon        = 'x'
    f_snowfracn    = 'x'
    f_keffn_top    = 'x'
    f_Tinz         = 'x'
    f_Sinz         = 'x'
    f_Tsnz         = 'x'
    f_fsurfn_ai    = 'x'
    f_fcondtopn_ai = 'x'
    f_fmelttn_ai   = 'x'
    f_flatn_ai     = 'x'
    f_fsensn_ai     = 'x'
/

&icefields_mechred_nml
    f_alvl         = 'm'
    f_vlvl         = 'm'
    f_ardg         = 'm'
    f_vrdg         = 'm'
    f_dardg1dt     = 'x'
    f_dardg2dt     = 'x'
    f_dvirdgdt     = 'x'
    f_opening      = 'x'
    f_ardgn        = 'x'
    f_vrdgn        = 'x'
    f_dardg1ndt    = 'x'
    f_dardg2ndt    = 'x'
    f_dvirdgndt    = 'x'
    f_krdgn        = 'x'
    f_aparticn     = 'x'
    f_aredistn     = 'x'
    f_vredistn     = 'x'
    f_araftn       = 'x'
    f_vraftn       = 'x'
/

&icefields_pond_nml
    f_apondn       = 'x'
    f_apeffn       = 'x'
    f_hpondn       = 'x'
    f_apond        = 'm'
    f_hpond        = 'm'
    f_ipond        = 'm'
    f_apeff        = 'm'
    f_apond_ai     = 'm'
    f_hpond_ai     = 'm'
    f_ipond_ai     = 'm'
    f_apeff_ai     = 'm'
/

&icefields_bgc_nml
    f_faero_atm    = 'x'
    f_faero_ocn    = 'x'
    f_aero         = 'x'  
    f_fbio         = 'm'
    f_fbio_ai      = 'm'
    f_zaero        = 'x'
    f_bgc_S        = 'm'
    f_bgc_N        = 'm'
    f_bgc_C        = 'x'
    f_bgc_DOC      = 'm'
    f_bgc_DIC      = 'x'
    f_bgc_chl      = 'x'
    f_bgc_Nit      = 'm'
    f_bgc_Am       = 'm'
    f_bgc_Sil      = 'm'
    f_bgc_DMSPp    = 'x'
    f_bgc_DMSPd    = 'x'
    f_bgc_DMS      = 'x'
    f_bgc_DON      = 'x'  
    f_bgc_Fe       = 'm'  
    f_bgc_hum      = 'm'   
    f_bgc_PON      = 'm'
    f_bgc_ml       = 'm'
    f_upNO         = 'm'
    f_upNH         = 'm'
    f_bTin         = 'm'
    f_bphi         = 'm' 
    f_iDi          = 'm'
    f_iki          = 'm'
    f_fbri         = 'm'  
    f_hbri         = 'm'
    f_zfswin       = 'm'
    f_bionet       = 'm'
    f_biosnow      = 'm'
    f_grownet      = 'm'
    f_PPnet        = 'm'
    f_algalpeak    = 'm'
    f_zbgc_frac    = 'm'
/

&icefields_drag_nml
    f_drag         = 'x'
    f_Cdn_atm      = 'x'
    f_Cdn_ocn      = 'x'
/

&icefields_fsd_nml
    f_fsdrad       = 'mh'
    f_fsdperim     = 'mh'
    f_afsd         = 'mh'
    f_afsdn        = 'mh'
    f_dafsd_newi   = 'h'
    f_dafsd_latg   = 'h'
    f_dafsd_latm   = 'h'
    f_dafsd_wave   = 'h'
    f_dafsd_weld   = 'h'
    f_wave_sig_ht  = 'x'
    f_aice_ww      = 'x'
    f_diam_ww      = 'x'
    f_hice_ww      = 'x'
/
'''
  
  #calculate istep0
  istep0 = int( (t0-dt.datetime(t0.year,1,1)).total_seconds() / deltat)
  tStamp = t0.strftime('%Y-%m-%d-00000');
 
  return s.format(year_init=t0.year, istep0=istep0, deltat=deltat, npt=npt, restart_fsd=restart_fsd, tStamp=tStamp)

def form_ww3_shel(t0, tf):
  s = '''$ -------------------------------------------------------------------- $
$ WAVEWATCH III shell input file                                       $
$ -------------------------------------------------------------------- $
$ Include ice and mud parameters only if IC1/2/3/4 used :
   C F     Ice parameter 1
   F F     Ice parameter 2
   F F     Ice parameter 3
   F F     Ice parameter 4
   F F     Ice parameter 5
   F F     Mud parameter 1
   F F     Mud parameter 2
   F F     Mud parameter 3
   F F     Water levels
   F F     Currents
   C F     Winds
   C F     Ice concentrations
   F       Assimilation data : Mean parameters
   F       Assimilation data : 1-D spectra
   F       Assimilation data : 2-D spectra
$
$ - Starting and ending time in yyyymmdd hhmmss format.
   {ymd0} {hms0}
   {ymdf} {hmsf}
$
$ Define output data ------------------------------------------------- $
  2
$ ---------------------------------------------------------------
$ Type 1 : Fields of mean wave parameters
   {ymd0} {hms0}  3600  {ymdf} {hmsf}
$ (1) Forcing Fields
  T
$ DPT CUR WND AST WLV ICE IBG D50 IC1 IC5
  T   T   T   T   F   T   F   F   T   F
$ (2) Standard mean wave Parameters
  T
$ HS  LM  T02 T0M1 T01 FP DIR SPR DP
  T   T   T   T   T   T   T   T   T
$ (3) Frequency-dependent parameters
  T
$ E3D TH1MF STH1MF TH2M STH2M WN
  F   F   F   F   F   F
$ (4) Spectral Partition Parameters
  T
$ PHS PTP PLP PDIR PSPR PWS TWS PNR
  T   T   T   T   T   T   T   T
$ (5) Atmosphere-waves layer
  T
$ UST CHA CGE FAW TAW TWA WCC WCF WCH WCM
  T   T   T   T   T   T   F   F   F   F
$ (6) Wave-Ocean layer
  T
$ SXY TWO BHD FOC TUS USS P2S USF P2L TWI FIC
  T   T   F   T   T   F   F   F   F   T   T
$ (7) Wave-bottom layer
  T
$ ABR UBR BED FBB TBB
  F   F   F   F   F
$ (8) Spectrum parameters
  T
$ MSS MSC WL02 AXT AYT AXY
  F   F   F   F   F   F
$ (9) Numerical diagnostics
  T
$ DTD FC  CFX CFD CFK
  T   T   F   F   F
$ (10) User defined (NOEXTR flags needed)
  F
$ LBRK ZBRK DMAX
$  T    T    T
$ ----------------------------------------------
$ Type 2 : Point output
   {ymd0} {hms0}   0  {ymdf} {hmsf}
$     5.0   5.0  'Outside   '
$     0.0   0.0  'STOPSTRING'
$ --------------------------------------------
$ Type 3 : Output along track.
   {ymd0} {hms0}   0  {ymdf} {hmsf}
$ --------------------------------------------
$ Type 4 : Restart files (no additional data required).
   {ymd0} {hms0}   86400  {ymdf} {hmsf}
$ --------------------------------------------
$ Type 5 : Boundary data (no additional data required).
   {ymd0} {hms0}   0  {ymdf} {hmsf}
$ ---------------------------------------------
$ Type 6 : Separated wave field data (dummy for now).
   {ymd0} {hms0}   0  {ymdf} {hmsf}
$ ---------------------------------------------
$ Type 7 : Coupling. (must be fully commented if not used with switch COU)
   {ymd0} {hms0}   300  {ymdf} {hmsf}
  N
$ Sent fields by ww3:
  TWI AHS FP 
$ Received fields by ww3:
  IC1 WND ICE
$ --- --------------------------------------
$ Homogeneous field data :
   'STP'
$ End of input file                                                    $
$ -------------------------------------------------------------------- $
'''

  #form timestrings
  ymd0 = t0.strftime('%Y%m%d'); hms0 = t0.strftime('%H%M%S')
  ymdf = tf.strftime('%Y%m%d'); hmsf = tf.strftime('%H%M%S')
 
  return s.format(ymd0=ymd0,hms0=hms0, ymdf=ymdf,hmsf=hmsf)

def form_oasis_namcouple(RUNTIME):
  s = '''#########################################################################
 $NFIELDS
             11
 $END
###########################################################################
 $RUNTIME
  {RUNTIME}
 $END
###########################################################################
 $NLOGPRT
  1 0
 $END
###########################################################################
 $STRINGS
strau:strav:wind strau:strav:wind 1 3600  0  atmoForce.nc INPUT
#
frain:fsnow:Tair:qair:Pair:rhoa frain:fsnow:Tair:qair:Pair:rhoa 1 3600  0  atmoForce.nc INPUT
#
swdvdr:swdvdf:swdidr:swdidf:lwd swdvdr:swdvdf:swdidr:swdidf:lwd 1 3600  0 atmoForce.nc INPUT
#
sst_ext:sss_ext:hmix_ext sst_ext:sss_ext:hmix_ext 1 86400 0 oceanForce.nc INPUT
#
uocn_ext:vocn_ext:qdp_ext uocn_ext:vocn_ext:qdp_ext  1 86400  0  oceanForce.nc INPUT
#
aice_bry:vice_bry:vsno_bry aice_bry:vice_bry:vsno_bry  1 86400  0  oceanForce.nc INPUT
#
aice:hi:floeDiam WW3__ICE:WW3__IC1:WW3__IC5 1 600 2 restart_ice2wave.nc EXPOUT
701311 1 697939 1 gbar ww3t LAG=+120
R  0  R  0
LOCTRANS MAPPING
INSTANT
rmp_gbar_to_ww3t_DISTWGT.nc
#
WW3_TWIU:WW3_TWIV:WW3_LBRK:WW3_ZBRK strwaveu:strwavev:lBreak:zBreak 1 600 2 restart_wave2ice.nc EXPOUT
697939 1 701311 1 ww3t gbar LAG=+300
R  0  R  0
LOCTRANS MAPPING
INSTANT
rmp_ww3t_to_gbar_DISTWGT.nc
#
WW3__U10:WW3__V10 WW3__U10:WW3__V10 1 3600  0  ww3WindForce.nc INPUT
#
alb:tice alb:tice 1 3600 0  oasis_cice_out.nc OUTPUT
gbar gbar
#
sst:uvel:vvel sst:uvel:vvel 1 3600  0  oasis_cice_out.nc OUTPUT
gbar gbar
#
 $END
'''
  return s.format(RUNTIME=RUNTIME)

def demo_write_nml():
  #
  t0 = dt.datetime(2019,12,6); deltat = 120.0; npt = 720; restart_fsd = '.false.'
  s = form_ice_in(t0, deltat, npt, restart_fsd)
  with open('ice_in', 'w') as f:
    f.write(s)

  #
  t0 = dt.datetime(2019,12,6); tf = dt.datetime(2019,12,7)
  s = form_ww3_shel(t0, tf)
  with open('ww3_shel.inp', 'w') as f:
    f.write(s)
  #
  #t0 = dt.datetime(2019,12,6); tf = dt.datetime(2019,12,7)
  tRun = (tf-t0).total_seconds(); tRun = int(tRun)
  s = form_oasis_namcouple(tRun)
  with open('namcouple', 'w') as f:
    f.write(s)

if __name__=='__main__':
  #pass
  demo_write_nml()

