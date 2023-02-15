#V3.30.18.00;_safe;_compile_date:_Sep 30 2021;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#_data_and_control_files: GOAPcodYEARMonDay.dat, Model19.1_Year.ctl
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
6 #_Nblock_Patterns
 2 4 1 1 1 1 #_blocks_per_pattern 
# begin and end years of blocks
 1996 2005 2006 2022
 1990 2004 2005 2006 2007 2016 2017 2022
 2017 2022
 2014 2016
 1976 1976
 1976 2006
#
# controls for all timevary parameters 
1 #_time-vary parm bound check (1=warn relative to base parm bounds; 3=no bound check); Also see env (3) and dev (5) options to constrain with base bounds
#
# AUTOGEN
 1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: P(y)=f(TVP,env_Zscore) w/ logit to stay in min-max;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  5: like 4 with logit transform to stay in base min-max
#_DevLinks(more):  21-25 keep last dev for rest of years
#
#_Prior_codes:  0=none; 6=normal; 1=symmetric beta; 2=CASAL's beta; 3=lognormal; 4=lognormal with biascorr; 5=gamma
#
# setup for M, growth, wt-len, maturity, fecundity, (hermaphro), recr_distr, cohort_grow, (movement), (age error), (catch_mult), sex ratio 
#_NATMORT
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=BETA:_Maunder_link_to_maturity
  #_no additional input for selected M option; read 1P per morph
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.5 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
2 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
1 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
 0.1 1.5 0.499189 -0.81 0.41 3 5 0 0 0 0 0 4 2 # NatM_uniform_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 0 50 6.3923 6.1252 99 0 1 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 70 130 99.4617 99.46 0.015 6 1 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0 1 0.188842 0.1966 0.03 6 1 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0 10 3.82037 0 0 0 10 0 0 0 0 0 0 0 # SD_young_Fem_GP_1
 0 20 7.42895 0 0 0 10 0 0 0 0 0 0 0 # SD_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -99 99 5.63096e-06 0 0 0 -3 0 0 0 0 0 0 0 # Wtlen_1_Fem_GP_1
 -99 99 3.1306 0 0 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 -99 99 53.7 0 0 0 -1 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -99 99 -0.273657 0 0 0 -1 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 -99 99 1 0 0 0 -1 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem_GP_1
 -99 99 0 0 0 0 -1 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem_GP_1
# Hermaphroditism
#  Recruitment Distribution  
#  Cohort growth dev base
 0.1 10 1 1 1 0 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
 -10 10 3 0 0 0 -5 0 0 0 0 0 0 0 # AgeKeyParm1
 -10 10 0 0 0 0 -10 0 0 0 0 0 6 2 # AgeKeyParm2
 -10 10 0 0 0 0 -10 0 0 0 0 0 6 2 # AgeKeyParm3
 -10 10 0 0 0 0 -1 0 0 0 0 0 0 0 # AgeKeyParm4
 -10 10 0.57 0 0 0 -1 0 0 0 0 0 0 0 # AgeKeyParm5
 -10 10 1.16 0 0 0 -1 0 0 0 0 0 0 0 # AgeKeyParm6
 -10 10 0 0 0 0 -1 0 0 0 0 0 0 0 # AgeKeyParm7
#  catch multiplier
#  fraction female, by GP
 1e-06 0.999999 0.5 0.5 0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#  M2 parameter for each predator fleet
#
# timevary MG parameters 
#_ LO HI INIT PRIOR PR_SD PR_type  PHASE
 0.1 2 0.874837 -0.81 0.41 3 6 # NatM_uniform_Fem_GP_1_BLK4repl_2014
 -9 9 0.391521 0 0 0 9 # AgeKeyParm2_BLK6repl_1976
 -9 9 -0.184608 0 0 0 9 # AgeKeyParm3_BLK6repl_1976
# info on dev vectors created for MGparms are reported with other devs after tag parameter section 
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
            10            20       13.1552             0             0             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
             0             1             1             1             0             0         -1          0          0          0          0          0          0          0 # SR_BH_steep
             0            10          0.44          0.44             0             0         -4          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             0             0         -3          0          0          0          0          0          5          1 # SR_regime
           -99            99             0             0             0             0         -1          0          0          0          0          0          0          0 # SR_autocorr
# timevary SR parameters
 -10 10 -0.523818 0 0 0 1 # SR_regime_BLK5add_1976
2 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1978 # first year of main recr_devs; early devs can preceed this era
2019 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase 
1 # (0/1) to read 13 advanced options
 1967 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 2 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1964.4 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1979.7 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2012.9 #_last_yr_fullbias_adj_in_MPD
 2017 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 0.9089 #_max_bias_adj_in_MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1967E 1968E 1969E 1970E 1971E 1972E 1973E 1974E 1975E 1976E 1977E 1978R 1979R 1980R 1981R 1982R 1983R 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017R 2018R 2019F 2020F 2021F 2022F 2023F 2024F 2025F 2026F 2027F 2028F 2029F 2030F 2031F 2032F 2033F 2034F 2035F 2036F
#  -0.449943 -0.302048 -0.310575 -0.240026 -0.116886 0.0120751 0.57627 0.159024 -0.297455 -0.214135 0.95527 0.115909 -0.0505858 0.292027 0.641131 0.622367 0.433308 0.302433 0.940895 0.137333 0.421361 0.352351 0.444418 0.60813 0.162637 0.0552266 -0.325384 -0.214222 0.0989621 -0.316777 -0.280733 -0.488254 -0.134835 0.042304 -0.364693 -0.74151 -0.543782 -0.359824 0.037633 0.464824 0.138027 0.457153 0.0812119 0.212659 0.422766 1.10052 0.641069 -0.416897 -0.581293 -0.527376 -0.966231 -0.660234 -1.28617 -0.739425 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#
#Fishing Mortality info 
0.3489 # F ballpark value in units of annual_F
-1999 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope midseason rate; 2=F as parameter; 3=F as hybrid; 4=fleet-specific parm/hybrid (#4 is superset of #2 and #3 and is recommended)
5 # max F (methods 2-4) or harvest fraction (method 1)
5  # N iterations for tuning in hybrid mode; recommend 3 (faster) to 5 (more precise if many fleets)
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 0
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#
# F rates by fleet x season
# Yr:  1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030 2031 2032 2033 2034 2035 2036
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# FshTrawl 0.00219026 0.0114252 0.0141297 0.0379011 0.0248807 0.0154197 0.0236247 0.0186424 0.0113289 0.0109809 0.0384011 0.0465244 0.0635516 0.153977 0.164735 0.171011 0.123814 0.102032 0.139555 0.172507 0.215776 0.222007 0.231251 0.178708 0.185513 0.158587 0.156221 0.147638 0.302981 0.294005 0.156897 0.225035 0.133008 0.170864 0.113079 0.12726 0.128519 0.179626 0.204604 0.185935 0.17408 0.0550868 0.0649757 0.045178 0.091178 0.116541 0.107346 0.109333 0.142798 0.163781 0.163781 0.163781 0.163781 0.163781 0.163781 0.163781 0.163781 0.163781 0.163781 0.163781
# FshLL 0.00741466 0.0350328 0.0474106 0.110176 0.0688238 0.0580151 0.07089 0.0433729 0.0457072 0.0715806 0.0207891 0.0104397 0.0102652 0.0184511 0.0204363 0.046197 0.0273683 0.0205711 0.0338618 0.0355033 0.0457241 0.0503266 0.0717986 0.0756343 0.0705901 0.110472 0.0999284 0.104373 0.0606752 0.0882747 0.131513 0.159217 0.169568 0.168432 0.15552 0.127307 0.0870458 0.128311 0.157896 0.131624 0.133425 0.0565731 0.0571508 0.0219493 0.0850431 0.0889106 0.0818956 0.083411 0.108942 0.12495 0.12495 0.12495 0.12495 0.12495 0.12495 0.12495 0.12495 0.12495 0.12495 0.12495
# FshPot 0 0 0 0 0 0 0 0 0 0 0.00186913 0.00409128 0.0010465 0.0161393 0.0330121 0.035571 0.035375 0.0322533 0.0574887 0.0489817 0.0449119 0.0636597 0.131318 0.131254 0.0584529 0.0669069 0.186681 0.22656 0.219927 0.241235 0.302717 0.326563 0.247784 0.300897 0.344656 0.254307 0.199384 0.290344 0.441631 0.512724 0.481601 0.14146 0.139594 0.0316868 0.159522 0.290006 0.267125 0.272068 0.355343 0.407558 0.407558 0.407558 0.407558 0.407558 0.407558 0.407558 0.407558 0.407558 0.407558 0.407558
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         4         1         0         0         0         0  #  Srv
         5         1         0         0         0         0  #  LLSrv
         6         0         0         0         0         0  #  IPHCLL
         7         0         0         0         0         0  #  ADFG
         8         3         0         0         0         0  #  SPAWN
         9         3         0         0         0         0  #  Seine
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -10            10      0.101348             0            99             0          1          0          0          0          0          0          0          0  #  LnQ_base_Srv(4)
           -10            10      0.115124             0            99             0          1        101          0          0          0          0          0          0  #  LnQ_base_LLSrv(5)
           -25            25             0             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_IPHCLL(6)
           -25            25             0             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_ADFG(7)
           -25            25      -5.39045             0            99             0         -1          0          0          0          0          0          0          0  #  LnQ_base_SPAWN(8)
            -1             5       2.24627             0            99             0         -1          0          0          0          0          0          0          0  #  Q_power_SPAWN(8)
           -50            25      -37.3057             0            99             0         -1          0          0          0          0          0          0          0  #  LnQ_base_Seine(9)
            -1             5       2.10009             0            99             0         -1          0          0          0          0          0          0          0  #  Q_power_Seine(9)
# timevary Q parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type     PHASE  #  parm_name
           -10            10       1.12292             0            99             0      5  # LnQ_base_LLSrv(5)_ENV_mult
# info on dev vectors created for Q parms are reported with other devs after tag parameter section 
#
#_size_selex_patterns
#Pattern:_0;  parm=0; selex=1.0 for all sizes
#Pattern:_1;  parm=2; logistic; with 95% width specification
#Pattern:_2;  parm=6; modification of pattern 24 with improved sex-specific offset
#Pattern:_5;  parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_11; parm=2; selex=1.0  for specified min-max population length bin range
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6;  parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8;  parm=8; double_logistic with smooth transitions and constant above Linf option
#Pattern:_9;  parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in length
#Pattern:_27; parm=special+3; cubic spline in length; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=special+3+2; cubic spline; like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 0 0 0 # 1 FshTrawl
 24 0 0 0 # 2 FshLL
 24 0 0 0 # 3 FshPot
 24 0 0 0 # 4 Srv
 24 0 0 0 # 5 LLSrv
 15 0 0 4 # 6 IPHCLL
 15 0 0 4 # 7 ADFG
 0 0 0 0 # 8 SPAWN
 0 0 0 0 # 9 Seine
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
 10 0 0 0 # 1 FshTrawl
 10 0 0 0 # 2 FshLL
 10 0 0 0 # 3 FshPot
 10 0 0 0 # 4 Srv
 10 0 0 0 # 5 LLSrv
 10 0 0 0 # 6 IPHCLL
 10 0 0 0 # 7 ADFG
 0 0 0 0 # 8 SPAWN
 0 0 0 0 # 9 Seine
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   FshTrawl LenSelex
            10           110       57.6551            50             0             0          1          0          1       1977       1989          3          2          2  #  Size_DblN_peak_FshTrawl(1)
           -20            10      -4.57915             0             0             0          2          0          0       1977       1989          0          2          2  #  Size_DblN_top_logit_FshTrawl(1)
           -10            10       5.09148             0             0             0          2          0          1       1977       1989          3          2          2  #  Size_DblN_ascend_se_FshTrawl(1)
           -10            10     -0.301675            10             0             0          2          0          1       1977       1989          3          2          2  #  Size_DblN_descend_se_FshTrawl(1)
         -1000       2.71828          -999           -10             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_FshTrawl(1)
           -10            10            10            10             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_FshTrawl(1)
# 2   FshLL LenSelex
            10            90       66.4253            50             0             0          1          0          1       1978       1989          3          2          2  #  Size_DblN_peak_FshLL(2)
           -20            10      -5.06337             0             0             0          2          0          0       1978       1989          0          2          2  #  Size_DblN_top_logit_FshLL(2)
           -10            10       5.14687             0             0             0          2          0          1       1978       1989          3          2          2  #  Size_DblN_ascend_se_FshLL(2)
             0            10            10            10             0             0         -2          0          0          0          0          0          2          2  #  Size_DblN_descend_se_FshLL(2)
         -1000       2.71828          -999           -10             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_FshLL(2)
           -10            10            10            10             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_FshLL(2)
# 3   FshPot LenSelex
            10            90       70.7282            50             0             0          1          0          0          0          0          0          3          2  #  Size_DblN_peak_FshPot(3)
           -20            10      -12.0643             0             0             0          2          0          0          0          0          0          3          2  #  Size_DblN_top_logit_FshPot(3)
           -10            10       5.01625             0             0             0          2          0          0          0          0          0          3          2  #  Size_DblN_ascend_se_FshPot(3)
             0            10        4.0558            10             0             0          2          0          0          0          0          0          0          0  #  Size_DblN_descend_se_FshPot(3)
         -1000       2.71828          -999           -10             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_FshPot(3)
           -10            10      0.284913            10             0             0          2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_FshPot(3)
# 4   Srv LenSelex
            10            90       60.2246            50             0             0          1          0          0          0          0          0          1          2  #  Size_DblN_peak_Srv(4)
           -20            10      -11.9482             0             0             0          2          0          0          0          0          0          1          2  #  Size_DblN_top_logit_Srv(4)
           -10            10       5.54831             0             0             0          2          0          0          0          0          0          1          2  #  Size_DblN_ascend_se_Srv(4)
             0            10        3.9945            10             0             0          5          0          0          0          0          0          1          2  #  Size_DblN_descend_se_Srv(4)
           -10       2.71828      -3.10266           -10             0             0          2          0          0          0          0          0          1          2  #  Size_DblN_start_logit_Srv(4)
           -10            10      -0.66555            10             0             0          5          0          0          0          0          0          1          2  #  Size_DblN_end_logit_Srv(4)
# 5   LLSrv LenSelex
            10            90       65.5449            50             0             0          1          0          0          0          0          0          0          0  #  Size_DblN_peak_LLSrv(5)
           -20            10       -12.478             0             0             0          2          0          0          0          0          0          0          0  #  Size_DblN_top_logit_LLSrv(5)
           -10            10       4.68122             0             0             0          2          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_LLSrv(5)
             0            10       4.64978            10             0             0          2          0          0          0          0          0          0          0  #  Size_DblN_descend_se_LLSrv(5)
         -1000       2.71828          -999           -10             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_LLSrv(5)
           -10            10     -0.466947            10             0             0          2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_LLSrv(5)
# 6   IPHCLL LenSelex
# 7   ADFG LenSelex
# 8   SPAWN LenSelex
# 9   Seine LenSelex
# 1   FshTrawl AgeSelex
# 2   FshLL AgeSelex
# 3   FshPot AgeSelex
# 4   Srv AgeSelex
# 5   LLSrv AgeSelex
# 6   IPHCLL AgeSelex
# 7   ADFG AgeSelex
# 8   SPAWN AgeSelex
# 9   Seine AgeSelex
#_No_Dirichlet parameters
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
            10            90       78.6311            50             0             0      1  # Size_DblN_peak_FshTrawl(1)_BLK2repl_1990
            10           100       94.5791            50             0             0      1  # Size_DblN_peak_FshTrawl(1)_BLK2repl_2005
            10            90       79.9048            50             0             0      1  # Size_DblN_peak_FshTrawl(1)_BLK2repl_2007
            10            90       77.2556            50             0             0      -1  # Size_DblN_peak_FshTrawl(1)_BLK2repl_2017
        0.0001             2           0.2           0.2           0.5             6      -5  # Size_DblN_peak_FshTrawl(1)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Size_DblN_peak_FshTrawl(1)_dev_autocorr
           -20            10      0.434816             0             0             0      2  # Size_DblN_top_logit_FshTrawl(1)_BLK2repl_1990
           -20            10      -4.93158             0             0             0      2  # Size_DblN_top_logit_FshTrawl(1)_BLK2repl_2005
           -20            10      -4.75885             0             0             0      2  # Size_DblN_top_logit_FshTrawl(1)_BLK2repl_2007
           -20            10      -4.88944             0             0             0      2  # Size_DblN_top_logit_FshTrawl(1)_BLK2repl_2017
           -10            10       6.02471             0             0             0      2  # Size_DblN_ascend_se_FshTrawl(1)_BLK2repl_1990
           -10            10       6.26352             0             0             0      2  # Size_DblN_ascend_se_FshTrawl(1)_BLK2repl_2005
           -10            10       6.24936             0             0             0      2  # Size_DblN_ascend_se_FshTrawl(1)_BLK2repl_2007
           -10            10       6.29422             0             0             0      2  # Size_DblN_ascend_se_FshTrawl(1)_BLK2repl_2017
        0.0001             2           0.2           0.2           0.5             6      -5  # Size_DblN_ascend_se_FshTrawl(1)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Size_DblN_ascend_se_FshTrawl(1)_dev_autocorr
             0            10       5.10604            10             0             0      2  # Size_DblN_descend_se_FshTrawl(1)_BLK2repl_1990
             0            10       5.15973            10             0             0      2  # Size_DblN_descend_se_FshTrawl(1)_BLK2repl_2005
             0            10       5.79906            10             0             0      2  # Size_DblN_descend_se_FshTrawl(1)_BLK2repl_2007
             0            10       5.26328            10             0             0      2  # Size_DblN_descend_se_FshTrawl(1)_BLK2repl_2017
        0.0001             2           0.2           0.2           0.5             6      -5  # Size_DblN_descend_se_FshTrawl(1)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Size_DblN_descend_se_FshTrawl(1)_dev_autocorr
            10            90       71.8699            50             0             0      1  # Size_DblN_peak_FshLL(2)_BLK2repl_1990
            10            90       68.4629            50             0             0      1  # Size_DblN_peak_FshLL(2)_BLK2repl_2005
            10            90       72.5472            50             0             0      1  # Size_DblN_peak_FshLL(2)_BLK2repl_2007
            10            90       71.9494            50             0             0      1  # Size_DblN_peak_FshLL(2)_BLK2repl_2017
        0.0001             2           0.2           0.2           0.5             6      -5  # Size_DblN_peak_FshLL(2)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Size_DblN_peak_FshLL(2)_dev_autocorr
           -20            10     -0.224231             0             0             0      2  # Size_DblN_top_logit_FshLL(2)_BLK2repl_1990
           -20            10      -5.10633             0             0             0      2  # Size_DblN_top_logit_FshLL(2)_BLK2repl_2005
           -20            10       -5.1146             0             0             0      2  # Size_DblN_top_logit_FshLL(2)_BLK2repl_2007
           -20            10      -4.69045             0             0             0      2  # Size_DblN_top_logit_FshLL(2)_BLK2repl_2017
           -10            10       5.36118             0             0             0      2  # Size_DblN_ascend_se_FshLL(2)_BLK2repl_1990
           -10            10       5.25191             0             0             0      2  # Size_DblN_ascend_se_FshLL(2)_BLK2repl_2005
           -10            10       5.44255             0             0             0      2  # Size_DblN_ascend_se_FshLL(2)_BLK2repl_2007
           -10            10       5.35388             0             0             0      2  # Size_DblN_ascend_se_FshLL(2)_BLK2repl_2017
        0.0001             2           0.2           0.2           0.5             6      -5  # Size_DblN_ascend_se_FshLL(2)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Size_DblN_ascend_se_FshLL(2)_dev_autocorr
             0            10            10            10             0             0      -2  # Size_DblN_descend_se_FshLL(2)_BLK2repl_1990
             0            10            10            10             0             0      -2  # Size_DblN_descend_se_FshLL(2)_BLK2repl_2005
             0            10            10            10             0             0      -2  # Size_DblN_descend_se_FshLL(2)_BLK2repl_2007
             0            10            10            10             0             0      -2  # Size_DblN_descend_se_FshLL(2)_BLK2repl_2017
            10            90       75.9646            50             0             0      1  # Size_DblN_peak_FshPot(3)_BLK3repl_2017
           -20            10       2.80531             0             0             0      2  # Size_DblN_top_logit_FshPot(3)_BLK3repl_2017
           -10            10       5.47946             0             0             0      2  # Size_DblN_ascend_se_FshPot(3)_BLK3repl_2017
            10            90        59.965            50             0             0      1  # Size_DblN_peak_Srv(4)_BLK1repl_1996
            10            90       56.5473            50             0             0      1  # Size_DblN_peak_Srv(4)_BLK1repl_2006
           -20            10      -4.15407             0             0             0      2  # Size_DblN_top_logit_Srv(4)_BLK1repl_1996
           -20            10      -3.73717             0             0             0      2  # Size_DblN_top_logit_Srv(4)_BLK1repl_2006
           -10            10        5.4719             0             0             0      2  # Size_DblN_ascend_se_Srv(4)_BLK1repl_1996
           -10            10       5.09389             0             0             0      2  # Size_DblN_ascend_se_Srv(4)_BLK1repl_2006
             0            10       4.21246            10             0             0      5  # Size_DblN_descend_se_Srv(4)_BLK1repl_1996
             0            10         4.195            10             0             0      5  # Size_DblN_descend_se_Srv(4)_BLK1repl_2006
           -10       2.71828      -2.09021           -10             0             0      2  # Size_DblN_start_logit_Srv(4)_BLK1repl_1996
           -10       2.71828      -2.01521           -10             0             0      2  # Size_DblN_start_logit_Srv(4)_BLK1repl_2006
             0            10            10            10             0             0      -5  # Size_DblN_end_logit_Srv(4)_BLK1repl_1996
             0            10            10            10             0             0      -5  # Size_DblN_end_logit_Srv(4)_BLK1repl_2006
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1)
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      1     1     1     4     2     0     0     0     0     0     0     0
#      1    15     2     6     2     0     0     0     0     0     0     0
#      1    16     3     6     2     0     0     0     0     0     0     0
#      2     4     4     5     1     0     0     0     0     0     0     0
#      3     2     5     0     0     1     1     0     0     0     0     0
#      5     1     6     2     2     0     0     1     1  1977  1989     3 -0.486887 -0.242494 0.0138779 0.590577 -0.40375 -0.238779 -0.00150576 -0.392678 0.736158 -0.22571 -0.054645 0.138429 0.567387
#      5     2    12     2     2     0     0     0     0     0     0     0
#      5     3    16     2     2     0     0     2     1  1977  1989     3 -0.319235 -0.00895258 -0.0142082 0.241316 -0.42071 0.0792314 -0.306249 -0.818384 0.462508 0.515134 -0.073697 -0.0731665 0.736553
#      5     4    22     2     2     0     0     3     1  1977  1989     3 3.03249e-07 1.6031e-07 3.00129e-08 3.55428e-07 -1.95775e-07 1.13557e-07 1.10029e-07 -3.73772e-07 2.63313e-07 1.58038e-07 3.83704e-08 -1.38391e-07 9.95e-08
#      5     7    28     2     2     0     0     4     1  1978  1989     3 -0.394629 -0.27587 0.312517 -0.762174 -0.719701 -0.665902 0.00819144 1.4492 1.27516 -0.0149766 -0.182897 -0.0626867
#      5     8    34     2     2     0     0     0     0     0     0     0
#      5     9    38     2     2     0     0     5     1  1978  1989     3 -0.176024 -0.138734 0.549014 -0.74311 -0.488345 -0.742542 0.0351698 0.99253 0.87756 0.0103257 -0.184667 0.0227802
#      5    10    44     2     2     0     0     0     0     0     0     0
#      5    13    48     3     2     0     0     0     0     0     0     0
#      5    14    49     3     2     0     0     0     0     0     0     0
#      5    15    50     3     2     0     0     0     0     0     0     0
#      5    19    51     1     2     0     0     0     0     0     0     0
#      5    20    53     1     2     0     0     0     0     0     0     0
#      5    21    55     1     2     0     0     0     0     0     0     0
#      5    22    57     1     2     0     0     0     0     0     0     0
#      5    23    59     1     2     0     0     0     0     0     0     0
#      5    24    61     1     2     0     0     0     0     0     0     0
     #
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
      1      8         0
      1      9         0
      1      6         0
      1      7         0
      1      4         0
      1      5         0
      4      1         1
      4      2         1
      4      3         1
      4      4         1
      4      5         1
      6      3         1
      6      4         1
      6      5         0
      4      6         0
      5      6         0
      6      6         0
      4      7         0
      5      7         0
      6      7         0
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 4 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 1 8 1 0 1
 1 6 1 0 1
 1 7 1 0 1
 1 9 1 0 1
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 #_CPUE/survey:_1
#  0 #_CPUE/survey:_2
#  0 #_CPUE/survey:_3
#  1 #_CPUE/survey:_4
#  1 #_CPUE/survey:_5
#  0 #_CPUE/survey:_6
#  0 #_CPUE/survey:_7
#  0 #_CPUE/survey:_8
#  0 #_CPUE/survey:_9
#  1 #_lencomp:_1
#  1 #_lencomp:_2
#  1 #_lencomp:_3
#  1 #_lencomp:_4
#  1 #_lencomp:_5
#  0 #_lencomp:_6
#  0 #_lencomp:_7
#  0 #_lencomp:_8
#  0 #_lencomp:_9
#  1 #_agecomp:_1
#  1 #_agecomp:_2
#  1 #_agecomp:_3
#  1 #_agecomp:_4
#  0 #_agecomp:_5
#  0 #_agecomp:_6
#  0 #_agecomp:_7
#  0 #_agecomp:_8
#  0 #_agecomp:_9
#  0 #_size-age:_1
#  0 #_size-age:_2
#  0 #_size-age:_3
#  1 #_size-age:_4
#  0 #_size-age:_5
#  0 #_size-age:_6
#  0 #_size-age:_7
#  0 #_size-age:_8
#  0 #_size-age:_9
#  1 #_init_equ_catch1
#  1 #_init_equ_catch2
#  1 #_init_equ_catch3
#  1 #_init_equ_catch4
#  1 #_init_equ_catch5
#  1 #_init_equ_catch6
#  1 #_init_equ_catch7
#  1 #_init_equ_catch8
#  1 #_init_equ_catch9
#  1 #_recruitments
#  1 #_parameter-priors
#  1 #_parameter-dev-vectors
#  1 #_crashPenLambda
#  0 # F_ballpark_lambda
0 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 # 0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
 # 0 0 # Growth: (1) growth pattern, (2) growth ages
 # 0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
 # -1 # list of bin #'s for selex std (-1 in first bin to self-generate)
 # -1 # list of ages for growth std (-1 in first bin to self-generate)
 # -1 # list of ages for NatAge std (-1 in first bin to self-generate)
999

