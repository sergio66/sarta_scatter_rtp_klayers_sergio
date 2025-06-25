# [sergio@strow-interact SymbolicLinks]$ ls -lt /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/
# total 3412
# lrwxrwxrwx 1 chepplew pi_strow      96 Mar 26 17:35 refl_therm.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/thermFactor_cris_hr_pbl_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     104 Mar 26 16:22 hdo.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_hdo_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     104 Mar 26 16:22 nh3.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_nh3_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     104 Mar 26 16:21 n2o.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_n2o_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     105 Mar 26 16:21 hno3.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_hno3_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     104 Mar 26 16:20 so2.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_so2_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     104 Mar 26 16:19 xnte_2x7term.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/iaa/CRIS_HR_PBL_xnte_415ppm_2x7term_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow      71 Jan 23 10:25 set2_ret.dat -> /home/chepplew/data/sarta/prod_2019/cris_hr/dec2018/dbase/Coef/set2.dat
# -rwxrwxr-x 1 chepplew pi_strow 3286848 Jan 22 11:42 set2_0p5.dat
# lrwxrwxrwx 1 chepplew pi_strow     111 Jan 22 09:54 co2.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_co2_5term_fowp_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     103 Jan 22 09:51 optran.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_optran_fmw_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     106 Jan 22 09:50 set7.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set7_fowp_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     106 Jan 22 09:49 set6.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set6_fowp_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     106 Jan 22 09:48 set5.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set5_fowp_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     106 Jan 22 09:48 set4.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set4_fcow_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     101 Jan 22 09:48 set3.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set3_fmw_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     102 Jan 22 09:47 set2.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set2_fowp_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     102 Jan 22 09:47 set1.dat -> /home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set1_fowp_g4.dat
# -rwxrwxr-x 1 chepplew pi_strow  176565 Jan 21 14:41 tunmlt_ones.txt
# -rwxrwxr-x 1 chepplew pi_strow    2063 Jan 21 14:40 fx_pbl.txt
# -rwxrwxr-x 1 chepplew pi_strow   16448 Jan 21 14:40 refprof_400ppm_pbl

########################################################################

echo "making PBL_100/CRIS_HR/2025/ links ...."
cd PBL_100/CRIS_HR/2025/

rm refl_therm.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/thermFactor_cris_hr_pbl_g4.dat           refl_therm.dat
rm hdo.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_hdo_all_g4.dat   hdo.dat
rm nh3.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_nh3_all_g4.dat   nh3.dat
rm n2o.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_n2o_all_g4.dat   n2o.dat
rm hno3.dat;         ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_hno3_all_g4.dat  hno3.dat
rm so2.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_PBL_R49_cutcoef_so2_all_g4.dat   so2.dat
rm xnte_2x7term.dat; ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/iaa/CRIS_HR_PBL_xnte_415ppm_2x7term_g4.dat   xnte_2x7term.dat
 
rm co2.dat;      ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_co2_5term_fowp_sun_g4.dat  co2.dat
rm optran.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_optran_fmw_g4.dat          optran.dat
rm set7.dat;     ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set7_fowp_sun_g4.dat       set7.dat
rm set6.dat;     ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set6_fowp_sun_g4.dat       set6.dat
rm set5.dat;     ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set5_fowp_sun_g4.dat       set5.dat
rm set4.dat;     ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set4_fcow_sun_g4.dat       set4.dat
rm set3.dat;     ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set3_fmw_g4.dat            set3.dat
rm set2.dat;     ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set2_fowp_g4.dat           set2.dat
rm set1.dat;     ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/fitc/r49/CRIS_HR_r49_cutcoef_set1_fowp_g4.dat           set1.dat
rm set2_ret.dat; ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2019/cris_hr/dec2018/dbase/Coef/set2.dat                                          set2_ret.dat   

rm set2_0p5.dat;        ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/set2_0p5.dat       .
rm tunmlt_ones.txt;     ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/tunmlt_ones.txt    .
rm fx_pbl.txt;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/fx_pbl.txt         .
rm refprof_400ppm_pbl;  ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/refprof_400ppm_pbl .

cd ../../../
echo "from /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac  test on chip using eg "
echo "     ../bin/jac_crisg4_hires_feb25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new_PBL fin=../TEST_PBL100_vs_AIRS100/cris_pbl_layering_example_clrsky.op.rtp fout=junk.rp.rtp"
