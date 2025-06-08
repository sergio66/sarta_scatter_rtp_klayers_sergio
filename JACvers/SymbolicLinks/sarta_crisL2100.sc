# sergio@strow-interact JACvers]$ ls -lt /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/dbase/Coef/
# total 196
# lrwxrwxrwx 1 chepplew pi_strow     91 Nov 14  2023 xnte_14term.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_xnte_v1_2x7term_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow    108 Nov 19  2022 set5.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/check2_CRIS_HR_r49_cutcoef_set5_fowp_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow    108 Nov 18  2022 set6.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/check2_CRIS_HR_r49_cutcoef_set6_fowp_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     98 Oct 12  2022 nte_7term.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_nte_all_7term.dat
# lrwxrwxrwx 1 chepplew pi_strow     94 Sep 14  2022 hdo.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_hdo_sw_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     95 Sep 14  2022 nh3.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_nh3_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     95 Sep 14  2022 n2o.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_n2o_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     96 Sep 14  2022 hno3.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_hno3_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     95 Sep 14  2022 so2.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_so2_all_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     87 Sep 13  2022 therm.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/thermFactor_cris_hr_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow    106 Sep 13  2022 co2.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_co2_5term_fowp_sun_g4.dat
# -rwxrwxr-x 1 chepplew pi_strow   2038 Sep 13  2022 fx.txt
# -rwxrwxr-x 1 chepplew pi_strow  16339 Sep 13  2022 refprof_400_tra
# -rwxrwxr-x 1 chepplew pi_strow 176565 Sep 13  2022 tunmlt_ones.txt
# lrwxrwxrwx 1 chepplew pi_strow    100 Sep 13  2022 optran.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_mergecoef_optran_fmw_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow    101 Sep 13  2022 set7.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set7_fowp_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow    101 Sep 13  2022 set4.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set4_fcow_sun_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     96 Sep 13  2022 set3.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set3_fmw_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     97 Sep 13  2022 set2.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set2_fowp_g4.dat
# lrwxrwxrwx 1 chepplew pi_strow     97 Sep 13  2022 set1.dat -> /home/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set1_fowp_g4.dat

########################################################################

echo "making L2_100/CRIS_HR/2022/ links ...."
cd L2_100/CRIS_HR/2022/

rm xnte_14term.dat;  ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_xnte_v1_2x7term_g4.dat                 xnte_14term.dat
rm nte_7term.dat;    ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_nte_all_7term.dat          nte_7term.dat
rm hdo.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_hdo_sw_g4.dat              hdo.dat
rm nh3.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_nh3_all_g4.dat             nh3.dat
rm n2o.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_n2o_all_g4.dat             n2o.dat
rm hno3.dat;         ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_hno3_all_g4.dat            hno3.dat
rm so2.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_R49_cutcoef_so2_all_g4.dat             so2.dat
rm therm.dat;        ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/thermFactor_cris_hr_g4.dat                     therm.dat 
rm co2.dat;          ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_co2_5term_fowp_sun_g4.dat  co2.dat

rm fx.txt;           ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/dbase/Coef/fx.txt          .
rm refprof_400_tra;  ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/dbase/Coef/refprof_400_tra .
rm tunmlt_ones.txt;  ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/dbase/Coef/tunmlt_ones.txt .


rm optran.dat; ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_mergecoef_optran_fmw_g4.dat            optran.dat
rm set7.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set7_fowp_sun_g4.dat           set7.dat
rm set6.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/check2_CRIS_HR_r49_cutcoef_set6_fowp_sun_g4.dat    set6.dat
rm set5.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/check2_CRIS_HR_r49_cutcoef_set5_fowp_sun_g4.dat    set5.dat
rm set4.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set4_fcow_sun_g4.dat           set4.dat
rm set3.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set3_fmw_g4.dat                set3.dat
rm set2.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set2_fowp_g4.dat               set2.dat
rm set1.dat;   ln -s /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/cris_hr/jul2022/fitc/r49/CRIS_HR_r49_cutcoef_set1_fowp_g4.dat               set1.dat

cd ../../../
echo "from /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac  test on chip using eg "
echo "     ../bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_newFeb2025 fin=../TEST_RTP/cris_1_100_12150.op.rtp fout=junk.rp.rtp"
