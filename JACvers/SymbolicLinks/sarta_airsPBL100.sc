# [sergio@strow-interact SymbolicLinks]$ ls -lt /home/sergio/asl/s1/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/dbase/Coef
# total 244
# lrwxrwxrwx 1 chepplew pi_strow     98 Mar 26 08:11 xnte_7term.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/iaa/AIRS_PBL_xnte_415ppm_le90_7term.dat
# lrwxrwxrwx 1 chepplew pi_strow     95 Mar 25 11:39 hdo.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_hdo_all.dat
# lrwxrwxrwx 1 chepplew pi_strow     95 Mar 25 11:38 nh3.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_nh3_all.dat
# lrwxrwxrwx 1 chepplew pi_strow     96 Mar 25 11:38 hno3.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_hno3_all.dat
# lrwxrwxrwx 1 chepplew pi_strow     95 Mar 25 11:37 so2.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_so2_all.dat
# lrwxrwxrwx 1 chepplew pi_strow     95 Mar 25 11:34 n2o.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_n2o_all.dat
# lrwxrwxrwx 1 chepplew pi_strow     95 Mar 25 11:26 xnte_2x7term.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/iaa/AIRS_PBL_xnte_415ppm_2x7term.dat
# lrwxrwxrwx 1 chepplew pi_strow     94 Jan 16 10:39 refltherm.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/refl_therm/thermFactor_airs_pbl_2834.dat
# -rwxrwxr-x 1 chepplew pi_strow   2063 Jan 16 10:29 fx_pbl.txt
# -rwxrwxr-x 1 chepplew pi_strow  16448 Jan 15 12:15 refprof_400ppm_pbl
# -rwxrwxr-x 1 chepplew pi_strow 224310 Jan 15 09:23 tunmlt_ones.txt
# lrwxrwxrwx 1 chepplew pi_strow     98 Jan 15 09:10 optran.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_optran_fmw.dat
# lrwxrwxrwx 1 chepplew pi_strow    106 Jan 15 09:10 co2.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_co2_5term_fowp_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow    101 Jan 15 09:10 set7.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set7_fowp_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow    101 Jan 15 09:10 set6.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set6_fowp_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow    101 Jan 15 09:09 set5.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set5_fowp_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow    101 Jan 15 09:09 set4.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set4_fcow_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow     96 Jan 15 09:09 set3.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set3_fmw.dat
# lrwxrwxrwx 1 chepplew pi_strow     97 Jan 15 09:09 set2.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set2_fowp.dat
# lrwxrwxrwx 1 chepplew pi_strow     97 Jan 15 09:08 set1.dat -> /home/chepplew/data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set1_fowp.dat

########################################################################

echo "making PBL_100/AIRS/2025/ links ...."
cd PBL_100/AIRS/2025/

rm xnte_7term.dat;      ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/iaa/AIRS_PBL_xnte_415ppm_le90_7term.dat    xnte_7term.dat
rm hdo.dat;             ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_hdo_all.dat       hdo.dat
rm nh3.dat;             ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_nh3_all.dat       nh3.dat
rm hno3.dat;            ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_hno3_all.dat      hno3.dat
rm so2.dat;             ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_so2_all.dat       so2.dat
rm n2o.dat;             ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/AIRS_PBL_R49_cutcoef_n2o_all.dat       n2o.dat
rm xnte_2x7term.dat;    ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/iaa/AIRS_PBL_xnte_415ppm_2x7term.dat       xnte_2x7term.dat
rm refltherm.dat;       ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/refl_therm/thermFactor_airs_pbl_2834.dat        refltherm.dat

rm fx_pbl.txt;          ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/dbase/Coef/fx_pbl.txt                                            fx_pbl.txt
rm refprof_400ppm_pbl;  ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/dbase/Coef/refprof_400ppm_pbl                                    refprof_400ppm_pbl
rm tunmlt_ones.txt;     ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/dbase/Coef/tunmlt_ones.txt                                       tunmlt_ones.txt

rm optran.dat;          ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_optran_fmw.dat          optran.dat
rm co2.dat;             ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_co2_5term_fowp_sun.dat  co2.dat
rm set7.dat;            ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set7_fowp_sun.dat       set7.dat
rm set6.dat;            ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set6_fowp_sun.dat       set6.dat  
rm set5.dat;            ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set5_fowp_sun.dat       set5.dat
rm set4.dat;            ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set4_fcow_sun.dat       set4.dat
rm set3.dat;            ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set3_fmw.dat            set3.dat
rm set2.dat;            ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set2_fowp.dat           set2.dat
rm set1.dat;            ln -s /home/sergio/asl/s1/chepplew//data/sarta/prod_2025/airs_pbl/jan2025a/fitc/r49/airs_pbl_r49_cutcoef_set1_fowp.dat           set1.dat
# 

cd ../../../
echo "from /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac  test on chip using eg "
echo "     ../bin/jac_airs_l1c_2834_cloudy_feb25_H2020_PBL fin=../TEST_PBL100_vs_AIRS100/cris_pbl_layering_example_clrsky.op.rtp fout=junk.rp.rtp"
