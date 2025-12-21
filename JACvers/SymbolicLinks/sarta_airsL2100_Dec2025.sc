# sergio@strow-interact JACvers]$ ls -lt /home/sergio/asl/s1/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/dbase/Coef/
# total 36
# lrwxrwxrwx 1 chepplew pi_strow   100 Jun  1  2023 xnte.dat -> /home/chepplew//data/sarta/prod_2022/airs_l1c/apr2021/fitc/r49/AIRS_L1C_R49_cutcoef_xnte_2x7term.dat
# -rwxrwxr-x 1 chepplew pi_strow 14760 Jun  1  2023 xnte_v02.dat
# lrwxrwxrwx 1 chepplew pi_strow   100 Oct 17  2022 nte_7term.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/apr2021/fitc/r49/AIRS_L1C_R49_cutcoef_nte_all_7term.dat
# lrwxrwxrwx 1 chepplew pi_strow    94 Aug 29  2022 hdo.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_R49_cutcoef_hdo_all.dat
# lrwxrwxrwx 1 chepplew pi_strow    91 Aug 27  2022 therm.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/thermFactor_airs_l1c_2834.dat
# lrwxrwxrwx 1 chepplew pi_strow   100 Aug 26  2022 nte.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_R49_cutcoef_nte_all_7term.dat
# lrwxrwxrwx 1 chepplew pi_strow    95 Aug 26  2022 hno3.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_R49_cutcoef_hno3_all.dat
# lrwxrwxrwx 1 chepplew pi_strow    94 Aug 26  2022 nh3.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_R49_cutcoef_nh3_all.dat
# lrwxrwxrwx 1 chepplew pi_strow    94 Aug 26  2022 so2.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_R49_cutcoef_so2_all.dat
# lrwxrwxrwx 1 chepplew pi_strow    94 Aug 26  2022 n2o.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_R49_cutcoef_n2o_all.dat
# lrwxrwxrwx 1 chepplew pi_strow    99 Aug 26  2022 optran.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_mergecoef_optran_fmw.dat
# lrwxrwxrwx 1 chepplew pi_strow    77 Aug 25  2022 tunmlt.txt -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/tunmlt_ones.txt
# -rwxrwxr-x 1 chepplew pi_strow 16338 Aug 25  2022 refprof_trace400
# -rwxrwxr-x 1 chepplew pi_strow  2038 Aug 25  2022 fx.txt
# lrwxrwxrwx 1 chepplew pi_strow   105 Aug 25  2022 co2.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_cutcoef_co2_5term_fowp_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow   100 Aug 25  2022 set7.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_cutcoef_set7_fowp_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow   100 Aug 25  2022 set5.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_cutcoef_set5_fowp_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow   100 Aug 25  2022 set6.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_cutcoef_set6_fowp_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow   100 Aug 25  2022 set4.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_cutcoef_set4_fcow_sun.dat
# lrwxrwxrwx 1 chepplew pi_strow    95 Aug 25  2022 set3.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_cutcoef_set3_fmw.dat
# lrwxrwxrwx 1 chepplew pi_strow    96 Aug 25  2022 set2.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_cutcoef_set2_fowp.dat
# lrwxrwxrwx 1 chepplew pi_strow    96 Aug 25  2022 set1.dat -> /home/chepplew/data/sarta/prod_2022/airs_l1c/jul2022/fitc/r49/AIRS_L1C_r49_cutcoef_set1_fowp.dat

########################################################################

echo "making L2_100/AIRS/2022 links .... copied from chepplew on yam"
cd L2_100/AIRS/2022

rm xnte.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/xnte.dat            xnte.dat
rm xnte_v02.dat;      ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/xnte_v02.dat        xnte_v02.dat
rm nte_7term.dat;     ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/nte_7term.dat       nte_7term.dat
rm therm.dat;         ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/therm.dat           therm.dat
rm hdo.dat;           ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/hdo.dat             hdo.dat
rm nte.dat;           ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/nte_7term.dat       nte.dat
rm hno3.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/hno3.dat            hno3.dat
rm nh3.dat;           ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/nh3.dat             nh3.dat
rm so2.dat;           ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/so2.dat             so2.dat
rm n2o.dat;           ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/n2o.dat             n2o.dat
rm optran.dat;        ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/optran.dat          optran.dat
rm tunmlt.txt;        ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/tunmlt_ones.txt     tunmlt.txt
rm refprof_trace400;  ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/refprof_trace400    refprof_trace400 
rm fx.txt;            ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/fx.txt              fx.txt
rm co2.dat;           ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/co2.dat             co2.dat
rm set7.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/set7.dat            set7.dat
rm set6.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/set6.dat            set6.dat
rm set5.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/set5.dat            set5.dat
rm set4.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/set4.dat            set4.dat
rm set3.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/set3.dat            set3.dat
rm set2.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/set2.dat            set2.dat
rm set1.dat;          ln -s /umbc/rs/pi_sergio/Yam_CHEPPLEW/SARTA_COEFFS/Prod2022/AIRSL1C/set1.dat            set1.dat

cd ../../../
echo "from /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac  test on chip using eg "
echo "     ../bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=../TEST_RTP/newdayx_1_100_12150.op.rtp fout=junk.rp.rtp"
