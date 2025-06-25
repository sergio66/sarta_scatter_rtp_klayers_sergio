# [sergio@strow-interact SymbolicLinks]$ ls /asl/rta/sarta_coef/Data_CrIS_oct16/Coef/
# ls: cannot access /asl/rta/sarta_coef/Data_CrIS_oct16/Coef/: No such file or directory
# [sergio@strow-interact SymbolicLinks]$ ls -lt /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/
# total 42804
# -rw-rwxr-- 1 chepplew pi_strow    16339 Dec  3  2019 refprof_400_tra
# -rw-rwxr-- 1 strow    pi_strow    15225 Jul 22  2018 profref_trace400
# -rw-rwxr-- 1 strow    pi_strow    16480 Jul 20  2018 profref_trace400_nh3
# -rw-rwxr-- 1 strow    pi_strow   909808 Jul 19  2018 nh3.dat
# -rw-rwxr-- 1 strow    pi_strow   176565 Dec 14  2016 tunmlt.txt
# -rw-rwxr-- 1 strow    pi_strow    89400 Dec 14  2016 therm.dat
# -rw-rwxr-- 1 strow    pi_strow   476720 Dec 14  2016 so2.dat
# -rw-rwxr-- 1 strow    pi_strow   290400 Dec 14  2016 set7.dat
# -rw-rwxr-- 1 strow    pi_strow  3824640 Dec 14  2016 set6.dat
# -rw-rwxr-- 1 strow    pi_strow   661200 Dec 14  2016 set5.dat
# -rw-rwxr-- 1 strow    pi_strow  2270016 Dec 14  2016 set4.dat
# -rw-rwxr-- 1 strow    pi_strow 12235968 Dec 14  2016 set3.dat
# -rw-rwxr-- 1 strow    pi_strow  3286848 Dec 14  2016 set2.dat
# -rw-rwxr-- 1 strow    pi_strow  6121088 Dec 14  2016 set1.dat
# -rw-rwxr-- 1 strow    pi_strow  9448408 Dec 14  2016 optran.dat
# -rw-rwxr-- 1 strow    pi_strow    11616 Dec 14  2016 nte_7term.dat
# -rw-rwxr-- 1 strow    pi_strow  1219328 Dec 14  2016 n2o.dat
# -rw-rwxr-- 1 strow    pi_strow   921120 Dec 14  2016 hno3.dat
# -rw-rwxr-- 1 strow    pi_strow     2038 Dec 14  2016 fx.txt
# -rw-rwxr-- 1 strow    pi_strow  1800288 Dec 14  2016 co2.dat

########################################################################

echo "making L2_100/IASI/2022 links ...."
cd L2_100/IASI/2022

rm refprof_400_tra;      ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/refprof_400_tra      .
rm profref_trace400;     ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/profref_trace400     .
rm profref_trace400_nh3; ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/profref_trace400_nh3 .
 
rm nte_7term.dat;        ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/nte_7term.dat .
 
rm nh3.dat;              ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/nh3.dat    .
rm so2.dat;              ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/so2.dat    .
rm n2o.dat;              ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/n2o.dat    .
rm hno3.dat;             ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/hno3.dat   .
rm co2.dat;              ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/co2.dat    .
 
rm tunmlt.txt;           ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/tunmlt.txt .
rm therm.dat;            ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/therm.dat  .
rm fx.txt;               ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/fx.txt     .

rm optran.dat;           ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/optran.dat .
rm set7.dat;             ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/set7.dat   .
rm set6.dat;             ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/set6.dat   .
rm set5.dat;             ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/set5.dat   .
rm set4.dat;             ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/set4.dat   .
rm set3.dat;             ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/set3.dat   .
rm set2.dat;             ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/set2.dat   .
rm set1.dat;             ln -s /home/sergio/asl/rta/sarta_coef/Data_CrIS_oct16/Coef/set1.dat   .

cd ../../../
## this is from /home/sergio/asl/rtp/iasi/iasi1/ecmwf/2018/iasi1_ecmwf_d20180121_clear.rtp_2
echo "from /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac  test on chip using eg "
echo "     ../bin/jac_sarta_iasi_jan25_H2020_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4 fin=../TEST_RTP/newdayx_1_100_12150_iasi.op.rtp_2 fout=junk.rp.rtp"

