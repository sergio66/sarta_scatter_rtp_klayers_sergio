fip = '/asl/s1/sergio/rtp/j1_ccast_hires/allfov//2025/01/08/cloudy_airs_l1c_ecm_sarta_baum_ice.2025.01.08.099.rtp';
fip = 'cloudy_airs_l1c_ecm_sarta_baum_ice.2025.01.08.099.rtp';

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD

toptsSARTA.klayers    = '/asl/packages/klayersV205/BinV201/klayers_airs';
toptsSARTA.sarta_cris = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new';       %% HITRAN 2016 spectroscopy with jacobians nanananana
toptsSARTA.sarta_cris = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new'; %% HITRAN 2016 spectroscopy with jacobians nanananana
toptsSARTA.sarta_airs = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';                        %% HITRAN 2020 spectroscopy with jacobians nanananana

toptsSARTA.PBL.klayers    = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_PBL/Feb05_2025/klayersV205/BinV201/klayers_pbl_wetwater_test';
toptsSARTA.PBL.sarta_cris = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_feb25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new_PBL';
toptsSARTA.PBL.sarta_airs = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_feb25_H2020_PBL';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_clear_test
disp('ret to continue'); pause

fip = 'cloudy_airs_l1c_ecm_sarta_baum_ice.2025.01.08.099.rtp';
do_cloud_test
disp('ret to continue'); pause

do_jac_test
disp('ret to continue'); pause

do_clear_test_airs
