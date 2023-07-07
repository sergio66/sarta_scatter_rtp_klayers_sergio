function sarta = read_sarta_nlte_coef()
%{
../src/rdcoef.f
../src/incFTC.f says FNCOFN = /asl/rta/sarta_database/Data_AIRS_may19/Coef/nte_7term.dat

../src/incFTC.f:465:       PARAMETER(MXCNTE = 200)

C      ------------
C      Read non-LTE
C      ------------
       IF (COFNTE) THEN
       IF (DEBUG) write(6,*) 'rdcoef:COFNTE=TRUE read coef file'

       OPEN(UNIT=IOUN,FILE=FNCOFN,FORM='UNFORMATTED',STATUS='OLD',
     $    IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOFN
          STOP
       ENDIF
C
       J=1
       DO I=1,MXCNTE
C         Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, (COEFN(IC,J),IC=1,NNCOEF)
C
C         Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLISTN(J)=ICHAN
             J=J + 1
          ENDIF
       ENDDO
       NCHNTE=J - 1
C
       CLOSE(IOUN)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fid = fopen('/asl/rta/sarta_database/Data_AIRS_may19/Coef/nte_7term.dat','r','ieee-be');
fid = fopen('/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/NLTE/prod_2022/airs_l1c/apr2021/fitc/r49/AIRS_L1C_R49_cutcoef_xnte_hmmmm__le90_7term.dat','r','ieee-be');
mxcnte = 200;
for ii = 1 : mxcnte
  flen = fread(fid, 1, 'integer*4');
  sarta.ichan(ii) = fread(fid,1,'integer*4');
  sarta.vchan(ii) = fread(fid,1,'real*4');
  junk = fread(fid,7,'real*4');
  sarta.cofn(ii,:) = junk;
  flen = fread(fid, 1, 'integer*4');  
end
fclose(fid);
