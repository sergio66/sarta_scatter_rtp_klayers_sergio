Rapid transmittance algorithm:

TGRS03mwrta.pdf, JGR06FM_VAL.pdf, igarss06.pdf: papers describing the algorithm & validation.
mwtran.f     high level subroutine to compute layer
              transmittances.  (3/29/06)
mwtran.inc   include file for mwtran (3/29/06)
bfield.f     magnetic field routine (needed for AMSU-A). (4/24/97)
shval3.f     subroutine called by bfield. (4/25/97)
opac2.f      calculation of opacity terms. called by mwtran,
              but it could be used independently.  (3/30/06)
vlint.f      subroutine called by opac2. (11/1/95)
getcoef.f    routine to read a coefficient file.  (7/10/06)
tb10.f       high-level subroutine to compute brightness temperature
               components from transmittances and level temperatures.
               (5/16/01)
tb11.f       high-level subroutine to compute brightness temperature
               components from transmittances and mean layer
               temperatures. (5/17/01)

When using tb10 or tb11, the surface emission and reflection and the
cosmic background contribution are to be added in the calling
program, e.g. the upward-propagating brightness temperature would be:
  TB = TD + E*( (1.-Refl)*Tsurf + Refl*(TR + ER*TBcosmic) )


Coefficient files for several microwave sounding instruments:

tr_amsua.dat AMSU-A
tr_amsub.dat AMSU-B
tr_mhs.dat   MHS
tr_atms.dat  ATMS
tr_mts1.dat  NAST-M 50-56 GHZ
tr_mts2s.dat NAST-M 118 GHZ
tr_mts3.dat  NAST-M 166/183 GHz
tr_mts4.dat  NAST-M 425 GHz

If using a case-sensitive system, unzip with the -L option.
