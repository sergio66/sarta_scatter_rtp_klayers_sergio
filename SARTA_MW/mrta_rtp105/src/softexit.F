c       coee = routine name
c       typeidx specified as XXYZZ
c
c       Y = step ID number
c
c       XX: retrieval module  |  ZZ: type of error
c       --- ----------------  |  --- -------------
c       00: not specified     |  00: normal exit (airsb)
c       01: RETAMSU           |  01: 
c       02: RETMHS            |  02: used a bad channel
c       03: RETSURF           |  03: option not implemented
c       04: RETTMPC           |  04: namelist error
c       05: RETWATC           |  05: I/O file error: open or read
c       06: RETOZON           |  06: decoding input file error
c       07: ETACHECK          |  07: dimension exceeded error
c       08: CLOUDHGT          |  08: illogical value obtained, code error
c       09: CLOUDMAT          |  09: highly singular matrix
c       10: RET_CO2           |  10: tqli -- too many iterations
c       11: RET_CH4           |  11: invalid iphys
c       12: RET_CO            |  12: file already exists (write only modes)
c       13: NOAA REG          |  13: invalid symmetric inverse
c       14: CO2_PER           !  14: problems in ERR2SLAB inverse


      subroutine soft_decode(unit, code, typeidx, errtype)
      implicit none

c     input variables
      integer*4 unit        ! print unit
      character*8 code      ! routine name
      integer*4   typeidx   ! coded error type

c     output variables
      integer*4   errtype   ! stop return code

c     local variables
      integer*4  MAXMOD, MAXERR
      parameter  (MAXMOD=16, MAXERR=16)
      integer*4   module, istep, i1000, i100, i10
      parameter    (i1000=1000, i100=100, i10=10)
      character*8 modnme(MAXMOD)
      character*20 errnme(MAXERR)

      data modnme/'RETAMSUC','RETMHS  ','RETSURFC','RETTMPC ',
     &            'RETWATC ','RETOZONC','ETACHECK','CLOUDHGT',
     &            'CLOUDMAT','RET_CO2 ','RET_CH4 ','RET_CO  ',
     &            'REGRESS ','CO2_PER ',' ??     ',' ??      '/


      data errnme/'                  ?? ','used a bad chl       ',
     &            'option does not exist','namelist error       ',
     &            'I/O file error       ','input file error     ',
     &            'exceeded dimension   ','illogical value      ',
     &            'highly singular matrx','tqli -- iterations   ',
     &            'invalid iphys        ','file already exists  ',
     &            'inverse of ncv matrix','inverse of [F*F^T]   ',
     &            '                  ?? ','                  ?? '/

        module  = 0
        errtype = typeidx
        if(typeidx.ge.i100) then
          istep  = typeidx/i100       ! XXY
          module = istep/i10          ! XX
          istep  = istep - i10*module ! Y
          errtype = typeidx - module*i1000 - istep*i100 ! ZZ
        endif

        if(errtype.ge.2.and.errtype.le.MAXERR) then
          if(module.ge.1.and.module.le.MAXMOD) then
            write(unit,210) modnme(module),istep,code,
     &           errtype, errnme(errtype)
          else
            write(unit,200) code, errtype, errnme(errtype)
          endif
        else
          write(unit,100) errtype
        endif
      return
  100 format('softexit: invalid exit code =',i6)
  200 format('softexit from ',a8,' error code =',i6,' (',a20,')')
  210 format('softexit from ',a8,'(',i1,') via ',a8,' error code =',
     &      i2,' (',a20,')')
      end

c     -----------------------------------------------------------------

      subroutine softexit(code,typeidx)
      implicit none
      integer*4 typeidx
      character*8 code

      integer*4   unit
      integer*4   errtype

      integer*4  MAXERR
      parameter  (MAXERR=16)

      close(20)

      if(typeidx.eq.0) then      ! all OK, but pass the value back

        stop 0

      else

        unit = 6
        call soft_decode(unit, code, typeidx, errtype)

        if(errtype.ge.1.and.errtype.le.MAXERR) then
          if(errtype.eq.1) then
            stop 1
          elseif(errtype.eq.2) then
            stop 2
          elseif(errtype.eq.3) then
            stop 3
          elseif(errtype.eq.4) then
            stop 4
          elseif(errtype.eq.5) then
            stop 5
          elseif(errtype.eq.6) then
            stop 6
          elseif(errtype.eq.7) then
            stop 7
          elseif(errtype.eq.8) then
            stop 8
          elseif(errtype.eq.9) then
            stop 9
          elseif(errtype.eq.10) then
            stop 10
          elseif(errtype.eq.11) then
            stop 11
          elseif(errtype.eq.12) then
            stop 12
          elseif(errtype.eq.13) then
            stop 13
          elseif(errtype.eq.14) then
            stop 14
          elseif(errtype.eq.15) then
            stop 15
          elseif(errtype.eq.16) then
            stop 16
          endif
        else
          stop 99
        endif
      endif
      close(6)

      return
      end

