        program trail
c        tries out program to calc. \Gamma_1 derivatives
        implicit real*8(a-h,o-z)
        open(1,file='gamder.in')
        open(3,file='gamder.out')

C       ASSUMES THAT THE EOS is in a file called eosfile
        do 100 i=1,20000
C       read(1,*,end=101) xr, xx,c,rho, p, ggam1, t, dela, x, z
        read(1,*,end=101) xr, t, p, rho, x, z, ggam1
        call gammder(t,p,rho,x,z,ggam1,pt1, gam1,gg, dgam1y,dgam1rho,
     1 dgam1p)
c        write(3,555)xr,p,pt1,ggam1,gam1,gg, dgam1p,dgam1rho, dgam1y
        write(3,555)xr,dgam1p/gam1,dgam1rho/gam1, dgam1y/gam1
100        continue
101        continue
555        format(1p9e25.15)

        stop
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     OPAL EQUATION OF STATE
C
c        Modified by SB to read OPAL 2005 EOS (Dec 2016)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C TAKES IN P, RHO, T, Z, X from file and calculate gamma derivatives
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE gammder (T,P,rho,X,Z,ggam1, pt1, gam1, gg, dgam1y, 
     1 dgam1rho, dgam1p)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT LOGICAL*4(L)
C
      PARAMETER ( IVARX = 25 )
      PARAMETER (CNVS = 0.434294481D0 )
      PARAMETER (ZERO = 0.0D0)
C
      COMMON/Eeos/ ESACT,EOS(10)
C
      IRAD=1
      IORDER=9
      T6=T/1.0D6
      P12=P/1.0D12
      IF(T6.LT.0.0020D0.OR.T6.GT.100.0D0) THEN
C        WRITE(*,*) 'OUT OF OPAL EOS TABLE'
         RETURN 
        END IF
      ZTAB=Z
C
C
cc        print*,' calling rhoofp'
         R=RHOOFP(X,ZTAB,T6,P12,IRAD)
c        print*, 'call done'
         IF(R.LE.-998.0D0)THEN
           WRITE(*,*)'OUT OF OPAL EOS TABLE',T,P
c            RETURN 1
         ENDIF
         D=R 
         D10=DLOG10(D)
c        find changes w.r.t. x, z, t, rho
c        w.r.t. X
        xr=x-0.05
        xu=x+0.05
c        print*, 'first ca;ll'
        call esac(xr, ztab, t6, rho, iorder, irad) 
c        print*, 'first call done'
        pt1l=1.d12*eos(1)
        chrh=eos(6)
        cht=eos(7)/1.d6
        cv=eos(4)/1.d6
        pt1=1.d12*eos(1)
        cp=cv+(pt1*cht*cht)/(rho*T*chrh)
        gg=(cp/cv)*chrh
        gam1xl=eos(8)
c        gam1xl=gg
        call esac(xu, ztab, t6, rho, iorder, irad) 
        pt1u=1.d12*eos(1)
        gam1xu=eos(8)
        chrh=eos(6)
        cht=eos(7)/1.d6
        cv=eos(4)/1.d6
        pt1=1.d12*eos(1)
        cp=cv+(pt1*cht*cht)/(rho*T*chrh)
        gg=(cp/cv)*chrh
c        gam1xu=gg
c-------------------------
c        w.r.t. T
        iorder=9
        t6l=t6*0.99d0
        call esac(x, ztab, t6l, rho, iorder, irad)
        gam1tl=eos(8)
        tt=t6l*1.d6
        chrh=eos(6)
        cht=eos(7)/1.d6
        cv=eos(4)/1.d6
        pt1=1.d12*eos(1)
        cp=cv+(pt1*cht*cht)/(rho*tt*chrh)
        gg=(cp/cv)*chrh
c        gam1tl=gg
        t6r=t6*1.01d0
        tt=t6r*1.d6
        call esac(x, ztab, t6r, rho, iorder, irad)
        gam1tu=eos(8)
        chrh=eos(6)
        cht=eos(7)/1.d6
        cv=eos(4)/1.d6
        pt1=1.d12*eos(1)
        cp=cv+(pt1*cht*cht)/(rho*tt*chrh)
        gg=(cp/cv)*chrh
c        gam1tu=gg
c---------------------
c        w.r.t rho
        rhol=rho*0.99d0
        call esac(x, ztab, t6, rhol, iorder, irad)
        gam1rl=eos(8)
        chrh=eos(6)
        cht=eos(7)/1.d6
        cv=eos(4)/1.d6
        pt1=1.d12*eos(1)
        cp=cv+(pt1*cht*cht)/(rhol*T*chrh)
        gg=(cp/cv)*chrh
c        gam1rl=gg

        rhor=rho*1.01d0
        call esac(x, ztab, t6, rhor, iorder, irad)
        gam1ru=eos(8)
        chrh=eos(6)
        cht=eos(7)/1.d6
        cv=eos(4)/1.d6
        pt1=1.d12*eos(1)
        cp=cv+(pt1*cht*cht)/(rhor*T*chrh)
        gg=(cp/cv)*chrh
c        gam1ru=gg
c------------------------------------
c        call eos normally
        call esac(x, ztab, t6, rho, iorder, irad)
        gam1=eos(8)

        chrh=eos(6)
        cht=eos(7)/1.d6
        cv=eos(4)/1.d6
        pt1=1.d12*eos(1)
        cp=cv+(pt1*cht*cht)/(rho*T*chrh)
        gg=(cp/cv)*chrh
        

        dpt=eos(7)*pt1/T
        dpr=eos(6)*pt1/rho
        dpx=(pt1u-pt1l)/0.1

        dgam1y=(gam1xu-gam1xl)/0.1
        dgam1t=(gam1tu-gam1tl)/(0.02*T)
        dgam1rho=(gam1ru-gam1rl)/(0.02*rho)

        dgam1p=pt1*dgam1t/dpt
        dgam1rho=rho*(dgam1rho-dgam1t*dpr/dpt)
        dgam1y=dgam1y+dgam1t*dpx/dpt



      RETURN
c  999 WRITE (20, *)'WARNNING... OPAL TBL FAIL'
c      RETURN 1
      END

      subroutine esac (xh,ztab,t6,r,iorder,irad, *)
        implicit real*8(a-h,o-z)
c..... The purpose of this subroutine is to interpolate 
c      the equation of state and its derivatives in X, T6, density
c        izi=0 recalulate table indices to use; =1 keep previous

c        xh=hydrogen mass fraction
c        ztab is metal fraction of the EOS5_data tables you are using.
c           included only for purpose of preventing mismatch
c        t6=T6=temperature in millions of degrees kelvin
c        r=rho=Rho=density(g/cm**3)
c..... to use esac insert common/eeos/ esact,eos(10) in the calling routine.
c      This common contains the interpolated EOS values for the EOS
c
c..... eos(i) are obtained from a quadradic interpolation at
c      fixed T6 at three values of Rho; followed by quadratic
c      interpolation along T6. Results smoothed by mixing
c      overlapping quadratics.
c         definitions:
c     
c            T6 is the temperature in units of 10**6 K
c 
c            rho is the density in grams/cc
c            R=Rho/T6**3

c            eos(1) is the pressure in megabars (10**12dyne/cm**2)
c not in use eos(2) is energy in 10**12 ergs/gm. Zero is zero T6
c not in use eos(3) is the entropy in units of energy/T6
c chg -      eos(4) is dE/dRHO at constant T6 
c chg eos(4) eos(5) is the specific heat, dE/dT6 at constant V.
c chg eos(5) eos(6) is dlogP/dlogRho at constant T6. 
c                   Cox and Guil1 eq 9.82
c chg eos(6) eos(7) is dlogP/dlogT6 at conxtant Rho.
c                   Cox and Guil1 eq 9.81
c chg eos(7) eos(8) is gamma1. Eqs. 9.88 Cox and Guili.
c chg eos(8) eos(9) is gamma2/(gaamma2-1). Eqs. 9.88 Cox and Guili
c chg eos(9) discarded (gamma3-1. Eqs 9.88 Cox and Guili)

c            iorder sets maximum index for eos(i);i.e., iorder=1
c                   gives just the pressure
c
c            irad  if =0 no radiation correction; if =1 adds radiation

c            index(i),i=1,10  sets order in which the equation of state
c            variables are stored in eos(i).  Above order corresponds
c            to block data statement:
c                 data (index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/. 
c            If you, for example, only want to return gamma1: set iorder=1 
c            and set: data (index(i),i=1,10)/8,2,3,4,5,6,7,1,9,10/
c
c
      save
      integer w
      real*8 moles
      parameter (mx=5,mv=10,nr=169,nt=197)
      character blank*1
      common/lreadco/itime
      common/eeeos/ epl(mx,nt,nr),xx(mx)
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),  
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos/zz(mx),nta(nr),iri(10),index(10)
      common/bbeos/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/eeos/esact,eos(mv)
      dimension frac(7)
      data aprop/83.14511d0/
      blank=' '
      if (iorder .gt. 9 ) then
         write (*,'(" iorder cannot exceed 9")')
      endif
      if ((irad .ne. 0) .and. (irad .ne. 1)) then
         write (*,'(" Irad must be 0 or 1")')
         stop
      endif

      xxi=xh
      ri=r
c
      slt=t6
      slr=r
c
      if(itime .ne. 12345678) then
        itime=12345678
        do i=1,10
          do j=1,10
            if  (index(i) .eq. j) iri(i)=j
          enddo
        enddo
        do  i=1,mx
          xx(i)=xa(i)
        enddo
c
c..... read the data files
        call readcoeos
        z=zz(1)
c        if (abs(ztab-z) .gt. 0.00001) then
c        write (*,'("requested z=",f10.6," EOSdata is for z=",f10.6)')
c     x    ztab,z      
c        stop
c        endif

        if(z+xh-1.d-6 .gt. 1 ) go to 61
      endif
c
c
c..... Determine T6,rho grid points to use in the
c      interpolation.
      if((slt .gt. t6a(1)).or.(slt .lt. t6a(nt))) go to 62
      if((slr .lt. rho(1)).or.(slr .gt. rho(nr))) go to 63
c
c
c
c     if (izi .eq. 0) then  ! freeze table indices if not 0
        ilo=2
        ihi=mx
    8   if(ihi-ilo .gt. 1) then
          imd=(ihi+ilo)/2
            if(xh .le. xa(imd)+1.d-7) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 8
        endif
        i=ihi
        mf=i-2
        mg=i-1
        mh=i
        mi=i+1
        mf2=mi
        if (xh .lt. 1.d-6) then
        mf=1
        mg=1
        mh=1
        mi=2
        mf2=1
        endif
        if((xh .le. xa(2)+1.d-7) .or. (xh .ge. xa(mx-2)-1.d-7)) mf2=mh
c
        ilo=2
        ihi=nr
   12     if(ihi-ilo .gt. 1) then
          imd=(ihi+ilo)/2
           if (slr .eq. rho(imd)) then
           ihi=imd
           go to 13
           endif
            if(slr .le. rho(imd)) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 12
          endif
   13     i=ihi
        l1=i-2
        l2=i-1
        l3=i
        l4=l3+1
        iqu=3
        if(l4 .gt. nr) iqu=2
c
        ilo=nt
        ihi=2
   11     if(ilo-ihi .gt. 1) then
          imd=(ihi+ilo)/2
           if (t6 .eq. t6list(1,imd)) then
           ilo=imd
           go to 14
           endif
            if(t6 .le. t6list(1,imd)) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 11
          endif
   14     i=ilo
        k1=i-2
        k2=i-1
        k3=i
        k4=k3+1
        ipu=3
        if (k4 .gt. nt) ipu=2
      if (k3 .eq. 0) then
      write (20,'(" ihi,ilo,imd",3i5)')
      endif

c     endif

c     check to determine if interpolation indices fall within
c     table boundaries.  choose largest allowed size.
      sum1=0.0
      sum2=0.0
      sum23=0.0
      sum33=0.0
      do m=mf,mf2
        do ir=l1,l1+1
          do it=k1,k1+1
            sum1=sum1+xz(m,1,it,ir)
          enddo
        enddo
        do ir=l1,l1+2
          do it=k1,k1+2
            sum2=sum2+xz(m,1,it,ir)
          enddo
        enddo
        if (ipu .eq. 3) then
          do ir=l1,l1+2
            do it=k1,k1+ipu
              sum23=sum23+xz(m,1,it,ir)
            enddo
          enddo
        else
          sum23=2.d+30
        endif
        if (iqu .eq. 3) then
          do ir=l1,l1+3
            do it=k1,k1+ipu
              sum33=sum33+xz(m,1,it,ir)
            enddo
          enddo
        else
          sum33=2.d+30
        endif
      enddo
      iq=2
      ip=2
      if (sum2 .gt. 1.d+30) then
        if (sum1 .lt. 1.d+25 ) then
          k1=k3-3
          k2=k1+1
          k3=k2+1
          l1=l3-3
          l2=l1+1
          l3=l2+1
          go to 15
            else
          go to 65
        endif
      endif
      if (sum23 .lt. 1.d+30) ip=3
      if (sum33 .lt. 1.d+30) iq=3

      if(t6 .ge. t6list(1,2)+1.d-7) ip=2
c      if(slr .le. rho(2)+1.d-7) iq=2
      if(slr .le. rho(2)+1.d-15) iq=2

      if((l3 .eq. nr) .or. (k3 .eq. nt)) then
         iq=2
         ip=2
      endif

   15 continue
      do 124 iv=1,iorder
      do 123 m=mf,mf2
 
      is=0
 
c__________
      do ir=l1,l1+iq
        do it=k1,k1+ip
        epl(m,it,ir)=xz(m,iv,it,ir)
        is=1
        enddo
      enddo
  123 continue
      if((zz(mg) .ne. zz(mf)) .or. (zz(mh) .ne. zz(mf))) then
        write(20,'("Z does not match Z in EOS5_data files you are"
     x ," using")')
        stop
      endif
      if(z .ne. zz(mf)) go to 66
      is=0
      iw=1
      do 45 ir=l1,l1+iq
        do it=k1,k1+ip
          if (mf2 .eq. 1) then
          esk(it,ir)=epl(mf,it,ir)
          go to 46
          endif
          esk(it,ir)=quadeos(is,iw,xh,epl(mf,it,ir),epl(mg,it,ir)
     x    ,epl(mh,it,ir),xx(mf),xx(mg),xx(mh))
          if(esk(it,ir) .gt. 1.d+20) then
          write(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir
     x    ,l3,k3,iq,ip
          write(*,'(3e12.4)')  (epl(ms,it,ir),ms=mf,mf+2)
          endif
          is=1
   46     continue
        enddo
   45 continue

      if (mi .eq. mf2) then  ! interpolate between quadratics
      is=0
      iw=1
       dixr=(xx(mh)-xh)*dfsx(mh)
      do 47 ir=l1,l1+iq
        do it=k1,k1+ip
          esk2(it,ir)=quadeos(is,iw,xh,epl(mg,it,ir),epl(mh,it,ir)
     x    ,epl(mi,it,ir),xx(mg),xx(mh),xx(mi))
          if(esk(it,ir) .gt. 1.d+20) then
          write(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir
     x    ,l3,k3,iq,ip
          write(*,'(3e12.4)')  (epl(ms,it,ir),ms=mg,mg+2)
          endif
          esk(it,ir)=esk(it,ir)*dixr+esk2(it,ir)*(1.-dixr)
          is=1
        enddo
   47 continue
   

      endif

      is=0
c
c..... completed X interpolation. Now interpolate T6 and rho on a
c      4x4 grid. (t6a(i),i=i1,i1+3),rho(j),j=j1,j1+3)).Procedure
c      mixes overlapping quadratics to obtain smoothed derivatives.
c
c
      call t6rinteos(slr,slt)
      eos(iv)=esact
  124 continue

      p0=t6*r
      eos(iri(1))=eos(iri(1))*p0   ! interpolated in p/po
      eos(iri(2))=eos(iri(2))*t6   ! interpolated in E/T6
      tmass=gmass(xh,z,moles,eground,fracz,frac)
      if (irad .eq. 1) then
c      call radsub (t6,r,moles,tmass,iorder)
      call radsub (t6,r,moles,tmass,iorder)
      else
c      eos(iri(4))=eos(iri(4))*moles*aprop/tmas
      eos(iri(5))=eos(iri(5))*moles*aprop/tmass
      endif
      return

   61 write(20,'(" Mass fractions exceed unity (61)")')
        WRITE(20,*)Z, XH
      stop
   62 write(20,'(" T6 outside of table range (62)")')
        WRITE(20,*)T6A(1), SLT, T6A(NT)
        RETURN 1
c      stop
   63 write(20,'(" LogR outside of table range (63)")')
        WRITE(20,*) RHO(1),SLR, RHO(NR)
      RETURN 1
c      stop
   65 write(20,'("T6/log rho in empty region od table (65)")')
      write (20,'("xh,t6,r=", 3e12.4)') xh,t6,r
         RETURN 1
c      stop
   66 write(20,'(" Z does not match Z in EOS5_data* files you are",
     . " using (66)")')
      write (20,'("mf,zz(mf)=",i5,e12.4)') mf,zz(mf)
      write (20,'("  iq,ip,k3,l3,xh,t6,r,z= ",4i5,4e12.4)')
     x ip,iq,k3,l3,xh,t6,r,z
      stop
      return
      end

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C MU
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     SUBROUTINE MU(T,P,D,X,Z,RMU,AMU,EMU,BETA)
C     IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT LOGICAL*4(L)
C     DIMENSION ATOMWT(4)
C     COMMON/COMP/XENV,ZENV,ZENVM,AMUENV,FXENV(12),XNEW,ZNEW,STOTAL,SENV
C     COMMON/CTLIM/ATIME(11),TCUT(5),TSCUT,TENV0,TENV1,TENV,TGCUT,LNEUTR
C     COMMON/CONST1/ CLN,CLNI,C4PI,C4PIL,C4PI3L,CC13,CC23
C     COMMON/CONST2/CGAS,CA3,CA3L,CSIG,CSIGL,CGL,CMKH,CMKHN
C     COMMON/CCOUT2/LDEBUG,LCORR,LMILNE,LTRACK,LSTORE,LSTPCH,NPUNCH
C     COMMON/LUOUT/ILAST,IDEBUG,ITRACK,ISHORT,IMILNE,IMODPT,ISTOR
C     DATA NZP1/12/
C     DATA ATOMWT/0.9921D0,0.24975D0,0.08322D0,0.4995D0/
C     DATA SCALET,SCALES/0.500D0, 2.000D0/
C SET UP FRACTIONAL ABUNDANCES
C     DFX1 = (X - XENV)
C     DFX12 = (Z - ZENV)
C     IF(DABS(DFX1) + DABS(DFX12) .LT. 1.0D-5) THEN
C USE ENVELOPE ABUNDANCES
C      AMU = AMUENV
C     ELSE
C      DFX1 = DFX1*ATOMWT(1)
C      DFX12 = DFX12*ATOMWT(3)
C      DFX4 = (XENV+ZENV -X-Z)*ATOMWT(2)
C ASSUME EXCESS Z(METALS) IS IN THE FORM OF CARBON(12)
C      AMU = AMUENV + DFX1 + DFX4 + DFX12
C     ENDIF
C
C     EE=((BETA*P)/(D*T*CGAS*AMU))-1.0D0
C     EMU=EE*AMU
C     RMU=CGAS*(AMU+EMU)
C     RETURN
C     END
C*************************************************************************

      subroutine t6rinteos(slr,slt)
        implicit real*8(a-h,o-z)
c     The purpose of this subroutine is to interpolate in T6 and rho
      save
      parameter (mx=5,mv=10,nr=169,nt=197)
      common/eeeos/ epl(mx,nt,nr),xx(mx)
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),  
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/bbeos/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/eeos/esact,eos(mv)
c
      iu=0
      is=0

      do kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quadeos(is,iw,slr,esk(kx,l1),esk(kx,l2),esk(kx,l3),
     x  rho(l1),rho(l2),rho(l3))
          if(iq. eq. 3) then
            iw=2
            q(iu)=quadeos(is,iw,slr,esk(kx,l2),esk(kx,l3),esk(kx,l4),
     x      rho(l2),rho(l3),rho(l4))
          endif
        is=1
      enddo
c
      is=0
      iw=1
c..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      esact=quadeos(is,iw,slt,h(1),h(2),h(3),t6a(k1),t6a(k2),t6a(k3))
        if (iq. eq. 3) then
c.....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          esactq=quadeos(is,iw,slt,q(1),q(2),q(3),t6a(k1),t6a(k2),
     x           t6a(k3))
        endif
        if(ip .eq. 3) then
c.....    eos(i) in lower-left 3x3.
          esact2=quadeos(is,iw,slt,h(2),h(3),h(4),t6a(k2),t6a(k3),
     x           t6a(k4))
c.....    eos(i) smoothed in left 3x4
          dix=(t6a(k3)-slt)*dfs(k3)
          esact=esact*dix+esact2*(1.-dix)
c       endif   ! moved to loc a
        if(iq .eq. 3) then
 
c.....     eos(i) in upper-right 3x3.
          esactq2=quadeos(is,iw,slt,q(2),q(3),q(4),t6a(k2),t6a(k3),
     x            t6a(k4))
          esactq=esactq*dix+esactq2*(1.-dix)
        endif
        endif  ! loc a
c
        if(iq .eq. 3) then
          dix2=(rho(l3)-slr)*dfsr(l3)
            if(ip .eq. 3) then
c.....        eos(i) smoothed in both log(T6) and log(R)
              esact=esact*dix2+esactq*(1-dix2)
            endif
        endif
        if (esact .gt. 1.d+15) then
          write(*,'("Interpolation indices out of range",
     x              ";please report conditions.")') 
          stop
        endif

      return
      end

c
c**********************************************************************
      subroutine readcoeos
        implicit real*8(a-h,o-z)
c..... The purpose of this subroutine is to read the data tables
      save
      parameter (mx=5,mv=10,nr=169,nt=197)
      LOGICAL*4 LOPALE,LOPALD
      CHARACTER*80 FOPALE
c      COMMON/OPALEOS/FOPALE,LOPALE,LOPALD,IOPALE

      real*8 moles
      character*1 blank
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),  
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos/zz(mx),nta(nr),iri(10),index(10)
      common/eeos/esact,eos(mv)
      common/eeeos/ epl(mx,nt,nr),xx(mx)
      common/eeeeos/moles(mx),xin(mx),tmass(mx),
     x    amu_M(mx,nr),alogNe(mx,nr),rhogr(mx,nr),frac(mx,6),
     x    alogr(nr,nt),icycuse(mx,nr)



      blank=' '


        if (itimeco .ne. 12345678) then
        do i=1,mx
          do j=1,mv 
            do k=1,nt
              do l=1,nr
                xz(i,j,k,l)=1.d+35
              enddo
            enddo
          enddo
        enddo
        itimeco=12345678
        endif
 
      close (2)
c..... read  tables
c       call system (' gunzip EOSdata')
       iopale=2
       open(IOPALE, FILE='eosfile')
 
 
      do 3 m=1,mx

      read (IOPALE,'(3x,f6.4,3x,f12.9,11x,f10.7,17x,f10.7)')
     x  xin(m),zz(m),moles(m),tmass(m)
      read (IOPALE,'(21x,e14.7,4x,e14.7,3x,e11.4,3x,e11.4,3x,e11.4,
     x 4x,e11.4)') (frac(m,i),i=1,6)
      read (IOPALE,'(a)') blank
      do 2 jcs=1,nr
      read (IOPALE,'(2i5,2f12.7,17x,e15.7)') numtot,icycuse(m,jcs)
     x  ,dum,dum,rhogr(m,jcs)
      if(numtot .ne. jcs) then
         write (*,'(" Data file incorrect: numtot,jcs= ",2i5)')
     x     numtot,jcs
         stop
      endif
      read(IOPALE,'(a)') blank
      read(IOPALE,'(a)') blank
      if (icycuse(m,jcs) .lt. nta(jcs)) then
         write (*,'("problem with data files: X=",f6.4," density=",
     x      e14.4)') xin(m), rhogr(m,jcs)
         stop
      endif
      do  i=1,icycuse(m,jcs)
      if (i .gt. nta(jcs)) then
         read (IOPALE,'(a)') blank
         go to 4
      endif
      read (IOPALE,'(f11.6,1x,f6.4,e11.4,2e13.6,2e11.3,5f10.6)')
     x     t6list(jcs,i),amu_M(jcs,i),alogNe(jcs,j),
     x     (xz(m,index(iv),i,jcs),iv=1,9)
    4 continue
      enddo
      read(IOPALE,'(a)') blank
      read(IOPALE,'(a)') blank
      read(IOPALE,'(a)') blank
    2 continue
      read(IOPALE,'(a)') blank
    3 continue
 
      do i=1,nt
         if(t6list(1,i) .eq. 0.0) then
            write(*,'("READCOEOS: Error:",i4,
     $           "-th T6 value is zero")') i
            stop
         endif
         t6a(i)=t6list(1,i)
      enddo
      do 12 i=2,nt
   12 dfs(i)=1./(t6a(i)-t6a(i-1))
      rho(1)=rhogr(1,1)
      do 13 i=2,nr
      rho(i)=rhogr(1,i)
   13 dfsr(i)=1./(rho(i)-rho(i-1))
      do i=2,mx
      dfsx(i)=1./(xx(i)-xx(i-1))
      enddo
      close(IOPALE)
c       call system ('gzip EOSdata')
      return
      end
c
c***********************************************************************
      real*8 function quadeos(ic,i,x,y1,y2,y3,x1,x2,x3)
        implicit real*8(a-h,o-z)
c..... this function performs a quadratic interpolation.
      save
      dimension  xx(3),yy(3),xx12(30),xx13(30),xx23(30),xx1sq(30)
     . ,xx1pxx2(30)
      xx(1)=x1
      xx(2)=x2
      xx(3)=x3
      yy(1)=y1
      yy(2)=y2
      yy(3)=y3
        if(ic .eq. 0) then
          xx12(i)=1./(xx(1)-xx(2))
          xx13(i)=1./(xx(1)-xx(3))
          xx23(i)=1./(xx(2)-xx(3))
          xx1sq(i)=xx(1)*xx(1)
          xx1pxx2(i)=xx(1)+xx(2)
        endif
      c3=(yy(1)-yy(2))*xx12(i)
      c3=c3-(yy(2)-yy(3))*xx23(i)
      c3=c3*xx13(i)
      c2=(yy(1)-yy(2))*xx12(i)-(xx1pxx2(i))*c3
      c1=yy(1)-xx(1)*c2-xx1sq(i)*c3
      quadeos=c1+x*(c2+x*c3)
      return
      end
c******************************************************************8
      real*8 function gmass(x,z,amoles,eground,fracz,frac)
        implicit real*8(a-h,o-z)
      dimension anum(6),frac(7), amas(7),eion(7)
      data (eion(i),i=1,6)/-3394.873554,-1974.86545,-1433.92718,
     x  -993.326315,-76.1959403,-15.29409/
      data (anum(i),i=1,6)/10.,8.,7.,6.,2.,1./
      xc=0.247137766
      xn=.0620782
      xo=.52837118
      xne=.1624188
      amas(7)=1.0079
      amas(6)=4.0026
      amas(5)=12.011
      amas(4)=14.0067
      amas(3)=15.9994
      amas(2)=20.179
      amas(1)=0.00054858
      fracz=z/(xc*amas(5)+xn*amas(4)+xo*amas(3)+xne*amas(2))
      xc2=fracz*xc
      xn2=fracz*xn
      xo2=fracz*xo
      xne2=fracz*xne 
      xh=x/amas(7)
      xhe=(1.-x -z)/amas(6)
      xtot=xh+xhe+xc2+xn2+xo2+xne2
      frac(6)=xh/xtot
      frac(5)=xhe/xtot
      frac(4)=xc2/xtot
      frac(3)=xn2/xtot
      frac(2)=xo2/xtot
      frac(1)=xne2/xtot
      eground=0.0
      amoles=0.0
      do i=1,6
      eground=eground+eion(i)*frac(i)
      amoles=amoles+(1.+anum(i))*frac(i)
      enddo
      anume=amoles-1.
      gmass=0.
      do i=2,7
      gmass=gmass+amas(i)*frac(i-1)
      enddo

      return
      end

c***********************************************************************
      subroutine radsub (t6,density,moles,tmass,iorder)
      implicit real*8(a-h,o-z)
      parameter (mx=5,mv=10,nr=169,nt=197)
      real*8 moles,k,molenak,Na
      common/eeos/esact,eos(mv)
      common/beos/zz(mx),nta(nr),iri(10),index(10)

      data Na/6.0221367d+23/, k/1.380658d-16/, unitf/0.9648530/, 
     x unitfold/0.965296/, c/2.9979245d+10/, sigma/5.67051d-5/
     x , sigmac/1.8914785d-15/, sigmacc/1.8914785d-3/, aprop/83.14510d0/

cPhysical constants
c       Na=6.0221367d+23
c       k =1.380658d-16 !   erg/degree K
c       Na*k=6.0221367d+23*1.380658d-16 erg/degree K=8.314510d+7 erg/degree K
c           =8.314510d+7*11604.5 erg/eV=0.9648575d+12 erg/eV
c           Define unitf= Na*k/d+12=0.9648575
c           unitf=0.9648530  ! obtained with latest physical constants
c           unitfold=0.965296   ! old units- still used in the EOS code
c           In these units energy/density is in units of Mb-CC/gm
c           Pressure is in units of d+12 bars=Mb
c       sigma is the Stefan-Boltzmann constant
c       sigma=5.67051d-5 !   erg /(s*cm**2*K**4)
c       c=2.99792458d+10 !   cm/sec

c     rat=sigma/c    ! dyne/(cm**2*K**4)

c     rat=rat*1.d+24  !   Convert degrees K to units 10**6 K (T replaced with T6)
      rat=sigmacc

      molenak=moles*aprop  ! Mb-cc/unit T6
      pr=4.d0/3.d0*rat*t6**4   ! Mb 
      er=3.d0*pr/density   ! Mb-cc/gm
      sr=4.d0/3.d0*er/t6   ! Mb-cc/(gm-unit T6)

c-----Calculate EOS without radiation correction
      pt=eos(iri(1))
      et=eos(iri(2))
      st=eos(iri(3))
      dedrho=eos(iri(4))
      chir=eos(iri(6))*eos(iri(1))/pt
      chitt=(eos(iri(1))*eos(iri(7)))/pt
      cvtt=(eos(iri(5))*molenak/tmass)
      if(iorder .ge. 5) gam3pt_norad=pt*chitt/(cvtt*density*t6)
      gam1t_norad=chir+chitt*gam3pt_norad
      if(iorder .ge. 5) gam2pt_norad=gam1t_norad/gam3pt_norad
c---- End  no radiation calculation

c---- Calculate EOS with radiation calculation
      pt=eos(iri(1))+pr
      et=eos(iri(2))+er
      st=eos(iri(3))+sr
      dedrho=eos(iri(4))-er/density
      chir=eos(iri(6))*eos(iri(1))/pt
      chitt=(eos(iri(1))*eos(iri(7))+4.*pr)/pt
c     gam1t(jcs,i)=(p(jcs,i)*gam1(jcs,i)+4./3.*pr)/pt(jcs,i)
c     gam2pt(jcs,i)=(gam2p(jcs,i)*p(jcs,i)+4.*pr)/pt(jcs,i)
c     gam3pt(jcs,i)=gam1t(jcs,i)/gam2pt(jcs,i)
      cvtt=(eos(iri(5))*molenak/tmass+4.*er/t6)
      gam3pt=pt*chitt/(cvtt*density*t6)                              ! DIRECT
      gam1t=chir+chitt*gam3pt !eq 16.16 Landau_Lifshitz (Stat. Mech) ! DIRECT
      gam2pt=gam1t/gam3pt                                            ! DIRECT

c     normalize cvt to 3/2 when gas is ideal,non-degenerate,
c     fully-ionized, and has no radiation correction
c     cvt=(eos(5)*molenak/tmass+4.*er/t6)
c    x  /molenak
      eos(iri(1))=pt
      eos(iri(2))=et
      eos(iri(3))=st
      eos(iri(4))=dedrho
      eos(iri(5))=cvtt
      eos(iri(6))=chir
      eos(iri(7))=chitt
c-----Add difference between EOS with and without radiation.  cvtt
c       calculation is not accurate enough to give accurate results using
c       eq. 16.16 Landau&Lifshitz (SEE line labeled DIRECT)
      eos(iri(8))=eos(iri(8))+gam1t-gam1t_norad
      eos(iri(9))=eos(iri(9))+gam2pt-gam2pt_norad
c     eos(iri(10))=eos(iri(10))+gam3pt-gam3pt_norad
      return
      end

c********************************************************************
      real*8 function rhoofp(x,ztab,t6,p,irad)
        implicit real*8(a-h, o-z)
      parameter (mx=5,mv=10,nr=169,nt=197)
      common/lreadco/itime
      common/aeos/  xz(mx,mv,nt,nr),  
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos/zz(mx),nta(nr),iri(10),index(10)
      common/eeos/esact,eos(mv)
      dimension nra(nt)
      data sigmacc/1.8914785d-3/
c--------------------------------------------------------------------
      data (nra(i),i=1,nt)/ 26*169, 168, 2*167, 166, 165, 164,
     x   2*163, 162, 161, 160, 2*159, 158, 2*130, 129, 2*128,
     x   2*126, 2*125, 124, 122, 121, 2*120, 2*119, 6*118, 2*117,
     x   2*116, 2*115, 4*114, 8*113, 22*111, 5*110, 6*109, 2*108,
     x   5*107, 2*106, 105, 2*104, 8*103, 16*102, 21*100, 9* 99,
     x   6* 98, 4* 97, 95, 94, 6* 87/
      rat=sigmacc
      pr=0.0
c      if(irad .eq. 1) pr=4./3.*rat*t6**4   ! Mb 
      pnr=p-pr

      if (itime .ne. 12345678) then
         XDBG = 0.5D0
         T6DBG = 1.0D0
         RDBG = 0.001D0
         IDBG = 1
         call esac (xdbg,ztab,T6DBG,RDBG,IDBG,irad, *999)
      endif

      ilo=2
      ihi=mx
 8    if(ihi-ilo .gt. 1) then
         imd=(ihi+ilo)/2
         if(x .le. xa(imd)+1.d-7) then
            ihi=imd
         else
            ilo=imd
         endif
         go to 8
      endif
      mlo=ilo

      ilo=nt
      ihi=2
 11   if(ilo-ihi .gt. 1) then
         imd=(ihi+ilo)/2
         if (t6 .eq. t6list(1,imd)) then
            ilo=imd
            go to 14
         endif
         if(t6 .le. t6list(1,imd)) then
            ihi=imd
         else
            ilo=imd
         endif
         go to 11
      endif
 14   klo=ilo
      
      pmax=xz(mlo,1,klo,nra(klo))*t6*rho(nra(klo))+
     + IRAD*4.D0/3.D0*RAT*T6**4
      pmin=xz(mlo,1,klo,1)*t6*rho(1)+
     + IRAD*4.D0/3.D0*RAT*T6**4
      if((pnr .gt. 1.25*pmax) .or. (pnr .lt. pmin)) then
      write (*,'(" The requested pressure-temperature not in table")')
c     stop
      write (*,'("pnr, pmax,pmin=",3e14.4)') pnr,pmax,pmin
      goto 999
      endif
      
      rhog1=rho(nra(klo))*pnr/pmax
      call esac (x,ztab,t6,rhog1,1,irad, *999)
      p1=eos(1)
      if(p1 .gt. pnr) then
         p2=p1
         rhog2=rhog1
         rhog1=0.2*rhog1
         if(rhog1 .lt. 1.d-14) rhog1=1.d-14
         call esac (x,ztab,t6,rhog1,1,IRAD, *999)
         p1=eos(1)
      else
         rhog2=5.d0*rhog1
         if(rhog2 .gt. rho(klo)) rhog2=rho(klo)
         call esac (x,ztab,t6,rhog2,1,IRAD, *999)
         p2=eos(1)
      endif

      icount=0
    1 continue
      icount=icount+1
      rhog3=rhog1+(rhog2-rhog1)*(pnr-p1)/(p2-p1)
      call esac (x,ztab,t6,rhog3,1,irad, *999)
      p3=eos(1)
      if (abs((p3-pnr)/pnr) .lt. 1.d-5) then
         rhoofp=rhog3
         return
      endif
      if (p3 .gt. pnr) then
         rhog2=rhog3
         p2=p3
         if (icount .lt. 11) go to 1
c        write (*,'("No convergence after 10 tries")')
         goto 999
c        stop
      else
         rhog1=rhog3
         p1=p3
         if (icount .lt. 11) go to 1
c        write (*,'("No convergence after 10 tries")')
        Go to 999
c        stop
      endif
 999  CONTINUE
      RHOOFP=-999.0D0
C     WRITE(20,'("FAIL TO FIND RHO")')
      RETURN
      end

c***********************************************************************
      block data
      implicit real*8 (a-h,o-z)
      parameter (mx=5,mv=10,nr=169,nt=197)
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),  
     .     t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     .     ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos/zz(mx),nta(nr),iri(10),index(10)
      data (xa(i),i=1,mx)/0.0D0,0.2D0,0.4D0,0.6D0,0.8D0/

      data (index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/
      data (nta(i),i=1,nr)/87*197,7*191,190,2*189,185,179,170,2*149,
     x    133,125,123,122,120,115,113,107,102,2*80,72,68,66,64,62,56,54,
     x    52,51,2*50,49,47,2*45,43,42,28*40,39,37,36,35,34,32,31,30,
     x    29,27,26/
      end
