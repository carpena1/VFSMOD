      subroutine runoff(P,CN,xIa,Q)
C---------------------------------------------------------------
C SCS runoff calculation
C---------------------------------------------------------------
      implicit double precision (a-h, o-z)

      S=25400.d0/CN-254.d0
      xIa=0.2d0*S
      if(P.gt.xIa)then
            Q=(P-xIa)**2.0d0/(P+0.8d0*S)
        else
            Q=0.d0
      endif

      return
      end

      subroutine calctc(pL,CN,Y,tc)
C---------------------------------------------------------------
C SCS concentration time calculation
C---------------------------------------------------------------
      implicit double precision (a-h, o-z)

      tc=pL**0.8d0*(1000.d0/CN-9.d0)**0.7d0/(4407.d0*dsqrt(Y))

      return
      end


      subroutine q_peak(A,Q,xIa,P,tc,j,qp,tp)
C---------------------------------------------------------------
C      SCS-TR55 peak flow calculation
C---------------------------------------------------------------
      implicit double precision (a-h, o-z)
      dimension ci(3,4,5)

      data ci/68.0317,-82.907,11.1619,144.547,-130.64,-55.230,-11.312,
     &16.6125,-43.015,-11.505,-64.177,65.9007,-74.693,105.222,-26.314,
     &-136.68,134.907,47.9565,12.1681,-16.337,50.4334,14.2182,85.7116,
     &-85.806,24.9255,-42.167,16.1126,41.8526,-45.773,-13.503,-6.5688,
     &6.4981,-19.740,-7.8919,-38.206,39.0036,-3.9797,6.7479,-2.9776,
     &-6.2829,6.585,2.1954,1.0577,-1.1784,3.2996,1.3836,6.7419,-6.8946,
     &2.5222,-0.8657,0.0456,2.3645,-0.6384,-0.2644,2.5021,-0.5476,
     &-0.3427,2.4007,-0.8899,0.2078/

c      do 10 j=1,4
c            do 10 i=1,3
c                  write(2,100) (ci(i,j,k),k=1,5)
c10      continue

c----rmc 04/20/03 - Fix for Q=0
      if(Q.le.0) then
        qp=0.d0
        tp=qp
        return
      endif
c----rmc 04/20/03 - end of fix for Q=0
      xIaP=xIa/P
c---- TR55 stablishes that if Ia/P is outside the range (0.1<Ia/P<0.5),
c-----use the limiting values
      if(xIa/P.gt.0.5d0)xIaP=0.5d0
      if(xIa/P.lt.0.1d0)xIaP=0.1d0

c--- Import Ia/P and storm type I,IA,II,III (j=1,4) ------------
      if(j.lt.5) then
         C0=ci(1,j,1)*xIaP**4+ci(1,j,2)*xIaP**3+ci(1,j,3)*xIaP**2+
     &     ci(1,j,4)*xIaP+ci(1,j,5)
         C1=ci(2,j,1)*xIaP**4+ci(2,j,2)*xIaP**3+ci(2,j,3)*xIaP**2+
     &     ci(2,j,4)*xIaP+ci(2,j,5)
         C2=ci(3,j,1)*xIaP**4+ci(3,j,2)*xIaP**3+ci(3,j,3)*xIaP**2+
     &     ci(3,j,4)*xIaP+ci(3,j,5)
        else
         C0=ci(1,3,1)*xIaP**4+ci(1,3,2)*xIaP**3+ci(1,3,3)*xIaP**2+
     &     ci(1,3,4)*xIaP+ci(1,3,5)
         C1=ci(2,3,1)*xIaP**4+ci(2,3,2)*xIaP**3+ci(2,3,3)*xIaP**2+
     &     ci(2,3,4)*xIaP+ci(2,3,5)
         C2=ci(3,3,1)*xIaP**4+ci(3,3,2)*xIaP**3+ci(3,3,3)*xIaP**2+
     &     ci(3,3,4)*xIaP+ci(3,3,5)
      endif



C---- Unit q peak, qp (m3/s) -----------------------------------
      qu=4.3046d0*10.d0**(C0+C1*dlog10(tc)+C2*(dlog10(tc))**2-6.d0)
      Fp=1.d0
      qp=qu*A*Q*Fp

c--- Unit hydrograph time to peak (min) ------------------------
      tp=0.127481d0*Q*A/qp/60.d0

100   format(5f9.3)
      return
      end

      subroutine unit_hyd(Q,A,qp,tp,D,tc,qp5,tp5,mref)
C---------------------------------------------------------------
C Unit NRCS hydrograph using Haan's equation (k=3.77)
C---------------------------------------------------------------
      implicit double precision (a-h, o-z)
      common/hydgph/u(10000,2),qh(10000,2)

      ck=3.77d0
      t=0.d0
c----Initialize u vector
      do 5 i=1,208
            do 5 j=1,2
            u(i,j)=0.d0
5     continue

c----rmc 09/30/11 obsolete from uh v2.x, scaling unit hydrograph now
c----uh v3.x hydrograph by convolution from synthetic excess hyetograph
c
c----rmc 04/20/03 - Fix for Q=0
c      if(Q.le.0) then
c          ttotal=D
c        else
c          ttotal=0.6375d0*Q*A/qp/60.d0
c      endif
c----rmc 04/20/03 - end of fix for Q=0
c      dt=ttotal/50.d0
c      do 10 i=0,51
c          t=i*dt
c          if(Q.le.0) then
c                qi=0.d0
c              else
c                qi=qp*(t/tp*dexp(1-t/tp))**ck
c          endif
c          rot(i+1)=t
c          roq(i+1)=qi
c          qdepth=qi*360.d0/A
c          write(2,100)t,qi,qdepth
c10      continue
c------ end of uh v2.x calculations----------------

c----rmc 09/15/11- New hydrograph calculation from unit hydrograph for
c----convolution (tabular) hidrograph calculation. Def:rain excess incr.(h)
c----,qp5(m3/s),tp5(h):peak t,q. Estimate time step for dimesionless unit
c----hydrograph as Def<>1/3.tp. Since tp=0.6tc+Def/2, if Def<>1/3tp then
c----Def<>0.24tc.
c-5min-  Def=5.d0/60.d0
      Def=0.24d0*tc
      write(2,205)Def*60.d0
205   format(' SCS ',f4.1,'-min unit hydrograph',/,
     1 ' time (h)    q(m3/s)   q(mm/h)',/,30('-'))
      tp5=0.6d0*tc+0.5d0*Def
      qp5=0.127481d0*A/(tp5*60.d0)
c--- SCS triangular hydrograph---
c     ttotal5=2.67d0*tp5
c----SCS adimensional unit hydrograph
      ttotal5=5.d0*tp5
c--------------------------------
      dt5=Def
      t5=0.d0
      i=0
      cqdepth5=0.d0
      do while (t5.le.ttotal5)
            t5=i*dt5
            if(Q.le.0.d0) then
                  qi5=0.d0
               else
                  qi5=qp5*(t5/tp5*dexp(1.d0-t5/tp5))**ck
            endif
            u(i+1,1)=t5
            u(i+1,2)=qi5
            qdepth5=qi5*360.d0/A
            cqdepth5=qdepth5+cqdepth5
            write(2,'(3f10.4)')(u(i+1,j),j=1,2),qdepth5
            i=i+1
      end do
c----rmc - mref, number of unit hydrograph steps needed in convolution
      mref=i
c----rmc - check unit hydrograph volumen <> 1
      unitq=cqdepth5*dt5
      if(dabs(1.d0-unitq).le.0.05d0) then
            write(2,120)'PASSED unit hydrograph check- V(mm)=',unitq
         else
            write(2,120)'FAILED unit hydrograph check- V(mm)=',unitq
      endif

c100   format(f9.2,2f10.4)
120   format(/,A36,f6.2)
      return
      end

      subroutine tab_hyd(Q,A,mref,nref,qp,nhyd)
C---------------------------------------------------------------
C Calculation of hydrograph by convolution (Chow, 1987) of SCS unit
C hydrograph and excess hyetograph
C   mref: number of unit hydrograph steps
C   nref: number of excess hyetograph steps
C---------------------------------------------------------------
      implicit double precision (a-h, o-z)
      common/hydgph/u(10000,2),qh(10000,2)
      common/rain/rfix,rti(10000),rfi(10000),rcum(10000,2),ref(10000,2),
     1 ncum

      write(2,205)mref,nref
      cqdepth5=0.d0
      do 40 i=1,mref
c            write(2,'(2f10.4)')(u(i,j),j=1,2)
            cqdepth5=u(i,2)*360.d0/A+cqdepth5
40    continue
      dt5=u(2,1)-u(1,1)
      unitq=cqdepth5*dt5
c      do 50 i=1,nref
c            write(2,'(2f10.4)')(ref(i,j),j=1,2)
c50    continue

c----Apply convolution of the u and ref values to obtained hydrograph
      Def=u(2,1)-u(1,1)
      qp=0.d0
      do 70 k=1,nref+mref-1
            qh(k,1)=ref(1,1)+(k-1)*Def
            qh(k,2)=0.d0
            do 60 i=1,k
                  qh(k,2)=qh(k,2)+ref(i,2)*u(k-i+1,2)
60          continue
c-----------For debugging, print full tabular hydrograph
c            write(*,'(300f10.4)')qh(k,1),qh(k,2)*360.d0/A,
c     &           (ref(i,2)*u(k-i+1,2)*360.d0/A,i=1,k)
c-----------End of tabular hydrograph printing
            if ( qh(k,2).gt.qp ) then
                 tpn=qh(k,1)
            end if
            qp=dmax1(qh(k,2),qp)

70    continue
      qdepth=qp*360.d0/A
      do 80 i=1,k-1
            write(2,'(300f10.4)')(qh(i,j),j=1,2),qh(i,2)*360.d0/A
80    continue
      nhyd=k-1
      write(2,1000)qp,qdepth
      write(2,1100)tpn,tpn*60.d0
      write(2,1200)nhyd


205   format('Number of unit hydrograph steps (n)  =',i5,/,
     1 'Number of excess hyetograph steps (m)=',i5,/,/,
     2 ' Final hydrograph (convolution)',/,
     3 ' time (h)    q(m3/s)   q(mm/h)',/,30('-'))
1000  format(/,'Peak flow=',f9.4,' m3/s= ',f9.4,' mm/h')
1100  format('Time to peak=',f8.2,' h= ',f8.2,' min')
1200  format('Number of final hydrograph steps (nhyd)  =',i5)

      return
      end
