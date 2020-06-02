*****&****************************************************************
***** partial wave analysis for J/psi --> K+ K- Pi0
***** date : March 27, 2005
      program pwa_050327
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      include './input/inc/const.in'
      include './input/inc/state.in'
*****&================================================================
c---- common blocks for fumili
** 01 parameters, which includes both strong interaction information
**    and coefficients for each partial wave
      common/a/faa(npa)
** 02 step length
      common/pl/fpl(npa)
** 03 lower limits
      common/al/amn(npa)
** 04 upper limits 
      common/au/amx(npa)
** 05 standard deviation
      common/sigma/sig0(npa)
** 06 covariance error matrix
**    common/z/errz(nzz)
** 07 total number of experimental point (including background)
      common/ned/nt,ns
** 08 indicator of observation
      common/exda/fex(nex)
** 09 start point
      common/uprmon/ijhr
*****&================================================================
c---- common block for amplitudes and cross section
** 01 momentum of final states, cross section at i-th experimental point
      common/check/ak1(4),ak2(4),ak3(4),prob
** 02 information from data (number, fianl state, momentum)
      common/data/akn(nd1,3,5)
** 03 cf->amplitude     kon->intermediate state
      common/cfnew/cf(nmode,2),kon(mmode),kon2
** 04 modula of partial wave amplitude
      common/trans/fu(mmode,mmode)
** 05 integrated of modula of partial wave amplitude
      common/coef/fun(mmode,mmode)
** 06 modula of partial wave amplitude for data point including background
      common/startup/funall(nd1,mmode,mmode),idcs
** 07 coefficient, d(sigma)/d(a), cross section
      common/cppp/cp(nmode),waa(npa),was
** 08 modula of coefficient
      common/pan/pa(nmode,nmode)
** 09 d(sigma)/d(a) without including background
      common/norm/wa1
** 10 data point
      common/x/x(3)
** 11 amplitude
      common/ppwa/ampx(nmode,nmode)
*****&================================================================
c---- variables for output
      real*8 vml,vmli,dam2,wti(nmode,nmode),xwt,sfx,sfy
      character hcoef(npa)*24
*****&================================================================
c---- start with general introduction
      goto 20
 10   write(*,500)
      stop
 20   write(*,510)
 500  format(/5x,'error : please choose a right parameter'/)
 510  format(/5x,50('*'),
     & /5x,'* Partail Wave Analysis for J/Psi ---> K+ K- Pi0 *',
     & /5x,'*',48x,'*',
     & /5x,'* please chose the run-stage parameter           *',
     & /5x,'* the availbale run -stage paramters are         *',
     & /5x,'* istep=0 ==> perform istep=1 and istep=2        *',
     & /5x,'* istep=1 ==> calculate fun(i,j) using MC sample *',
     & /5x,'* istep=2 ==> optimize fractions of resonances   *',
     & /5x,'* istep=3 ==> output weight of MC sample         *',
     & /5x,'*',48x,'*',
     & /5x,'* now please chose the run-stage parameter       *',
     & /5x,50('*'),/)
*****&================================================================
c---- choice of run_stage of program (see 510)
      open(1,file='./start/step.in',status='old')
      read(1,*,err=10) istep
      close(1)
C      print *,'Input which step do you want to do ?'
C      read *,istep
      write(*,520) istep
      if (istep.lt.0) goto 10
      if (istep.gt.4) goto 10
 520  format(/5x,'ok ! the stage is : ',i3/)
*****&================================================================
c---- mass,width,form-of-breit-wigner,name of resonance
      open(2,file='./input/res/res.in',status='old')
      do i1=1,nresx
        read(2,530) ip,xmass(i1),xwidth(i1),kff(i1),hres(i1)
      enddo
      close(2)
      if (nbgxx.gt.0) then
        nnnn=nresx+1
        open(3,file='./input/res/bak/bg.in',status='old')
        do i1=nnnn,nmode
          read(3,540) hres(i1)
        enddo
        close(3)
      endif
      hres(mmode)='background'
 530  format(2x,i3,f13.5,f13.5,15x,i3,2x,a10,1x)
 540  format(2x,a10)
*****&================================================================
c---- coefficients
      if (istep.eq.3) then
        open(3,file='./output/out/para.out',status='old')
      else
        open(3,file='./input/res/para.in',status='old')
      endif
      do i1=1,npa
         read(3,550) ip,faa(i1),fpl(i1),hcoef(i1),amn(i1),amx(i1)
      enddo
      close(3)
 550  format(1x,i3,f10.4,f9.3,2x,a24,f10.3,f10.3)
*****&================================================================
c---- initialization of parameters
      do i1=1,mmode
      do i2=1,mmode
        fu(i1,i2)=0.0
      enddo
      enddo
      ijhr=0
      idcs=0
      kon2=0
      call metric
      call conv(cp)
*****&================================================================
      if (istep.eq.3) then
        write(*,560)
        open(4,file='./output/out/fun.out',status='old')
        do i1=1,mmode
        if (kon(i1).eq.1) then
          do i2=1,mmode
          if (kon(i2).eq.1) then
            read(4,*) fun(i1,i2)
          endif
          enddo
        endif
        enddo
        close(4)
        call eval
        write(*,570)
        stop
      endif
 560  format(/5x,51('*'),
     &  /8x,'istep=3',8x,'scale MC Sample with weight !',
     &  /5x,51('*'),/)
 570  format(/5x,51('*'),
     &  /5x,'* this step has finished (*^_^*)                  *'
     &  /5x,'* result has been written into output/out/reweight.wt   *'
     &  /5x,51('*'),/)
*****&================================================================
      if (istep.le.1) then
        write(*,580)
        call integ(1)
        open(5,file='./output/out/fun.out',status='unknown')
        do i1=1,mmode
        if (kon(i1).eq.1) then
          do i2=1,mmode
          if (kon(i2).eq.1) then
            write(5,*) fun(i1,i2)
          endif
          enddo
        endif
        enddo
        close(5)
        write(*,590)
        if (istep.eq.0) then
          istep=2
          goto 40
        endif
        stop
      endif
 580  format(/5x,51('*'),
     &  /8x,'istep=1',8x,'begin to calculate fun(i,j) !',
     &  /5x,51('*'),/)
 590  format(/5x,51('*'),
     &  /5x,'* this step has finished (*^_^*)                  *'
     &  /5x,'* result has been written into output/out/fun.out *'
     &  /5x,51('*'),/)
*****&================================================================
      if (istep.eq.2) then
        write(*,700)
        open(4,file='./output/out/fun.out',status='old')
        do i1=1,mmode
        if (kon(i1).eq.1) then
          do i2=1,mmode
          if (kon(i2).eq.1) then
            read(4,*) fun(i1,i2)
          endif
          enddo
        endif
        enddo
        close(4)
      endif
 40   continue
      do i1=1,mmode
      if (kon(i1).eq.1) then
        do i2=1,mmode
        if (kon(i2).eq.1) then
          fun(i1,i2)=fun(i1,i2)*1.0
        else
          fun(i1,i2)=0.0
        endif
        enddo
      endif
      enddo
 700  format(/5x,51('*'),
     &  /8x,'istep=2',8x,'begin to optimize coefficients !',
     &  /5x,51('*'),/)
*****&================================================================
***   read four-momentum of final states of data events
c--   i2=1,2,3    : K+, K-, pi0
c--   i3=1,2,3,4  : Px, Py, P, E
      write(*,710)
      i5=ndt
      open(10,file='./data/PWA_dt.dat',status='old')
      do i1=1,i5
      do i2=1,3
        read(10,*) (akn(i1,i2,i3),i3=1,4)
      enddo
      enddo
      close(10)
      if (nbg1.gt.0) then
        i4=i5+1
        i5=i5+nbg1
        open(11,file='./data/bg1.dat',status='old')
        do i1=i4,i5
        do i2=1,3
          read(11,*) (akn(i1,i2,i3),i3=1,4)
        enddo
        enddo
        close(11)
      endif
      if (nbg2.gt.0) then
        i4=i5+1
        i5=i5+nbg2
        open(12,file='./data/bg2.dat',status='old')
        do i1=i4,i5
        do i2=1,3
          read(12,*) (akn(i1,i2,i3),i3=1,4)
        enddo
        enddo
        close(12)
      endif
      if (nbg3.gt.0) then
        i4=i5+1
        i5=i5+nbg3
        open(13,file='./data/bg3.dat',status='old')
        do i1=i4,i5
        do i2=1,3
          read(13,*) (akn(i1,i2,i3),i3=1,4)
        enddo
        enddo
        close(13)
      endif
 710  format(/5x,'now reading four-momentum of final states',
     &       /5x,'please be patience ! ....................',/)
*****&================================================================
***   calculate amplitude at each data point
      write(*,720)
      do i1=1,ntot
        x(1)=i1
        temp=funct(x)
        do i2=1,nmode
        if (kon(i2).eq.1) then
          do i3=1,nmode
          if (kon(i3).eq.1) then
            funall(i1,i2,i3)=fu(i2,i3)
          endif
          enddo
        endif
        enddo
        if (kon(mmode).eq.1) funall(i1,mmode,mmode)=fu(mmode,mmode)
      enddo
      idcs=1
 720  format(/5x,'now calculating amplitude at each data point',
     &       /5x,'please be patience ! .......................',/)
*****&================================================================
***   using maximum likelihood fit to optimize coefficients
      write(*,730)
      it=1 ! printing control
      ns=1 ! used by subroutine sgz(m,z)
      n1=3
      n2=2     ! times for successful fit before changing parameter
      n3=150   ! iteraction limit
      epfu=0.1 ! tolerance to control termination criterion
      nt=ntot
      do i1=1,ntot
        fex(i1)=i1
      enddo
      call likelm(ss,npa,n1,n2,n3,epfu,aka,ala,it,ms)
 730  format(/5x,'now optimizing coefficients of all amplitude',
     &       /5x,'please be patience ! .......................',/)
*****&================================================================
***   output the optimized coefficients
      write(*,740)
      open (4,file='./output/out/para.out',status='unknown')
      open (5,file='./output/out/coef.err',status='unknown')
      do i1=1,npa
        write(4,550) i1,faa(i1),fpl(i1),hcoef(i1),amn(i1),amx(i1)
        if (fpl(i1).lt.0.0) then
          write(5,750) i1,sig0(i1)
        else
          write(5,760) i1,sig0(i1)
        endif
      enddo
      close(4)
      close(5)
      vml=ss
c      vml=0.0
c      do i1=1,ntot
c        x(1)=i1
cc---  normalized differential cross section
c        dam2=funct(x)
c        if (i1.le.ndt) then
c          vmli=+dlog(dam2)
c        elseif(i1.gt.ndt.and.i1.le.(ndt+nbg1)) then 
c          vmli=-1.80*dlog(dam2)
c        elseif(i1.gt.(ndt+nbg1).and.i1.le.(ndt+nbg1+nbg2)) then 
c          vmli=-0.00*dlog(dam2)
c        elseif(i1.gt.(ndt+nbg1+nbg2).and.i1.le.(ntot)) then 
c          vmli=-0.00*dlog(dam2)
c        endif
c        vml=vml-vmli
c      enddo
 740  format(/5x,'now computing the value of maximum likelihood',
     &       /5x,'corresponding to the optimized coefficients',
     &       /5x,'and outputing the result',/)
 750  format(1x,i3,1x,'fixed',1x,f12.5)
 760  format(1x,i3,7x,f12.5)
*****&===============================================================
***   output of fitting result
      open(5,file='./output/out/like.out',status='unknown')
      write(5,800) ms
c--- normal return on finding minimum
      if (ms.eq.+1) write(5,810)
c--- no further decrease in objective function
      if (ms.eq.-1) write(5,820)
c--- error matrix can not be inverted
      if (ms.eq.-2) write(5,830)
c--- iteraction limit is exceeded
      if (ms.eq.-3) write(5,840)
c--- objective function is -ln(f) with f<0
      if (ms.eq.-4) write(5,850)
c--- maximum-likelihood,cross-section,data,signal,MC-sample
      write(5,860) vml,was,ntot,nsig,nmc
 800  format(7x,57('*')
     $     /7x,'*',24x,'summary',24x,'*'
     $     /7x,'*',55x,'*'
     &     /7x,'* the optimization situation is :',4x,i4,15x,'*')
 810  format(7x,'* normal return (*^_^*)',33x,'*')
 820  format(7x,'* not decrease in objective function more',15x,'*')
 830  format(7x,'* error matrix can not be inverted',22x,'*')
 840  format(7x,'* iteration limit is exceeded',27x,'*')
 850  format(7x,'* -ln(f) with objective function f<0',19x,'*')
 860  format(7x,'* maximum likelihood estimation: ',f20.3,3x,'*',
     &     /7x,'* value of total denominator   : ',f20.3,3x,'*',
     &     /7x,'* (data) signal + background   : ',i20,  3x,'*',
     &     /7x,'*',8x,'signal',i12,4x,'MC',i14,9x,'*',
     &     /7x,'*',55x,'*',
     $     /7x,57('*'),///3x,'i',3x,'j',3x,
     &     'percent',7x,'events',3x,'(meson-1) x (meson-2)'/)
*****&===============================================================
***   output of information on resonance
      do i1=1,nmode
      do i2=1,nmode
        wti(i1,i2)=0.0
      enddo
      enddo
      do i1=1,nmode
      if (kon(i1).eq.1) then
        do i2=1,nmode
        if (kon(i2).eq.1) then
          wti(i1,i2)=pa(i1,i2)*fun(i1,i2)
        endif
        enddo
      endif
      enddo
      xwt=0.0
      sfx=100.0/dble(nmc)/was
      sfy=dble(nsig)/dble(nmc)/was
      do i1=1,nmode
      if (kon(i1).eq.1) then
        write(5,*) ''
        do i2=1,nmode
        if (kon(i2).eq.1) then
          if (i1.eq.i2) then
            tmp=wti(i1,i2)
            xwt=xwt+tmp
            write(5,870) i1,i2,tmp*sfx,tmp*sfy,hres(i1),hres(i2)
          elseif (i1.lt.i2) then
            tmp=wti(i1,i2)+wti(i2,i1)
            xwt=xwt+tmp
            write(5,870) i1,i2,tmp*sfx,tmp*sfy,hres(i1),hres(i2)
          endif
        endif
        enddo
      endif
      enddo
      write(5,880)
      open(6,file='./output/out/res.out',status='unknown')
      do i1=1,nmode
        if (kon(i1).eq.1) then
          tmp=wti(i1,i1)
          write(5,890) i1,hres(i1),tmp*sfx,tmp*sfy
          write(6,900) i1,hres(i1),1,tmp,tmp*sfx,tmp*sfy
        else
          write(6,900) i1,hres(i1),0,0.0,0.0,0.0
        endif
      enddo
      if (kon(mmode).eq.1) then
        tmp=faa(npa)*fun(mmode,mmode)
        xwt=xwt+tmp
        write(5,890) mmode,hres(mmode),tmp*sfx,tmp*sfy
        write(6,900) mmode,hres(mmode),1,tmp,tmp*sfx,tmp*sfy
      else
        write(6,900) mmode,hres(mmode),0,0.0,0.0,0.0
      endif
      write(5,910) xwt*sfx,xwt*sfy
      write(6,920) xwt
      close(5)
      close(6)
 870  format(2(1x,i3),2x,f8.3,2x,f11.2,2x,a10,' & ',a10)
 880  format(//70('=')//5x,'no.',7x,'meson',3x,'percent',7x,'events'/)
 890  format(5x,i3,2x,a10,1x,":",f8.3,1x,":",f11.2)
c 900  format(1x,i3,2x,a10,2x,i3,2x,f20.3,2x,f8.3,2x,f11.2)
 900  format(1x,i3,2x,a10,1x,":",i3,1x,":",f20.3,1x,":",f8.3,2x,f11.2)
 910  format(/5x,'sum',14x,f8.3,2x,f11.2)
 920  format(23x,f20.3)
*****&================================================================
      return
      end
*****&****************************************************************


*****&****************************************************************
**    coefficient or weight of partial wave
*****&================================================================
      subroutine conv(cp)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      include './input/dat/res.dat'
      include './input/inc/state.in'
      common/a/faa(npa)
      common/pl/fpl(npa)
      common/cfnew/cf(nmode,2),kon(mmode),kon2
      complex*16 cp(nmode),ci
      integer mk1,mk2,mk3,mk4,mr1,mr3,mbg
      real*8 phi
*****&================================================================
      ci=cmplx(0.0,1.0)
c      do i1=1,2*nkvxx,2
c         i3=i1*2
c         i2=i3-1
c         i4=i1+1
c         phi=faa(i3)*0.0174532925 ! pi/180=0.0174532925199
c         cp(i1)=faa(i2)*(cos(phi)+ci*sin(phi))
c         cp(i4)=cp(i1)
c      enddo
      do i1=1,2*nkvxx,1
         i3=i1*2
         i2=i3-1
         i4=i1+1
c         phi=faa(i3)*0.0174532925 ! pi/180=0.0174532925199
         phi=faa(i3)
         cp(i1)=faa(i2)*(cos(phi)+ci*sin(phi))
c         cp(i4)=cp(i1)
      enddo
      j1=2*nkvxx+1
      j2=2*nkvxx+nrvxx+navxx+nbgxx
      do i1=j1,j2
         i3=i1*2
         i2=i3-1
c         phi=faa(i3)*0.0174532925 ! pi/180=0.0174532925199
         phi=faa(i3)
         cp(i1)=faa(i2)*(cos(phi)+ci*sin(phi))
      enddo
      if (kon2.ne.0) return
*****&================================================================
      kon2=1;mk1=0;mk2=0;mk3=0;mk4=0;mr1=0;mr3=0;ma2=0;ma4=0;mbg=0
      do i1=1,2*nkvxx
         i3 = i1*2
         i2 = i3-1
         if (((fpl(i2).lt.0.0).and.(fpl(i3).lt.0.0)).or.
     &        ((faa(i2).eq.0.0).and.(faa(i3).eq.0.0)))    then
            kon(i1)=0
         else
            kon(i1)=1
            
            kr=2*nkv1m
            if (i1.le.kr) then
               mk1=mk1+1
            else
               kr=kr+2*nkv2p
               if (i1.le.kr) then
                  mk2=mk2+1
               else
                  kr=kr+2*nkv3m
                  if (i1.le.kr) then
                     mk3=mk3+1
                  else
                     mk4=mk4+1 
                  endif
               endif
            endif
         endif
      enddo
      mk1 = mk1/2
      mk2 = mk2/2
      mk3 = mk3/2
      mk4 = mk4/2 

      j1=2*nkvxx+1
      j2=2*nkvxx+nrvxx+navxx+nbgxx
      do i1=j1,j2
         i3=i1*2
         i2=i3-1
         if (((fpl(i2).lt.0.0).and.(fpl(i3).lt.0.0)).or.
     &        ((faa(i2).eq.0.0).and.(faa(i3).eq.0.0))) then
            kon(i1)=0
         else
            kon(i1)=1
            kr=2*nkvxx+nrv1m
            if (i1.le.kr) then
               mr1=mr1+1
            else
               kr=kr+nrv3m
               if (i1.le.kr) then
                  mr3=mr3+1
               else
                  kr=kr+nav2p
                  if(i1.le.kr) then
                     ma2=ma2+1 
                  else
                     kr=kr+nav4p
                     if(i1.le.kr) then
                        ma4=ma4+1
                     else
                        mbg=mbg+1
                     endif
                  endif
               endif
            endif
         endif
      enddo
      kon(mmode)=1
      if ((fpl(npa).lt.0.0).or.(faa(npa).eq.0.0)) kon(mmode)=0
*****&================================================================
      write(*,100)
      do i1=1,nmode
         if (kon(i1).eq.1) write(*,110) i1,hres(i1)
      enddo
      kr=mk1+mk2+mk3+mk4+mr1+mr3+ma2+ma4+mbg
      write(*,120) kr,mk1,mk2,mk3,mk4,mr1,mr3,ma2,ma4,mbg
 100  format(/5x,'the resonances considered in this fitting are',/)
 110  format(15x,i5,10x,a10)
 120  format(/5x,'summary:  D*  D2*  D3*  D4*   X1   X3   X2,  X4   bg',
     &       /5x,i4,2x,9(i5))
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)


*****&****************************************************************
***   differential cross section
*****&================================================================
      function dcs(pa)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      include './input/inc/const.in'
      include './input/inc/tensor.in'
      include './input/inc/state.in'
      common/x/x(3)
      common/a/faa(npa)
      common/cppp/cp(nmode),waa(npa),was
      common/check/ak1(4),ak2(4),ak3(4),prob
      common/data/akn(nd1,3,5)
      common/cfnew/cf(nmode,2),kon(mmode),kon2
      common/trans/fu(mmode,mmode)
      common/coef/fun(mmode,mmode)
      common/startup/funall(nd1,mmode,mmode),idcs
      common/ppwa/ampx(nmode,nmode)
      real*8 pa(nmode,nmode),wa,ampx
*****&================================================================
      if (idcs.eq.1) then
         wa=0.0
         ii=x(1)
         do i1=1,nmode
            if (kon(i1).eq.1) then
               do i2=1,nmode
                  if (kon(i2).eq.1) then
                     fu(i1,i2)=funall(ii,i1,i2)
                     wa=wa+pa(i1,i2)*fu(i1,i2)
                  endif
               enddo
            endif
         enddo
         fu(mmode,mmode)=funall(ii,mmode,mmode)
         if (wa.lt.0.0) wa=0.0
         dcs=wa
         return
      endif
*****&================================================================
** momentum and mass of intermediate states
** 1 : positive charged   ak1,sp1 --> Ks    akx1,sx1 --> K*+-
** 2 : negative charged   ak2,sp2 --> K+-   akx2,sx2 --> K*0
** 3 : neutral            ak3,sp3 --> pi+-  akx3,sx3 --> X0
      do i=1,4
         akx1(i)=ak3(i)+ak1(i)
         akx2(i)=ak2(i)+ak3(i)
         akx3(i)=ak1(i)+ak2(i)
      enddo
      sp1 = scalar(ak1, ak1)  
      sp2 = scalar(ak2, ak2)  
      sp3 = scalar(ak3, ak3)  
      sx1 = scalar(akx1,akx1) 
      sx2 = scalar(akx2,akx2) 
      sx3 = scalar(akx3,akx3) 
*****&================================================================
** Eq.(13) in EPJA16p537(03) [hep-ph/0211457]
      q2jx1 = q2abc(sj,sx1,sp2)
      q2jx2 = q2abc(sj,sx2,sp1)
      q2jx3 = q2abc(sj,sx3,sp3)
      q2x1 = q2abc(sx1,sp1,sp3)
      q2x2 = q2abc(sx2,sp2,sp3)
      q2x3 = q2abc(sx3,sp1,sp2)
** Eq.(2) in EPJA16p537(03) [hep-ph/0211457]
c---- g~_uv = {-1,-1,-1,0} for J/psi
      do i=1,4
         do j=1,4
            delx1(i,j) = del(i,j)-akx1(i)*akx1(j)/sx1 
            delx2(i,j) = del(i,j)-akx2(i)*akx2(j)/sx2 
            delx3(i,j) = del(i,j)-akx3(i)*akx3(j)/sx3 
         enddo
      enddo
** Eq.(10) in EPJA16p537(03) [hep-ph/0211457]
      do i=1,3                 
         t1jx1(i) = -ak2(i)    
         t1jx2(i) = -ak1(i)    
         t1jx3(i) = -ak3(i)    
      enddo
      do i=1,4
         t1xx1(i)=0.0
         t1xx2(i)=0.0
         t1xx3(i)=0.0
         do j=1,4
            t1xx1(i)=t1xx1(i)+delx1(i,j)*del(j,j)*(ak3(j)-ak1(j))
            t1xx2(i)=t1xx2(i)+delx2(i,j)*del(j,j)*(ak2(j)-ak3(j))
            t1xx3(i)=t1xx3(i)+delx3(i,j)*del(j,j)*(ak1(j)-ak2(j))
         enddo
      enddo
** term (r.r) in Eq.(11) Eq.(12) in EPJA16p537(03) [hep-ph/0211457]
      if ((nresx-nkv1m-nrv1m).gt.0)  then
         r2jx1= scalar(t1jx1,t1jx1) 
         r2jx2= scalar(t1jx2,t1jx2) 
         r2jx3= scalar(t1jx3,t1jx3) 
         r2xx1= scalar(t1xx1,t1xx1) 
         r2xx2= scalar(t1xx2,t1xx2) 
         r2xx3= scalar(t1xx3,t1xx3) 
      endif
c---- 1- K*
      if(nkv1m.gt.0) then
** Eq.(14) in EPJA16p537(03) [hep-ph/0211457]
         b1jx1 = bwbar1(q2jx1,fud1)*bwbar1(q2x1,fud2)
         b1jx2 = bwbar1(q2jx2,fud1)*bwbar1(q2x2,fud2)
** Eq.(31) in EPJA16p537(03) [hep-ph/0211457]
         do i1=1,2
            temp1=0.0
            temp2=0.0
            do i2=1,3
               if (i2.ne.i1) then
                  do i3=1,3
                     tmp=eps(i1,i2,i3,4)
                     if (tmp.ne.0.0) then
                        temp1=temp1+tmp*amj*t1jx1(i2)*t1xx1(i3)
                        temp2=temp2+tmp*amj*t1jx2(i2)*t1xx2(i3)
                     endif
                  enddo
               endif
            enddo

            do i2=1,(2*nkv1m),2
               if (kon(i2).eq.1) then
                  i3=i2+1
                  if (kff(i2).eq.1) then
                     cf(i2,i1)=temp1*b1jx1*cbwp1(i2,sx1)
                     cf(i3,i1)=temp2*b1jx2*cbwp1(i3,sx2)
                  elseif (kff(i2).eq.2) then
                     cf(i2,i1)=temp1*b1jx1*cbwp2(i2,1,sx1,sp1,sp3)
                     cf(i3,i1)=temp2*b1jx2*cbwp2(i3,1,sx2,sp2,sp3)
                  endif
               endif
            enddo
** for background of Ks pi, K pi
c            if ( kon(23).eq.1 ) then
c               cf(23,i1)=temp1*b1jx1+temp2*b1jx2
c            endif
         enddo
      endif
c---- 2+ K*
      if (nkv2p.gt.0) then
** Eq.(11) in EPJA16p537(03) [hep-ph/0211457]
         do i1=1,3
            do i2=1,3
               t2jx1(i1,i2) = 0.0
               t2jx2(i1,i2) = 0.0
               t2xx1(i1,i2) = 0.0
               t2xx2(i1,i2) = 0.0
            enddo
         enddo

         b2jx1 = 0.0
         b2jx2 = 0.0

         do i1=1,3
            do i2=1,3
               t2jx1(i1,i2)=t1jx1(i1)*t1jx1(i2)-r2jx1*del(i1,i2)/3.0
               t2jx2(i1,i2)=t1jx2(i1)*t1jx2(i2)-r2jx2*del(i1,i2)/3.0
               t2xx1(i1,i2)=t1xx1(i1)*t1xx1(i2)-r2xx1*delx1(i1,i2)/3.0
               t2xx2(i1,i2)=t1xx2(i1)*t1xx2(i2)-r2xx2*delx2(i1,i2)/3.0
            enddo
         enddo
** Eq.(15) in EPJA16p537(03) [hep-ph/0211457]
         b2jx1 = bwbar2(q2jx1,fud1)*bwbar2(q2x1,fud2)
         b2jx2 = bwbar2(q2jx2,fud1)*bwbar2(q2x2,fud2)
         do i1=1,2
            temp1=0.0
            temp2=0.0
            do i2=1,3
               if (i2.ne.i1) then
                  do i3=1,3
                     tmp=eps(i1,i2,i3,4)
                     if (tmp.ne.0.0) then
                        do i4=1,3
                           temp1=temp1+tmp*amj*t2jx1(i2,i4)*t2xx1(i3,i4)
                           temp2=temp2-tmp*amj*t2jx2(i2,i4)*t2xx2(i3,i4)
                        enddo
                     endif
                  enddo
               endif
            enddo
            i4=2*nkv1m+1
            i5=2*(nkv1m+nkv2p)
            do i2=i4,i5,2
               if (kon(i2).eq.1) then
                  i3=i2+1
                  if (kff(i2).eq.1) then
                     cf(i2,i1)=temp1*b2jx1*cbwp1(i2,sx1)
                     cf(i3,i1)=temp2*b2jx2*cbwp1(i3,sx2)
                  elseif (kff(i2).eq.2) then
                     cf(i2,i1)=temp1*b2jx1*cbwp2(i2,2,sx1,sp1,sp3)
                     cf(i3,i1)=temp2*b2jx2*cbwp2(i3,2,sx2,sp2,sp3)
                  endif
               endif
            enddo
** for background of Ks pi, K+- pi 
c            if ( kon(30).eq.1 ) then
c               cf(30,i1)=temp1*b2jx1 + temp2*b2jx2
c            endif
         enddo
      endif
c---- 3- K*
      if (nkv3m.gt.0) then
** Eq.(12) in EPJA16p537(03) [hep-ph/0211457]
         do i1=1,3
            do i2=1,3
               do i3=1,3
                  t3jx1(i1,i2,i3)=t1jx1(i1)*t1jx1(i2)*t1jx1(i3)
     &                 -r2jx1/5.0*
     $                 ( del(i1,i2)*t1jx1(i3)
     &                 + del(i2,i3)*t1jx1(i1)
     &                 + del(i3,i1)*t1jx1(i2) )
                  t3jx2(i1,i2,i3)=t1jx2(i1)*t1jx2(i2)*t1jx2(i3)
     &                 -r2jx2/5.0*
     $                 ( del(i1,i2)*t1jx2(i3)
     &                 + del(i2,i3)*t1jx2(i1)
     &                 + del(i3,i1)*t1jx2(i2) )
                  t3xx1(i1,i2,i3)=t1xx1(i1)*t1xx1(i2)*t1xx1(i3)
     &                 -r2xx1/5.0*
     $                 ( delx1(i1,i2)*t1xx1(i3)
     &                 + delx1(i2,i3)*t1xx1(i1)
     &                 + delx1(i3,i1)*t1xx1(i2) )
                  t3xx2(i1,i2,i3)=t1xx2(i1)*t1xx2(i2)*t1xx2(i3)
     &                 -r2xx2/5.0*
     $                 ( delx2(i1,i2)*t1xx2(i3)
     &                 + delx2(i2,i3)*t1xx2(i1)
     &                 + delx2(i3,i1)*t1xx2(i2) )
               enddo
            enddo
         enddo
** Eq.(16) in EPJA16p537(03) [hep-ph/0211457]
         b3jx1 = bwbar3(q2jx1,fud1)*bwbar3(q2x1,fud2)
         b3jx2 = bwbar3(q2jx2,fud1)*bwbar3(q2x2,fud2)
** Eq.(32) in EPJA16p537(03) [hep-ph/0211457]
         do i1=1,2
            temp1=0.0
            temp2=0.0
            do i2=1,3
               if (i2.ne.i1) then
                  do i3=1,3
                     tmp=eps(i1,i2,i3,4)
                     if (tmp.ne.0.0) then
                        do i4=1,3
                           do i5=1,3
                              temp1=temp1+tmp*amj*t3jx1(i2,i4,i5)
     $                             *t3xx1(i3,i4,i5)
                              temp2=temp2+tmp*amj*t3jx2(i2,i4,i5)
     $                             *t3xx2(i3,i4,i5)
                           enddo 
                        enddo
                     endif
                  enddo
               endif
            enddo
            i4=2*(nkv1m+nkv2p)+1
            i5=2*(nkv1m+nkv2p+nkv3m)
            do i2=i4,i5,2
               if (kon(i2).eq.1) then
                  i3=i2+1
                  if (kff(i2).eq.1) then
                     cf(i2,i1)=temp1*b3jx1*cbwp1(i2,sx1)
                     cf(i3,i1)=temp2*b3jx2*cbwp1(i3,sx2)
                  elseif (kff(i2).eq.2) then
                     cf(i2,i1)=temp1*b3jx1*cbwp2(i2,3,sx1,sp1,sp3)
                     cf(i3,i1)=temp2*b3jx2*cbwp2(i3,3,sx2,sp2,sp3)
                  endif
               endif
            enddo
** for background of Ks pi, K+- pi 
c            if ( kon(27).eq.1 ) then
c               cf(27,i1)=temp1*b3jx1 + temp2*b3jx2
c            endif
         enddo
      endif
c---- 4+ K*
      if (nkv4p.gt.0) then
** Eq.(37) in PRD48p1229(93)
        do i1=1,3
           do i2=1,3
              do i3=1,3
                 do i4=1,3
                    t4jx1(i1,i2,i3,i4)=t1jx1(i1)*t1jx1(i2)*t1jx1(i3)
     $                   *t1jx1(i4)-r2jx1/7.0*
     $                   ( del(i1,i2)*t1jx1(i3)*t1jx1(i4)
     &                   + del(i1,i3)*t1jx1(i2)*t1jx1(i4)
     &                   + del(i1,i4)*t1jx1(i2)*t1jx1(i3)
     &                   + del(i2,i3)*t1jx1(i1)*t1jx1(i4)
     &                   + del(i2,i4)*t1jx1(i1)*t1jx1(i3)
     &                   + del(i3,i4)*t1jx1(i1)*t1jx1(i2) )
     &                   +(r2jx1**2)/35.0*
     $                   ( del(i1,i2)*del(i3,i4)
     &                   + del(i1,i3)*del(i2,i4)
     &                   + del(i1,i4)*del(i2,i3) )
                    t4jx2(i1,i2,i3,i4)=t1jx2(i1)*t1jx2(i2)*t1jx2(i3)
     $                   *t1jx2(i4)-r2jx2/7.0*
     $                   ( del(i1,i2)*t1jx2(i3)*t1jx2(i4)
     &                   + del(i1,i3)*t1jx2(i2)*t1jx2(i4)
     &                   + del(i1,i4)*t1jx2(i2)*t1jx2(i3)
     &                   + del(i2,i3)*t1jx2(i1)*t1jx2(i4)
     &                   + del(i2,i4)*t1jx2(i1)*t1jx2(i3)
     &                   + del(i3,i4)*t1jx2(i1)*t1jx2(i2) )
     &                   +(r2jx2**2)/35.0*
     $                   ( del(i1,i2)*del(i3,i4)
     &                   + del(i1,i3)*del(i2,i4)
     &                   + del(i1,i4)*del(i2,i3) )
                    t4xx1(i1,i2,i3,i4)=t1xx1(i1)*t1xx1(i2)*t1xx1(i3)
     $                   *t1xx1(i4)-r2xx1/7.0*
     $                   ( delx1(i1,i2)*t1xx1(i3)*t1xx1(i4)
     &                   + delx1(i1,i3)*t1xx1(i2)*t1xx1(i4)
     &                   + delx1(i1,i4)*t1xx1(i2)*t1xx1(i3)
     &                   + delx1(i2,i3)*t1xx1(i1)*t1xx1(i4)
     &                   + delx1(i2,i4)*t1xx1(i1)*t1xx1(i3)
     &                   + delx1(i3,i4)*t1xx1(i1)*t1xx1(i2) )
     &                   +(r2xx1**2)/35.0*
     $                   ( delx1(i1,i2)*delx1(i3,i4)
     &                   + delx1(i1,i3)*delx1(i2,i4)
     &                   + delx1(i1,i4)*delx1(i2,i3) )
                    t4xx2(i1,i2,i3,i4)=t1xx2(i1)*t1xx2(i2)*t1xx2(i3)*
     $                   t1xx2(i4)-r2xx2/7.0*
     $                   ( delx2(i1,i2)*t1xx2(i3)*t1xx2(i4)
     &                   + delx2(i1,i3)*t1xx2(i2)*t1xx2(i4)
     &                   + delx2(i1,i4)*t1xx2(i2)*t1xx2(i3)
     &                   + delx2(i2,i3)*t1xx2(i1)*t1xx2(i4)
     &                   + delx2(i2,i4)*t1xx2(i1)*t1xx2(i3)
     &                   + delx2(i3,i4)*t1xx2(i1)*t1xx2(i2) )
     &                   +(r2xx2**2)/35.0*
     $                   ( delx2(i1,i2)*delx2(i3,i4)
     &                   + delx2(i1,i3)*delx2(i2,i4)
     &                   + delx2(i1,i4)*delx2(i2,i3) )
                 enddo
              enddo
           enddo
        enddo
** Eq.(17) in EPJA16p537(03) [hep-ph/0211457]
        b4jx1 = bwbar4(q2jx1,fud1)*bwbar4(q2x1,fud2)
        b4jx2 = bwbar4(q2jx2,fud1)*bwbar4(q2x2,fud2)
        do i1=1,2
           temp1=0.0
           temp2=0.0
           do i2=1,3
              if (i1.ne.i2) then
                 do i3=1,3
                    tmp=eps(i1,i2,i3,4)
                    if (tmp.ne.0.0) then
                       do i4=1,3
                          do i5=1,3
                             do i6=1,3
                                temp1=temp1+tmp*amj*t4jx1(i2,i4,i5,i6)
     $                               *t4xx1(i3,i4,i5,i6)
                                temp2=temp2-tmp*amj*t4jx2(i2,i4,i5,i6)
     $                               *t4xx2(i3,i4,i5,i6)
                             enddo
                          enddo
                       enddo
                    endif
                 enddo
              endif
           enddo
           i4=2*(nkv1m+nkv2p+nkv3m)+1
           i5=2*nkvxx
           do i2=i4,i5,2
              if (kon(i2).eq.1) then
                 i3=i2+1
                 if (kff(i2).eq.1) then
                    cf(i2,i1)=temp1*b4jx1*cbwp1(i2,sx1)
                    cf(i3,i1)=temp2*b4jx2*cbwp1(i3,sx2)
                 elseif (kff(i2).eq.2) then
                    cf(i2,i1)=temp1*b4jx1*cbwp2(i2,4,sx1,sp1,sp3)
                    cf(i3,i1)=temp2*b4jx2*cbwp2(i3,4,sx2,sp2,sp3)
                 endif
              endif
           enddo
        enddo
      endif
c---- 1- rho
      if (nrv1m.gt.0) then
** Eq.(14) in EPJA16p537(03) [hep-ph/0211457]
         b1jx3 = bwbar1(q2jx3,fud1)*bwbar1(q2x3,fud2)
** Eq.(28) in EPJA16p537(03) [hep-ph/0211457]
         do i1=1,2
            temp3=0.0
            do i2=1,3
               if (i2.ne.i1) then
                  do i3=1,3
                     tmp=eps(i1,i2,i3,4)
                     if (tmp.ne.0.0) then
                        temp3=temp3+tmp*amj*t1jx3(i2)*t1xx3(i3)
                     endif
                  enddo
               endif
            enddo
            i4=2*nkvxx+1
            i5=2*nkvxx+nrv1m
            do i2=i4,i5
              IF(I2.NE.23) THEN
                if (kon(i2).eq.1) then
                  if (kff(i2).eq.1) then
                    cf(i2,i1)=temp3*b1jx3*cbwp1(i2,sx3)
                  elseif (kff(i2).eq.2) then
                    cf(i2,i1)=temp3*b1jx3*cbwp2(i2,1,sx3,sp1,sp2)
                  endif
                endif
              ENDIF
            enddo
** for background of Ks K+- 
            if ( kon(23).eq.1 ) then
              cf(23,i1)=temp3*b1jx3
            endif
         enddo
      endif
c---- 3- rho
      if (nrv3m.gt.0) then
** Eq.(12) in EPJA16p537(03) [hep-ph/0211457]
         do i1=1,3
            do i2=1,3
               do i3=1,3
                  t3jx3(i1,i2,i3)=t1jx3(i1)*t1jx3(i2)*t1jx3(i3)
     &                 -r2jx3/5.0*
     $                 ( del(i1,i2)*t1jx3(i3)
     &                 + del(i2,i3)*t1jx3(i1)
     &                 + del(i3,i1)*t1jx3(i2) )
                  t3xx3(i1,i2,i3)=t1xx3(i1)*t1xx3(i2)*t1xx3(i3)
     &                 -r2xx3/5.0*
     $                 ( delx3(i1,i2)*t1xx3(i3)
     &                 + delx3(i2,i3)*t1xx3(i1)
     &                 + delx3(i3,i1)*t1xx3(i2) )
               enddo
            enddo
         enddo
** Eq.(16) in EPJA16p537(03) [hep-ph/0211457]
         b3jx3=bwbar3(q2jx3,fud1)*bwbar3(q2x3,fud2)
** Eq.(29) in EPJA16p537(03) [hep-ph/0211457]
         do i1=1,2
            temp3=0.0
            do i2=1,3
               if (i2.ne.i1) then
                  do i3=1,3
                     tmp=eps(i1,i2,i3,4)
                     if (tmp.ne.0.0) then
                        do i4=1,3
                           do i5=1,3
                              temp3=temp3+tmp*amj*t3jx3(i2,i4,i5)
     $                             *t3xx3(i3,i4,i5)
                           enddo
                        enddo
                     endif
                  enddo
               endif
            enddo
            i4=2*nkvxx+nrv1m+1
            i5=2*nkvxx+nrvxx
            do i2=i4,i5
               IF(I2.NE.27)  THEN
                  if (kon(i2).eq.1) then
                     if (kff(i2).eq.1) then
                        cf(i2,i1)=temp3*b3jx3*cbwp1(i2,sx3)
                     elseif (kff(i2).eq.2) then
                        cf(i2,i1)=temp3*b3jx3*cbwp2(i2,3,sx3,sp1,sp2)
                     endif
                  endif
               ENDIF
            enddo
** for background of Ks K+- 
            if ( kon(27).eq.1 ) then
               cf(27,i1)=temp3*b3jx3 + cf(27,i1)
            endif
         enddo
      endif
c---- 2+ a2
      if (nav2p.gt.0) then
** Eq.(11) in EPJA16p537(03) [hep-ph/0211457]
         do i1=1,3
            do i2=1,3
               t2jx3(i1,i2)=t1jx3(i1)*t1jx3(i2)-r2jx3*del(i1,i2)/3.0
               t2xx3(i1,i2)=t1xx3(i1)*t1xx3(i2)-r2xx3*delx3(i1,i2)/3.0
            enddo
         enddo
** Eq.(15) in EPJA16p537(03) [hep-ph/0211457]
         b2jx3 = bwbar2(q2jx3,fud1)*bwbar2(q2x3,fud2)
         do i1=1,2
            temp3=0.0
            do i2=1,3
               if (i2.ne.i1) then
                  do i3=1,3
                     tmp=eps(i1,i2,i3,4)
                     if (tmp.ne.0.0) then
                        do i4=1,3
                           temp3=temp3-tmp*amj*t2jx3(i2,i4)*t2xx3(i3,i4)
                        enddo
                     endif
                  enddo
               endif
            enddo
            i4=2*nkvxx+nrvxx+1
            i5=2*nkvxx+nrvxx+nav2p
            do i2=i4,i5
              IF(I2.NE.30) THEN
                if (kon(i2).eq.1) then
                  if (kff(i2).eq.1) then
                    cf(i2,i1)=temp3*b2jx3*cbwp1(i2,sx3)
                  elseif (kff(i2).eq.2) then
                    cf(i2,i1)=temp3*b2jx3*cbwp2(i2,2,sx3,sp1,sp2)
                  endif
                endif
              ENDIF
            enddo
** for background of Ks K+- 
            if ( kon(30).eq.1 ) then
               cf(30,i1)=temp3*b2jx3
            endif
         enddo
      endif
c---- 4+ a
      if (nav4p.gt.0) then
** Eq.(37) in PRD48p1229(93)
         do i1=1,3
            do i2=1,3
               do i3=1,3
                  do i4=1,3
                     t4jx3(i1,i2,i3,i4)=t1jx3(i1)*t1jx3(i2)*t1jx3(i3)
     $                    *t1jx3(i4)-r2jx3/7.0*
     $                    ( del(i1,i2)*t1jx3(i3)*t1jx1(i4)
     &                    + del(i1,i3)*t1jx3(i2)*t1jx1(i4)
     &                    + del(i1,i4)*t1jx3(i2)*t1jx1(i3)
     &                    + del(i2,i3)*t1jx3(i1)*t1jx1(i4)
     &                    + del(i2,i4)*t1jx3(i1)*t1jx1(i3)
     &                    + del(i3,i4)*t1jx3(i1)*t1jx1(i2) )
     &                    +(r2jx3**2)/35.0*
     $                    ( del(i1,i2)*del(i3,i4)
     &                    + del(i1,i3)*del(i2,i4)
     &                    + del(i1,i4)*del(i2,i3) )
                     t4xx3(i1,i2,i3,i4)=t1xx3(i1)*t1xx1(i2)
     $                    *t1xx3(i3)*t1xx3(i4)-r2xx3/7.0*
     $                    ( delx3(i1,i2)*t1xx3(i3)*t1xx3(i4)
     &                    + delx3(i1,i3)*t1xx3(i2)*t1xx3(i4)
     &                    + delx3(i1,i4)*t1xx3(i2)*t1xx3(i3)
     &                    + delx3(i2,i3)*t1xx3(i1)*t1xx3(i4)
     &                    + delx3(i2,i4)*t1xx3(i1)*t1xx3(i3)
     &                    + delx3(i3,i4)*t1xx3(i1)*t1xx3(i2) )
     &                    +(r2xx3**2)/35.0*
     $                    ( delx3(i1,i2)*delx3(i3,i4)
     &                    + delx3(i1,i3)*delx3(i2,i4)
     &                    + delx3(i1,i4)*delx3(i2,i3) )
                  enddo
               enddo
            enddo
         enddo
** Eq.(17) in EPJA16p537(03) [hep-ph/0211457]
         b4jx3 = bwbar4(q2jx3,fud1)*bwbar4(q2x3,fud2)
         do i1=1,2
            temp4=0.0
            do i2=1,3
               if (i1.ne.i2) then
                  do i3=1,3
                     tmp=eps(i1,i2,i3,4)
                     if (tmp.ne.0.0) then
                        do i4=1,3
                           do i5=1,3
                              do i6=1,3
                                 temp4=temp4-tmp*amj*t4jx3(i2,i4,i5,i6)
     $                                *t4xx3(i3,i4,i5,i6)
                              enddo
                           enddo
                        enddo
                     endif
                  enddo
               endif
            enddo
            i4=2*nkvxx+nrvxx+nav2p+1
            i5=nresx
            do i2=i4,i5
               if (kon(i2).eq.1) then
                  if (kff(i2).eq.1) then
                     cf(i2,i1)=temp4*b4jx3*cbwp1(i2,sx3)
                  elseif (kff(i2).eq.2) then
                     cf(i2,i1)=temp4*b4jx3*cbwp2(i2,4,sx3,sp1,sp2)
                  endif
               endif
            enddo
         enddo
      endif

c---- background function
      if (nbgxx.gt.0) then
         nnnn=nbgxx
         do i2=1,nnnn
            i3=nresx+i2
            if (kon(i3).eq.1) then
               cf(i3,1)=cbgxx(i2,sx1)
               cf(i3,2)=0.0
            endif
         enddo
      endif
c---- amplitude of each partial wave, and their interference
      do i1=1,nmode
         do i2=1,nmode
            ampx(i1,i2)=0.0
         enddo
      enddo
      wa=0.0
      ii=x(1)
      do i1=1,nmode
         if (kon(i1).eq.1) then
            do i2=1,nmode
               if (kon(i2).eq.1) then
                  cw=cmplx(0.0,0.0)
                  do i3=1,2
                     cw=cw+cf(i1,i3)*conj(cf(i2,i3))/2.0
                  enddo
                  if (i1.le.i2) fu(i1,i2)=+dreal(cw)
                  if (i1.gt.i2) fu(i1,i2)=-dimag(cw)
                  ampx(i1,i2)=pa(i1,i2)*fu(i1,i2)
                  wa=wa+ampx(i1,i2)
               endif
            enddo
         endif
      enddo
c--- the total cross section must be a positive number
      if (wa.gt.0.0) goto 200
 100  continue
      do i1=1,nmode
         do i2=1,nmode
            ampx(i1,i2)=0.0
         enddo
      enddo
      wa=0.0
 200  continue
      fu(mmode,mmode)=spxx(dble(1.0))
      dcs=wa
      return
      end
*****&====== (*^_^*) (*^_^*)

*****&****************************************************************
***   Blatt-Weisskopf centrifugal barrier factor
***   refer: Eq.(13)---Eq.(17) in hep-ph/0211457
*****&================================================================
**  Eq.(13) in EPJA16p537(03) [hep-ph/0211457]
      function q2abc(sa,sb,sc)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
c      real*8 q2abc,sa,sb,sc
      q2abc = ((sa+sb-sc)**2)/(4.0*sa)-sb
      return
      end
**  Eq.(14) in EPJA16p537(03) [hep-ph/0211457]
      function bwbar1(q2,qz2)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      real*8 bwbar1,q2,qz2
      bwbar1 = sqrt( 2.0/(q2+qz2) )
      return
      end
**  Eq.(15) in EPJA16p537(03) [hep-ph/0211457]
      function bwbar2(q2,qz2)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      real*8 bwbar2,q2,qz2
      bwbar2 = sqrt( 13.0/(q2**2+3.0*q2*qz2+9.0*(qz2**2)) )
      return
      end
**  Eq.(16) in EPJA16p537(03) [hep-ph/0211457]
      function bwbar3(q2,qz2)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      real*8 bwbar3,q2,qz2
      bwbar3 = sqrt( 277.0/(q2**3+6.0*(q2**2)*qz2+45.0*q2*(qz2**2)
     &       + 225.0*(qz2**3)) )
      return
      end
**  Eq.(17) in EPJA16p537(03) [hep-ph/0211457]
      function bwbar4(q2,qz2)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      real*8 bwbar4,q2,qz2
      bwbar4 = sqrt( 12746.0/(q2**4+10.0*(q2**3)*qz2
     &       + 135.0*(q2**2)*(qz2**2)+1575.0*q2*(qz2**3)
     &       + 11025.0*(qz2**4)) )
      return
      end
*****&====== (*^_^*) (*^_^*)

*****&================================================================
      function funct(x)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      common/check/ak1(4),ak2(4),ak3(4),prob
      common/data/akn(nd1,3,5)
      common/a/faa(npa)
      common/cppp/cp(nmode),waa(npa),was
      common/pan/pa(nmode,nmode)
      common/cfnew/cf(nmode,2),kon(mmode),kon2
      common/trans/fu(mmode,mmode)
      common/coef/fun(mmode,mmode)
      common/gggg/fc
      dimension x(3)
*****&================================================================
      ii=x(1)
      if (ii.gt.1) goto 60
***   calculate the coefficients of partial wave, and the normalized
***   factor (i.e. the total cross section from MC Sample)
      call conv(cp)
      was=0.0
      do i1=1,nmode
      if (kon(i1).eq.1) then
        do i2=1,nmode
        if (kon(i2).eq.1) then
          cw=cp(i1)*conj(cp(i2))
          if (i1.eq.i2) pa(i1,i2)=cw
          if (i1.lt.i2) pa(i1,i2)=2.0*dreal(cw)
          if (i1.gt.i2) pa(i1,i2)=2.0*dimag(cw)
          was=was+pa(i1,i2)*fun(i1,i2)
        endif
        enddo
      endif
      enddo
      fc=faa(npa)   ! this quantity is used by function "grad"
      if (kon(mmode).eq.1) was=was+faa(npa)*fun(mmode,mmode)
      was=was/nmc
 60   continue
***   get four momentum of final state at the i-th experimental point
      do i1=1,4
        ak1(i1)=akn(ii,1,i1)
        ak2(i1)=akn(ii,2,i1)
        ak3(i1)=akn(ii,3,i1)
      enddo
***   calculate cross section at the i-th experimental point,
**    i.e. numerator of objective function
      temp1=dcs(pa)
      if (kon(mmode).eq.1) then
        temp2=faa(npa)*fu(mmode,mmode)
      else
        temp2=0.0
      endif
      prob =(temp1+temp2)
      funct=prob/was
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)

*****&================================================================
      subroutine eval
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      include './input/dat/res.dat'
      common/cfnew/cf(nmode,2),kon(mmode),kon2
      common/cppp/cp(nmode),waa(npa),was
      common/coef/fun(mmode,mmode)
      common/pan/pa(nmode,nmode)
*****&================================================================
      do i1=1,nmode
      if (kon(i1).eq.1) then
        do i2=1,nmode
        if (kon(i2).eq.1) then
          cw=cp(i1)*conj(cp(i2))
          if (i1.eq.i2) pa(i1,i2)=cw
          if (i1.lt.i2) pa(i1,i2)=2.0*dreal(cw)
          if (i1.gt.i2) pa(i1,i2)=2.0*dimag(cw)
        endif
        enddo
      endif
      enddo
      call integ(3)
      return
      end
*****&====== (*^_^*) (*^_^*)

*****&================================================================
      subroutine integ(iph)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      include './input/inc/const.in'
      common/check/ak1(4),ak2(4),ak3(4),prob
      common/a/faa(npa)
      common/cfnew/cf(nmode,2),kon(mmode),kon2
      common/trans/fu(mmode,mmode)
      common/coef/fun(mmode,mmode)
      common/pan/pa(nmode,nmode)
      common/ppwa/ampx(nmode,nmode)
*****&===============================================================
      if (iph.eq.1) then
        do i1=1,nmode
        do i2=1,nmode
          pa(i1,i2)=0.0
        enddo
        enddo
        do i1=1,mmode
        do i2=1,mmode
          fun(i1,i2)=0.0
        enddo
        enddo
c---- read four momentum of final states
        open(1,file='./data/PWA_mc.dat',status='old')
        do i1=1,nmc
          read(1,*) (ak1(i2),i2=1,4)
          read(1,*) (ak2(i2),i2=1,4)
          read(1,*) (ak3(i2),i2=1,4)
c---- calculate the integrated partial wave amplitudes
          temp=dcs(pa)
          do i2=1,mmode
          if (kon(i2).eq.1) then
            do i3=1,mmode
            if (kon(i3).eq.1) then
              fun(i2,i3)=fun(i2,i3)+fu(i2,i3)
            endif
            enddo
          endif
          enddo
        enddo
        close(1)
        return
      endif
*****&================================================================
c---  output of weight of MC Sample
      if (iph.eq.3) then
        open(10,file='./output/out/reweight.wt',status='unknown')
        open(11,file='./output/out/mc_int.wt',status='unknown')
        open(20,file='./data/PWA_mc.dat',
     & status='old')
        do i1=1,nmc
          read(20,*) (ak1(i2),i2=1,4)
          read(20,*) (ak2(i2),i2=1,4)
          read(20,*) (ak3(i2),i2=1,4)
          if (i1.eq.1) then
            was=0.0
            do i3=1,nmode
            if (kon(i3).eq.1) then
              do i4=1,nmode
              if (kon(i4).eq.1) then
                was=was+pa(i3,i4)*fun(i3,i4)
              endif
              enddo
            endif
            enddo
            was=was+faa(npa)*fun(mmode,mmode)
            wsf=dble(nsig)/was
          endif
          temp=wsf*dcs(pa)
          wabg=wsf*faa(npa)*fu(mmode,mmode)
          wate=temp+wabg
          write(10,*) wate,wabg

          do i3=1,nmode
          if (kon(i3).eq.1) then
            do i4=1,nmode
            if (kon(i4).eq.1) then
          if (i3.eq.i4) write(11,*)wsf*ampx(i3,i4)
          if (i3.lt.i4) write(11,*)wsf*(ampx(i3,i4)+ampx(i4,i3))
            endif
            enddo
          endif
          enddo
        enddo
        close(20)
        close(10)
      endif
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)


*****&================================================================
***   called by fumili for every event to form gradients
      subroutine arithm(y)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      common/x/x(3)
      common/df/df(npa)
      common/endflg/endflg,na,indflg(5)
      common/a/a(npa)
      common/pl/pl(npa)
      common/au/amx(npa)
      common/al/amn(npa)
      data rp/1.e-14/
*****&================================================================
***   unperturbed y
      y=funct(x)
***   run through all variables
      do i1=1,na
        df(i1)=0.0
        if (pl(i1).gt.0.0) then
          ai=a(i1)
***   hi=perturbation
          hi=0.01*pl(i1)
          pi=rp*abs(ai)
          if (hi.le.pi) hi=pi
          a(i1)=ai+hi
***   form gradient and restore parameters
          df(i1)=grad(x,hi,i1)
          a(i1)=ai
        endif
      enddo
      return
      end
*****&====== (*^_^*) (*^_^*)
*****&================================================================
***   form gradients (called by subroutine arithm)
***   note : this subroutine must not disturb cp or pa
c---  much of this subroutine is a carbon copy of funt
c---  ider=index of parameter in faa being varied
      function grad(x,step,ider)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      common/a/faa(npa)
      common/cfnew/cf(nmode,2),kon(mmode),kon2
      common/data/akn(nd1,3,5)
      common/check/ak1(4),ak2(4),ak3(4),prob
      common/trans/fu(mmode,mmode)
      common/coef/fun(mmode,mmode)
      common/norm/wa1
      common/cppp/cpn(nmode),waa(npa),was
      common/gggg/fc
      common/pan/pa(nmode,nmode)
      common/pad/pag(nmode,nmode)
      dimension cp(nmode),cdr(nmode),x(3)
      ii=x(1)
      call conv(cp)
      do i1=1,nmode
        if (kon(i1).eq.1) then
          cdr(i1)=(cp(i1)-cpn(i1))/step
        else
          cdr(i1)=cmplx(0.0,0.0)
        endif
      enddo
      ff=(faa(npa)-fc)/step
      dprob=0.0
      do i1=1,nmode
      if (kon(i1).eq.1) then
        do i2=1,nmode
        if (kon(i2).eq.1) then
          cw=cdr(i1)*conj(cpn(i2))+cpn(i1)*conj(cdr(i2))
          if (i1.eq.i2) pag(i1,i2)=cw
          if (i1.lt.i2) pag(i1,i2)=2.0*dreal(cw)
          if (i1.gt.i2) pag(i1,i2)=2.0*dimag(cw)
          dprob = dprob+pag(i1,i2)*fu(i1,i2)
        endif
        enddo
      endif
      enddo
      dprob=dprob+ff*fu(mmode,mmode)
      if (ii.ge.2)  goto 60
***   form perturbed denominator at first entry and store in waa
      wa1=0.0
      do i1=1,nmode
      if (kon(i1).eq.1) then
        do i2=1,nmode
        if (kon(i2).eq.1) then
          wa1=wa1+pag(i1,i2)*fun(i1,i2)
        endif
        enddo
      endif
      enddo
      waa(ider)=(wa1+ff*fun(mmode,mmode))/nmc
***   form gradient
 60   continue
      grad=dprob/was-prob/(was**2)*waa(ider)
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)


*****&****************************************************************
**  Breit-Wigner propagator
*****&================================================================
***** B.Zou and D.Bugg, Eur.Phys.J.A16,537(2003)
      function cbwp1(ipx,ss)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      include './input/dat/res.dat'
      include './input/inc/state.in'
      integer ipx
      complex*16 cbwp1,ci
      real*8 ss,xm0,xt0
      ci=cmplx(0.0,1.0)
      xm0=xmass(ipx)
      xt0=xwidth(ipx)
      cbwp1=1.0/(ss-xm0**2+ci*xm0*xt0)
      return
      end
***** J.H.Kuhn and A.Santamaria, Z.Phys.C48,445(1990)
      function cbwp2(ipx,jjx,ss,sm1,sm2)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      include './input/dat/res.dat'
      include './input/inc/state.in'
      integer ipx,jjx,ikx
      real*8 ss,xm0,xt0,sm0,sm1,sm2,sps,spx
      complex*16 cbwp2,ci,cqs,cqx,cts
      ikx=2*jjx+1
      ci=cmplx(0.0,1.0)
      xm0=xmass(ipx)
      xt0=xwidth(ipx)
      sm0=xm0**2
      sps=q2abc(ss, sm1,sm2)
      spx=q2abc(sm0,sm1,sm2)
      cqs=sqrt( abs(sps) )
      cqx=sqrt( abs(spx) )
      if (sps.lt.0.0) cqs=ci*cqs
      if (spx.lt.0.0) cqx=ci*cqx
      cts=xt0*(sm0/ss)*((cqs/cqx)**ikx)
      cbwp2=1.0/(ss-sm0+ci*sqrt(ss)*cts)
      return
      end
***** background
      include './input/bgform/bgxx.f'
      include './input/bgform/spxx.f'
*****&================================================================

*****&================================================================
**  scalar production of two vector in 4D
      function scalar(a,b)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      real*8  scalar,a(4),b(4)
      scalar = a(4)*b(4)-a(1)*b(1)-a(2)*b(2)-a(3)*b(3)
      return
      end
*****&====== (*^_^*) (*^_^*)

*****&================================================================
**  form complex*16 conjugate
      function conj(cx)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      complex*16 cx,conj
      conj = cmplx(dreal(cx),-dimag(cx))
      return
      end
*****&====== (*^_^*) (*^_^*)


*****&================================================================
** set up delta(i,j) and epsilon(i,j,k,l)
      subroutine metric
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      include'./input/inc/const.in'
      character cht*30
c---- unit of imaginary number
      ci=cmplx(0.0,1.0)
c---- metric tensor
      do i=1,4
      do j=1,4
        del(i,j)=0.0
        if (i.eq.j) del(i,j)=-1.0
        do k=1,4
        do l=1,4
          eps(i,j,k,l)=0.0
        enddo
        enddo
      enddo
      enddo
      del(4,4)=1.0
c---- 4-rank asymmetric tensor
      eps(1,2,3,4) =  1.0
      eps(1,2,4,3) = -1.0
      eps(1,3,2,4) = -1.0
      eps(1,3,4,2) =  1.0
      eps(1,4,2,3) =  1.0
      eps(1,4,3,2) = -1.0
      eps(2,1,3,4) = -1.0
      eps(2,1,4,3) =  1.0
      eps(2,3,1,4) =  1.0
      eps(2,3,4,1) = -1.0
      eps(2,4,1,3) = -1.0
      eps(2,4,3,1) =  1.0
      eps(3,1,2,4) =  1.0
      eps(3,1,4,2) = -1.0
      eps(3,2,1,4) = -1.0
      eps(3,2,4,1) =  1.0
      eps(3,4,1,2) =  1.0
      eps(3,4,2,1) = -1.0
      eps(4,1,2,3) = -1.0
      eps(4,1,3,2) =  1.0
      eps(4,2,1,3) =  1.0
      eps(4,2,3,1) = -1.0
      eps(4,3,1,2) = -1.0
      eps(4,3,2,1) =  1.0
c---- mass of J/psi, K+, K0, pi+, pi0
      amj  = 4.680    ! mass of J/psi
      amkp = 1.86483  ! mass of K+
      amkz = 1.86965 ! mass of K0
      ampip= 0.139570 ! mass of pi+
      ampiz= 0.135    ! mass of pi0
      sj   = amj**2
      skp  = amkp**2
      skz  = amkz**2
      spip = ampip**2
      spiz = ampiz**2
c--- centrifugal barrier factors
      open(1,file='./input/dat/bar.dat',status='old')
      read(1,10) cht, fud1
      read(1,10) cht, fud2
      read(1,10) cht, fud3
      close(1)
 10   format(4x,a30,f15.10)
      return
      end
*****&====== (*^_^*) (*^_^*)


*****&================================================================
      subroutine fumili(s,m,n1,n2,n3,eps,akappa,alambd,it,mc)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      call fumixx(s,m,n1,n2,n3,eps,akappa,alambd,it,mc,1)
      return
      end
      subroutine likelm(s,m,n1,n2,n3,eps,akappa,alambd,it,mc)
      implicit double precision (a,b,d-h,o-z)
      implicit complex*16 (c)
      call fumixx(s,m,n1,n2,n3,eps,akappa,alambd,it,mc,2)
      return
      end
*****&================================================================
      subroutine fumixx(s,m,n1,n2,n3,eps,akappa,alambd,it,mc,nlog)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      common/z/z(nzz)
      common/g/g(npa)
      common/a/a(npa)
      common/pl/pl0(npa)
      common/sigma/sigma(npa)
      common/plu/pl(npa)
      common/r/r(npa)
      common/da/da(npa)
      common/endflg/endflg,na,indflg(5)
      common/au/amx(npa)
      common/al/amn(npa)
      common/z0/z0(nzz)
      data rp/1.e-14/
      goto (201,202),nlog
 201  continue
      indflg(3)=0
 1    if (it.ge.0) write(*,84)
      call vfll(amx,amn,npa,1.D35)
      nn2=0
      n=m
      fixflg=0.0
      endflg=0.0
      indflg(2)=0
      ifix1=0.
      fi=0.0
      nn3=0
      do 2 i=1,n
      r(i)=0.0
      if (eps.gt.0.0) sigma(i)=0.0
 2    pl(i)=pl0(i)
c-----start new iteration
 3    nn1=1
      t1=1.0
c-----repeat iteration with smaller step
 4    s=0.0
      n0=0
      do 7 i=1,n
      g(i)=0.0
      if (pl0(i)) 7,7,5
 5    n0=n0+1
      if (pl(i)) 7,7,6
 6    pl0(i)=pl(i)
 7    continue
      nn0=n0*(n0+1)/2
      if (nn0.lt.1) goto 9
      do 8 i=1,nn0
      z(i)=0.0
 8    continue
 9    na=m
      indflg(1)=0
c-----calculate objective function
      call sgz (m,s)
      sp=rp*abs(s)
      if (nn0.lt.1) goto 11
      do 10 i=1,nn0
      z0(i)=z(i)
 10   continue
 11   if (nn3) 19,19,12
 12   if (nn1-n1) 13,13,19
 13   t=2.0*(s-olds-gt)
      if (indflg(1)) 16,14,16
 14   if (abs(s-olds).le.sp.and.-gt.le.sp) goto 19
      if (0.59*t+gt) 19,15,15
 15   t=-gt/t
      if (t-0.25) 16,17,17
 16   t=0.25
 17   gt=gt*t
      t1=t1*t
      nn2=0
      do 18 i=1,n
      if (pl(i).le.0.0) go to 18
      a(i)=a(i)-da(i)
      pl(i)=pl(i)*t
      da(i)=da(i)*t
      a(i)=a(i)+da(i)
 18   continue
      nn1=nn1+1
      goto 4
c-----remove contribution of fixed parameters from z
 19   if (indflg(1).eq.0) goto 20
      endflg=-4.0
      goto 85
 20   k1=1
      k2=1
      i1=1
      do 30 i=1,n
      if (pl0(i)) 30,30,21
 21   if (pl(i).eq.0.0) pl(i)=pl0(i)
      if (pl(i)) 23,23,24
 22   pl(i)=0.0
 23   k1=k1+i1
      go to 29
 24   if (a(i).ge.amx(i).and.g(i).lt.0.0) go to 22
      if (a(i).le.amn(i).and.g(i).gt.0.0) go to 22
      do 28 j=1,i
      if (pl0(j)) 28,28,25
 25   if (pl(j)) 27,27,26
 26   z(k2)=z0(k1)
      k2=k2+1
 27   k1=k1+1
 28   continue
 29   i1=i1+1
 30   continue
c-----invert z
      i1=1
      l=i1
      do 32 i=1,n
      if (pl(i)) 32,32,31
 31   r(i)=z(l)
      i1=i1+1
      l=l+i1
 32   continue
      n0=i1-1
      call mconv (n0)
      if (indflg(1)) 33,34,33
 33   indflg(1)=0
      indflg(2)=1
      go to 49
 34   continue
c-----calculate theoretical step to minimum
      i1=1
      do 41 i=1,n
      da(i)=0.0
      if (pl(i)) 41,41,35
 35   l1=1
      do 40 l=1,n
      if (pl(l)) 40,40,36
 36   if (i1-l1) 37,37,38
 37   k=l1*(l1-1)/2+i1
      go to 39
 38   k=i1*(i1-1)/2+l1
 39   da(i)=da(i)-g(l)*z(k)
      l1=l1+1
 40   continue
      i1=i1+1
 41   continue
c-----check for parameters on boundary
      afix=0.0
      ifix=0
      i1=1
      l=i1
      do 47 i=1,n
      if (pl(i)) 47,47,42
 42   sigi=sqrt(abs(z(l)))
      r(i)=r(i)*z(l)
      if (eps) 44,44,43
 43   sigma(i)=sigi
 44   if   ((a(i).lt.amx(i).or.da(i).le.0.0)
     & .and.(a(i).gt.amn(i).or.da(i).ge.0.0)) goto 46
      akap=abs(da(i)/sigi)
      if (akap-afix) 46,46,45
 45   afix=akap
      ifix=i
      ifix1=i
 46   i1=i1+1
      l=l+i1
 47   continue
      if (ifix) 48,50,48
 48   pl(ifix)=-1.0
 49   fixflg=fixflg+1.0
      fi=0.0
c-----repeat calculation of theoretical step after fixing each parameter
      go to 19
c-----calculate step correction factor
 50   alambd=1.0
      akappa=0.0
      imax=0
      do 60 i=1,n
      if (pl(i)) 60,60,51
 51   bm=amx(i)-a(i)
      abi=a(i)+pl(i)
      abm=amx(i)
      if (da(i)) 52,52,53
 52   bm=a(i)-amn(i)
      abi=a(i)-pl(i)
      abm=amn(i)
 53   bi=pl(i)
      if (bi-bm) 55,55,54
 54   bi=bm
      abi=abm
 55   if (abs(da(i))-bi) 58,58,56
 56   al=abs(bi/da(i))
      if (alambd-al) 58,58,57
 57   imax=i
      aimax=abi
      alambd=al
 58   akap=abs(da(i)/sigma(i))
      if (akap-akappa) 60,60,59
 59   akappa=akap
 60   continue
c-----calculate new corrected step
      gt=0.0
      amb=1.e18
      if (alambd) 62,62,61
 61   amb=0.25/alambd
 62   continue
      do 67 i=1,n
      if (pl(i)) 67,67,63
 63   if (nn2-n2) 66,66,64
 64   if (abs(da(i)/pl(i))-amb) 66,65,65
 65   pl(i)=4.0*pl(i)
      t1=4.0
 66   da(i)=da(i)*alambd
      gt=gt+da(i)*g(i)
 67   continue
c-----check if minimum attained and set exit mode
      if (-gt.gt.sp.or.t1.ge.1.0.or.alambd.ge.1.0) goto 68
      endflg=-1.0
 68   if (endflg) 85,69,69
 69   if (akappa-abs(eps)) 70,75,75
 70   if (fixflg) 72,71,72
 71   endflg=1.0
      go to 85
 72   if (endflg) 85,77,73
 73   if (ifix1) 85,85,76
 74   if (fi-fixflg) 76,76,77
 75   if (fixflg) 74,76,74
 76   fi=fi+1.0
      endflg=0.0
 85   if(endflg.eq.0.0.and.nn3.ge.n3) endflg=-3.0
      if(endflg.gt.0.0.and.indflg(2).gt.0) endflg=-2.0
      call monito (s,m,nn3,it,eps,gt,akappa,alambd)
      if (endflg) 83,79,83
c-----check if fixing on bound is correct
 77   endflg=1.0
      fixflg=0.0
      ifix1=0
      do 78 i=1,m
 78   pl(i)=pl0(i)
      indflg(2)=0
      go to 19
c-----next iteration
 79   endflg=0.0
      do 80 i=1,n
      a(i)=a(i)+da(i)
 80   continue
      if (imax) 82,82,81
 81   a(imax)=aimax
 82   olds=s
      nn2=nn2+1
      nn3=nn3+1
      go to 3
 83   mc=endflg
      return
 202  continue
      indflg(3)=1
      goto 1
 84   format(/3x,61('*'),/3x,'*',59x,'*',/3x,'*',
     &  5x,'function minimisation by subroutine fumili/likelm',5x,'*',
     & /3x,'*',17x,'in the following print-out',16x,'*',
     & /3x,'*',59x,'*',
     & /3x,'*   s = value of objective function',25x,'*',
     & /3x,'*  ec = expected change in s during next iteration',10x,'*',
     & /3x,'*  kappa  = estimated distance to minimum',19x,'*',
     & /3x,'*  lambda = step length modifier',28x,'*',
     & /3x,'*',59x,'*',/3x,61('*'),/)
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)

*****&================================================================
***  inverts the positive definite packed symmetric matrix z by the
***  square-root method
      subroutine mconv(n)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      common/z/z(nzz)
      common/plu/pl0(npa)
      common/r/r(npa)
      common/endflg/endflg,na,indflg(5)
c-----maximum real*8 number and 10.*maximum relative precision on cdc6000
      data am,rp/1.e35,1.e-14/
      if (n.lt.1) return
      aps=sqrt(am/float(n))
      ap=1.0/(aps*aps)
      ir=0
      do 11 i=1,n
 1    ir=ir+1
      if (pl0(ir)) 1,1,2
 2    ni=i*(i-1)/2
      ii=ni+i
      k=n+1
      if (z(ii).le.rp*abs(r(ir)).or.z(ii).le.ap) goto 19
      z(ii)=1./sqrt(z(ii))
      nl=ii-1
 3    if (nl-ni) 5,5,4
 4    z(nl)=z(nl)*z(ii)
      if (abs(z(nl)).ge.aps) goto 16
      nl=nl-1
      goto 3
 5    if (i-n) 6,12,12
 6    k=k-1
      nk=k*(k-1)/2
      nl=nk
      kk=nk+i
      d=z(kk)*z(ii)
      c=d*z(ii)
      l=k
 7    ll=nk+l
      li=nl+i
      z(ll)=z(ll)-z(li)*c
      l=l-1
      nl=nl-l
      if (l-i) 9,9,7
 8    ll=nk+l
      li=ni+l
      z(ll)=z(ll)-z(li)*d
 9    l=l-1
      if (l) 10,10,8
 10   z(kk)=-c
      if (k-i-1) 11,11,6
 11   continue
 12   do 14 i=1,n
      do 14 k=i,n
      nl=k*(k-1)/2
      ki=nl+i
      d=0.0
      do 13 l=k,n
      li=nl+i
      lk=nl+k
      d=d+z(li)*z(lk)
      nl=nl+l
 13   continue
      ki=k*(k-1)/2+i
      z(ki)=d
 14   continue
 15   return
 16   k=i+nl-ii
      ir=0
      do 18 i=1,k
 17   ir=ir+1
      if (pl0(ir)) 17,17,18
 18   continue
 19   pl0(ir)=-2.0
      r(ir)=0.0
      indflg(1)=ir
      go to 15
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)


*****&================================================================
***   sets up s (objective function),
***           g (gradient of s),
***           z (approximate covariance matrix)
      subroutine sgz (m,s)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      common/g/g(npa)
      common/z/z(nzz)
      common/a/a(npa)
      common/df/df(npa)
      common/ned/nt,ns
      common/exda/exda(nex)
      common/x/x(3)
      common/pl/pl(npa)
      common/endflg/endf6g,na,indflg(5)
      common/uprmon/ikk
      we=1.0
      k=nt
      k2=1
      do 12 l1=1,k
      k1=k2
      nx=ns-2
      if (indflg(3)) 1,2,1
 1    nx=ns
      k1=k1-2
 2    continue
      do 3 i=1,nx
      ki=k1+1+i
      x(i)=exda(ki)
 3    continue
      call arithm (y)
      if(ikk.eq.2) then
      write(*,100) x(1),y,exda(k1),exda(k1+1)
 100  format(5x,4e15.4)
      endif
      if (indflg(3)) 4,6,4
 4    if (y) 13,13,5
c-----maximum liklehood
 5    continue

c      if(l1.gt.ndt.and.l1.le.ntot)             we=-0.6
      if(l1.gt.ndt.and.l1.le.(ndt+nbg1))             we=-1.80
      if(l1.gt.(ndt+nbg1).and.l1.le.(ndt+nbg1+nbg2)) we=-0.00
      if(l1.gt.(ndt+nbg1+nbg2).and.l1.le.(ntot))     we=-0.00
**      if(l1.gt.Nt1.and.l1.le.nt)    we=-1.0
      s=s-log(y)*we
      y=-y
      sig=y
      go to 7
c-----chi squared
 6    sig=exda(k2+1)
      y=y-exda(k1)
      s=s+((y/sig)**2)/2.
 7    continue
      n=0
      do 9 j=1,m
      if (pl(j)) 9,9,8
 8    n=n+1
      df(n)=df(j)/sig
      g(j)=g(j)+df(n)*(y/sig)*we
 9    continue
      l=1
      if (n.lt.1) go to 11
      do 10 i=1,n
      do 10 j=1,i
      z(l)=z(l)+df(i)*df(j)*we
      l=l+1
 10   continue
 11   k2=k2+ns
 12   continue
      return
c-----  -ve or zero y in maximum liklehood
 13   indflg(1)=1
      s=1.e10
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)


*****&================================================================
      subroutine monito(s,m,nn3,it,eps,gt,akappa,alambd)
      implicit double precision (a,b,d-h,o-z)
      include './input/dat/num.dat'
      include './input/dat/res.dat'
      common/a/a(npa)
      common/sigma/sigma(npa)
      common/r/r(npa)
      common/pl/pl(npa)
      common/plu/pl0(npa)
      common/endflg/endflg,na,indflg(5)
      common/au/amx(npa)
      common/al/amn(npa)
      common/uprmon/ikk
      if (it)    11, 3,1
 1    if (nn3)    4, 4,2
 2    if (nm)     3, 4,4
 3    if (endflg) 4,12,4
c---  printer carriage control
 4    i1=6
      if (m.gt.6)  i1=5
      if (m.gt.12) i1=4
      if (m.gt.23) i1=1
c---  non-ansi carriage control suppressed on ibm
      if (i1.gt.1) i1=0
      print 19,i1,nn3,s,gt,akappa,alambd
      if (ikk.ge.1) write(*,191)
      do 10 i=1,m
        if (pl(i))  9,9,5
 5      if (pl0(i)) 8,7,6
 6      continue
        if(ikk.ge.1) write(*,20) i,a(i),sigma(i),r(i)
        goto 10
 7      continue
        if(ikk.ge.1) write(*,21) i,a(i),sigma(i),r(i)
        goto 10
 8      if (pl0(i).ge.-1.0) goto 7
        if(ikk.ge.1) write(*,22) i,a(i)
        goto 10
 9      continue
        if(ikk.ge.1) write(*,23) i,a(i)
 10   continue
 11   nm=-it
 12   nm=nm+1
      if (endflg) 13,14,14
 13   i=-endflg
      goto (15,16,17,18), i
 14   return
 15   write(*,24)
      goto  14
 16   write(*,25)
      goto  14
 17   write(*,26)
      goto  14
 18   write(*,27)
      goto  14
 19   format(/1x,i5,' **********',2x,'iteration no.',i3,
     &       /4x,'s= ',f12.3,3x,'ec = ',f10.3,3x,'kappa= ',f9.3,
     &        3x,'lambda= ',f8.3)
 191  format(/1x,6x,9hparameter,6x,9hparameter,9x,8hstandard,
     &        8x,11hcorrelation/1x,8x,6hnumber,9x,5hvalue,11x,
     &        9hdeviation,9x,6hfactor/)
 20   format(1x,8x,i3,4x,3(5x,e12.5))
 21   format(1x,8x,i3,4x,3(5x,e12.5),' parameter on boundary')
 22   format(1x,8x,i3,9x,e12.5,5x,' infinite error estimated')
 23   format(1x,8x,i3,9x,e12.5,5x,' this parameter fixed')
 24   format(1x,'minimisation terminated as no further decrease',
     &       1x,'in s is obtainable',/)
 25   format(1x,'minimisation terminated as infinite errors',
     &       1x,'estimated',/)
 26   format(1x,'minimisation terminated as iteration limit',
     &       1x,'reached',/)
 27   format(1x,'minimisation terminated as negative or zero',
     &       1x,'y encountered as logarithmic arguement',/)
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)

*****&================================================================
****  upper-limit and lower-limit of parameters
      subroutine vfll(amx,amn,n,edd)
      implicit complex*16 (c)
      implicit double precision (a,b,d-h,o-z)
      include './input/dat/res.dat'
      dimension amx(npa),amn(npa)
      common/pl/pl(npa)
      do i=1,n
        if (amx(i).lt.amn(i)) then
          write(*,100) i
          goto 10
        endif
        if (amx(i).eq.0.0.and.amn(i).eq.0.0) then
          amx(i)=+edd
          amn(i)=-edd
        endif
 10     continue
      enddo
c      do i=1,n
c        if (amx(i).eq.0.0) write (*,200) i
c        if (amn(i).eq.0.0) write (*,300) i
c      enddo
 100  format(5x,'maximum less than minimum for parametr',i3)
 200  format(5x,'warning : maximum is zero for parametr',i3)
 300  format(5x,'warning : minimum is zero for parametr',i3)
*****&================================================================
      return
      end
*****&====== (*^_^*) (*^_^*)
