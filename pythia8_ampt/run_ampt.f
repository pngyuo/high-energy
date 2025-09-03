      subroutine ini_ampt(pzA, pzB, KSEED)

      implicit double precision (a-h, o-z)

      real pzA, pzB

c     initialization value for parton cascade:
      common /para2/ xmp, xmu, alpha, rscut2, cutof2
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common /rndm3/ iseedp

cc      SAVE /SOFT/
      common/anim/nevent,isoft,isflag,izpc

      COMMON /AREVT/ IAEVT, IARUN, MISS

      real pslimit
c     parton coalescence radii in case of string melting:
      common/coal/dpcoal,drcoal,ecritl,xmbcut,pslimit,imbcut,drbmRatio

      real efrm,epsiPz,epsiPt,PZPROJ,PZTARG
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG

clin-9/2018 common:
      real rdstard,rbstarb,rtstart,rrhopi,romrho0,rkstark,pttrig
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
     1     ,rdstard,rbstarb,rtstart,rrhopi,romrho0,rkstark,ihfdcy

      real ARPAR1, ARINT1
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)

      !pert deuteron not used but needed to keep ART running
      common /para8/ idpert,npertd,idxsec

      real DT
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      common/oscar1/iap,izp,iat,izt
      common/oscar2/FRAME,amptvn
      common/resdcy/NSAV,iksdcy

      integer NSEED
      COMMON/RNDF77/NSEED

      character*25 amptvn


      SAVE 

      !NTMAX=2 !2 for negligible hadron cascade
      DT=0.2

      !ioscar=0 only zpc.dat, ampt.dat produced
      !ioscar>1 primordial_pythia.dat produced
      !ioscar=3 parton-collisionHistory.dat produced
      ioscar=3

      !isoft=1 string hadronization, only gluon cascade
      !isoft=4 string melting, coalescence
      !isoft=2 string hadronization, breakup diquark/resonance with beam
      !remnant, q+g parton cascade, only valid at parton level at this
      !moment
      !isoft=4

      !used for internal consistancy check
      PZPROJ=abs(real(pzA))
      PZTARG=abs(real(pzB))
      epsiPz=MAX(PZPROJ*0.01,0.5)
      epsiPt=MAX(PZPROJ*0.001,0.1)
      iseedp=8
      print*,'PZROJ=',PZROJ,' epsiPz=',epsiPz,' pzB=',pzB

      iksdcy=0
      iphidcy=1
      ipi0dcy=0

      !dum par initialized for deuteron routines only to switch it off
      idpert=0
      npertd=1
      idxsec=1

      !xmu=2.265d0 !3mb
      !xmu=3.203d0 !1.5mb
      !xmu=99997.162d0 !0.1mb
      !alpha=0.33d0

      !scattering angle distribution in zpc
      !D=0: forward angle parton scattering, 100: isotropic distribution
      izpc=0

      NSEED=MOD((2*KSEED+1), 1000000)
      CALL SRAND(NSEED)

      print*,'HIJING random number seed:',NSEED
      print*,'mu=',xmu,' alpha=',alpha,' isoft=',isoft,' ntmax=',ntmax

      imbcut=0
      xmbcut=0d0
      !drbmRatio=0.72 ! in test value to be used with pythia8
      !drbmRatio=0.61d0 !lin&he new coal tune
      drbmRatio=0.53d0 !new value from chao's fit to data 20181207
      !drbmRatio=0.40d0 !in test value with pythia8 no FSI
clin-8/2018:
      OPEN (97, FILE = 'ana/primordial-hf.dat', STATUS = 'UNKNOWN')
c     1) This sets the ratio of (all primordial D*)/(all primordial D)
c     ~= primordial-D*0/primordial-D0; data suggest value within (0.96,3):
      !rdstard=1.0 !default value
      rdstard=2.0
c     2) similar for the B*/B ratio & T*/T ratio:
c     Note that thermal model constrains all these 6 ratios within (0, 3):
      rbstarb=2.0
      rtstart=2.0
c     3) Implement heavy resonance decays at last timestep when ihfdcy=1:
      ihfdcy=1
clin-9/2018 set rrhopi,romrho0,rkstark:
c     4) They set respectively ratios of (primordial rho0)/(primordial pi0),
c     (primordial omega)/(primordial eta),(primordial K*)/(primordial K),
c     data suggest value within (0.27,0.35), 1/1.4, (0.43,0.54) respectively:
      !rrhopi=0.30
      rrhopi=0.30 !in test value to be used with pythia8
      romrho0=0.90
      rkstark=0.50



c     AMPT momentum and space info at freezeout:
      OPEN (16, FILE = 'ana/ampt.dat', STATUS = 'UNKNOWN')
      open(14,file='ana/zpc.dat',status='unknown')

      !before tranport to ART
      OPEN (33, FILE = 'ana/beforeART.dat', STATUS = 'UNKNOWN')
      OPEN (32, FILE = 'ana/primordial_pythia.dat', STATUS = 'UNKNOWN')


      call ARTSET
      CALL INIZPC

      if(isoft.eq.1) then
         amptvn = '1.31t3 (Default)'
      elseif(isoft.eq.4) then
         amptvn = '2.31t3 (StringMelting)'
      else
         amptvn = 'Test-Only'
      endif
      WRITE(6,50) amptvn
 50   FORMAT(' '/
     &11X,'##################################################'/1X,
     &10X,'#      AMPT (A Multi-Phase Transport) model      #'/1X,
     &10X,'#          Version ',a25,                  '     #'/1X,
     &10X,'#               09/16/2018                       #'/1X,
     &10X,'##################################################'/1X,
     &10X,' ')

c     string formation time: to be used in ART
      ARPAR1(1) = 0.7



      return
      end


      subroutine run_ampt(ievt, ncoll, bval)

      implicit double precision (a-h, o-z)

      PARAMETER (MAXSTR=150001)
      PARAMETER (MAXPTN=400001)
      PARAMETER (MAXIDL=4001)

      real EATT, PATT
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)

      COMMON /PARA1/ MUL
cc      SAVE /PARA1/
      COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)  
      COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)

      real GXAR,GYAR,GZAR,FTAR,PXAR,PYAR,PZAR,PEAR,XMAR
       COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
cc      SAVE /ARPRC/

      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)

      common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1       xstrg0(MAXPTN),ystrg0(MAXPTN),
     2       xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)


      ! in zpc prc2 variables do not have 5 label 
      double precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
      integer ityp5
        common /prec2/gx5(MAXPTN),gy5(MAXPTN),gz5(MAXPTN),ft5(MAXPTN),
     &       px5(MAXPTN), py5(MAXPTN), pz5(MAXPTN), e5(MAXPTN),
     &       xmass5(MAXPTN), ityp5(MAXPTN)
cc     SAVE /prec2/

      real PXSG,PYSG,PZSG,PESG,PMSG
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)

c     7/20/01: use double precision
c     otherwise sometimes beta>1 and gamma diverge in lorenz():
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/

      real GXN,GYN,GZN,FTN,PXN,PYN,PZN,EEN,XMN
clin-4/26/01 lepton and photon info:
      COMMON /NOPREC/ NNOZPC, ITYPN(MAXIDL),
     &       GXN(MAXIDL), GYN(MAXIDL), GZN(MAXIDL), FTN(MAXIDL),
     &       PXN(MAXIDL), PYN(MAXIDL), PZN(MAXIDL), EEN(MAXIDL),
     &       XMN(MAXIDL)
cc      SAVE /NOPREC/

      !ART parameters
      real ARPAR1, ARINT1
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)

      common/anim/nevent,isoft,isflag,izpc

      common /para7/ ioscar,nsmbbbar,nsmmeson


      !icolln output number of collisions in zpc
      common /ilist6/ t, iopern, icolln
cc    SAVE /ilist6/

clin-6/22/01:
      real bimp
      common /lastt/itimeh,bimp
      COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP

      COMMON /AREVT/ IAEVT, IARUN, MISS


      IAEVT=ievt           
      NINTHJ=ncoll
      bimp=bval

      !print*,'NATT=', NATT

      !primordial pythia data output
      if(ioscar.gt.1) then
      write(32,*) ievt,' ',1,' ',NATT,bimp,' ',0,' ',0,' ',0,' ',0,' ',
     & 0,NINTHJ,0.0d0
      do i=1,NATT
         write(32,*) ITYPAR(i),' ',PXAR(i),' ',PYAR(i),' ',PZAR(i),' ',
     &  PEAR(i),' ',GXAR(i),' ',GYAR(i),' ',GZAR(i),' ',
     &  FTAR(i)
      enddo
      endif


      CALL HTOP

      nsp=0
      nst=0
      nsg=natt
      NSI=NSG

      CALL ZPCMN

      !record the number of collisions in zpc for output
      NELT=icolln

      !print*,'MUL=', MUL

c.....transfer data back from ZPC ready for coalescence
        DO 1011 I = 1, MAXSTR
           DO 1012 J = 1, 3
              K1SGS(I, J) = 0
              K2SGS(I, J) = 0
              PXSGS(I, J) = 0d0
              PYSGS(I, J) = 0d0
              PZSGS(I, J) = 0d0
              PESGS(I, J) = 0d0
              PMSGS(I, J) = 0d0
              GXSGS(I, J) = 0d0
              GYSGS(I, J) = 0d0
              GZSGS(I, J) = 0d0
              FTSGS(I, J) = 0d0
 1012      CONTINUE
 1011   CONTINUE

        NQP=0
        NQPB=0
        NGLU=0

c        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
        if(ioscar.ge.1) 
     &  WRITE (14, 210) IAEVT, MISS, MUL, bimp, 0, 0, 0, 0

        DO 1013 I = 1, MUL
           IITYP=ITYP5(I)
           NSTRG = LSTRG1(I)
           NPART = LPART1(I)
           K2SGS(NSTRG, NPART) = ITYP5(I)
           PXSGS(NSTRG, NPART) = PX5(I)
           PYSGS(NSTRG, NPART) = PY5(I)
           PZSGS(NSTRG, NPART) = PZ5(I)
           PMSGS(NSTRG, NPART) = XMASS5(I)
clin-7/20/01 E5(I) does no include the finite parton mass XMASS5(I), 
c     so define it anew:
c           PESGS(NSTRG, NPART) = E5(I)
c           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0) 
c     1          write(91,*) 'a',PX5(i),PY5(i),XMASS5(i),PZ5(i),E5(i)
           E5(I)=dsqrt(PX5(I)**2+PY5(I)**2+PZ5(I)**2+XMASS5(I)**2)
           PESGS(NSTRG, NPART) = E5(I)
c           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0) 
c     1          write(91,*) 'b: new E5(I)=',E5(i)
clin-7/20/01-end
           GXSGS(NSTRG, NPART) = GX5(I)
           GYSGS(NSTRG, NPART) = GY5(I)
           GZSGS(NSTRG, NPART) = GZ5(I)
           FTSGS(NSTRG, NPART) = FT5(I)

           if(IITYP.eq.21) then
             NGLU=NGLU+1
c             print*,'G px=',PX5(I),' py=',PY5(I),' pz=',PZ5(I),' rap=',
     &  -0.5*log((E5(I)-PZ5(I))/((E5(I)+PZ5(I))))
           else if(IITYP.gt.0) then
             NQP=NQP+1
c             print*,'Q px=',PX5(I),' py=',PY5(I),' pz=',PZ5(I),' rap=',
     &  -0.5*log((E5(I)-PZ5(I))/((E5(I)+PZ5(I))))
           else if(IITYP.lt.0) then
             NQPB=NQPB+1
c             print*,'QB px=',PX5(I),' py=',PY5(I),' pz=',PZ5(I),' rap=',
     &  -0.5*log((E5(I)-PZ5(I))/((E5(I)+PZ5(I))))
           endif

        if(ioscar.ge.1) 
     &     write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1       GX5(I), GY5(I), GZ5(I), FT5(I)
 1013   CONTINUE      
 211    format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
 210    format(3I8,f10.4,4I5)

c      print*,'total:',MUL,' Q:',NQP,' Qb:',NQPB,' g:',NGLU

      if(isoft.ne.4) return


      NATT=0
      EATT=0.
      call ptoh

      !no zpc particles gamma, e need to be filled, then call art
      do 1014 I=1,nnozpc
              NATT=NATT+1
              KATT(NATT,1)=ITYPN(I)
              PATT(NATT,1)=PXN(I)
              PATT(NATT,2)=PYN(I)
              PATT(NATT,3)=PZN(I)
              PATT(NATT,4)=EEN(I)
              EATT=EATT+EEN(I)
              GXAR(NATT)=GXN(I)
              GYAR(NATT)=GYN(I)
              GZAR(NATT)=GZN(I)
              FTAR(NATT)=FTN(I)
              ITYPAR(NATT)=ITYPN(I)
              PXAR(NATT)=PXN(I)
              PYAR(NATT)=PYN(I)
              PZAR(NATT)=PZN(I)
              PEAR(NATT)=EEN(I)
              XMAR(NATT)=XMN(I)
 1014 continue

      pxtmp=0
      pytmp=0
      pztmp=0
      petmp=0

      if(ioscar.ge.1) 
     & write(33,*) ievt,' ',1,' ',NATT,bimp,' ',0,' ',0,' ',0,' ',0,' ',
     & NELT,NINTHJ,0.0d0

      do i=1,NATT
        if(ioscar.ge.1) 
     &   write(33,*) ITYPAR(i),' ',PXAR(i),' ',PYAR(i),' ',PZAR(i),' ',
     &  PEAR(i),' ',GXAR(i),' ',GYAR(i),' ',GZAR(i),' ',
     &  FTAR(i)

        pxtmp=pxtmp+PXAR(i)
        pytmp=pytmp+PXAR(i)
        pztmp=pztmp+PXAR(i)
        petmp=petmp+PEAR(i)
      enddo

      !4-mom conservation test
c      print*,'momentum conservation test right after coal px:',pxtmp,
c     & ' py:',pytmp,' pz:',pztmp,' E:',petmp

      !art parameter need to read in before call art ini
      !number of tracks
      IAINT2(1) = NATT  

c.....ART initialization and run
      CALL ARINI
      CALL ARINI2(1)
      CALL ARTMN

      return
      end
