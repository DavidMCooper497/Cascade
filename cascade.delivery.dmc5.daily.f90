!*==AA0001.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
 module f77kinds
!
! This module was generated automatically by MKKIND.F90
!	  (supplied with the plusFORT toolkit from Polyhedron Software).
! It maps common non-standard Fortran 77 types to Fortran 90 kinds.
! Because different Fortran 90 compilers use different kind
! numbers, MKKIND.F90 should be compiled and run with the
! compiler you intend to use.
!
		integer,parameter :: &
			I1KIND  = 1,		& ! INTEGER*1
			I2KIND  = 2,		& ! INTEGER*2
			I4KIND  = 4,		& ! INTEGER*4
			L1KIND  = 1,		& ! LOGICAL*1
			L4KIND  = 4,		& ! LOGICAL*4
			R4KIND  = 4,		& ! REAL*4
			R8KIND  = 8,		& ! REAL*8
			DPKIND  = 8,		& ! DOUBLE PRECISION
			CX8KIND = 4,		& ! COMPLEX*8
			CX16KIND= 8 	  ! COMPLEX*16
		end module f77kinds
		
module i_common
      use f77kinds
      implicit none
!
! PARAMETER definitions
!
      integer,parameter :: mhru=1500,mstate=350,mlayer=3,mdet=30,    &
                         & msoil=30,mluse=30,mls=1000,mtime=6000,       &
                         & mchksite=20,mpsed=4,mpest=10,mlscape=30,     &
                         & ljopt=mtime*mdet,lvopt=15,liwopt=1,          &
                         & lwopt=6*lvopt+ljopt*lvopt+2*ljopt+           &
                         & lvopt*(lvopt-1)/2
!
! COMMON /FNAMES/
!
      character(1) :: aintopt,answer
		character(6) :: chireg
      character(4),dimension(3) :: anamechar
      character(20),dimension(30) :: anamedet
      character(20) :: maincat
      character(40) :: paramdir
      character(200) :: paramfile,fname1,findexname,findaccum,fout
      character(20),dimension(20) :: pestname
      character(180) :: rmain,rout
      character(60) :: root,rsub
!      character(60),dimension(20) :: subcat
!
! COMMON /IFLOW/
!
      integer :: i,evapskip,iagr,idailyevap,iday,iday1,igap,igeol,ihost,  &
               & ihostdom,ilsno,imbout,imonth,imonth1,ioptwt,ioutfbyhru,&
               & ihru,ihruno,isopt,istopping,isubstep,itime,ixcentr,    &
               & iycentr,iyear,iyear1,jsoil,julrainstart,julstart,mbchk,&
               & nag,nagluse,ncall,ndet,ndetmax,nfdet,nhru,nite,nlayer, &
               & nlm1,nlscape,nluse,nmorecs,noriv,nparam,npest,         &
               & npestarea,nrain,nredstate,nregion,nsdet,nsoil,nsubcat, &
               & ntime,rainskip,isims,idum,idd,idettot,isedtot,ifixtot,  &
					& k,ncdet,isubcat,jhru,id1,id2,id3,iparam,iflag, &
					& iflagopt,jdet,kdet,ksdet,kfdet,j,llay,idet,ideeper,ilscape, &
					& mm,ipcount,ipest,id,jlscape,itype,jco,mdum1,mdum2,jmon,iostat, &
					& jtype,infile,ico,imon,iyr,ihr,imin,isec,ksubcat,jjday, &
					& idh,kch,mch,it,ivy,ivm,ivd,ifl,ndays,isvdate,jeco,kl, &
					& iorgd,iorgm,iorgy,ndd,ifvmax,ipr,icmon,icday,icyr,igf
      integer,dimension(150) :: iaag,icag
      integer,dimension(mpest,150,12) :: iaccump
      integer,dimension(mtime) :: id2v,iyearv,imonthv,idayv,ihrv,iminv,isecv
      integer,dimension(50) :: idetindex
      integer,dimension(30) :: idoparam
      integer,dimension(mlscape) :: msedindex
      integer,dimension(12) :: ileap,jday=(/31,28,31,30,31,30,31,31,30,31,30,31/)
 !     integer,dimension(mchksite,mhru) :: indhru
      integer,dimension(mstate,mstate) :: inextstate
      integer,dimension(mhru) :: iregconstr,jhruarea,nflchk
      integer,dimension(mstate) :: istatereal,ntrans
      integer,dimension(mchksite) :: len
      integer,dimension(100) :: jea,js
      integer,dimension(200,mdet) :: mindex
      integer,dimension(200) :: nacc
      integer,dimension(3) :: ndacc
!
! COMMON /ISTSTEUPD/
!
      integer,dimension(mstate) :: jstarea
!
! COMMON /NMCH/
!
      character(1) :: aimbout,qopt
      character(60),dimension(mdet) :: anamedetl,bnamedet
!
! COMMON /RFLOW/
!
      real(R8KIND),dimension(mlscape) :: aksed,aovf,csedinit,phost,plusev,&
                & psoilv,rho
      real(R8KIND),dimension(mlscape,12,mdet) :: accum
      character(80) :: ach
      real(R8KIND) :: aend,aestart,aetot,aovfroot,arstart,astart,       &
                    & betavert,bovfroot,capinf,cbase,covf,cs,dg,dsmd,   &
                    & ovfin,psub,ps_correction,purb,sedmean,sm,soilarea,&
                    & soilareat,sumsed,thiscondsat,thiscsedinit, &
                    & thisdepth,ovlflow,    &
                    & thistemp,thiswcr,thiswcs,totae,validdata,vsed,wtpe, &
						  & thisrainr,thisper,urbmult,sumins,thisaer,respe, &
						  & amaxev,suminss,som,rhod,theta,sumpest,tprime,akh, &
						  & akp,sum1,sum2,gamma,vertm,available,potsnew, &
						  & sovf,beta,alatm,beta1,sssss,aroot,broot,val,arend,aeend, &
						  &  aday,sdconc,ttotsum,stepmx,eta,vertfl,potdrain
		real(R8KIND),dimension(60,40) :: agarea
      real(R8KIND),dimension(150) :: agarearev
      real(R8KIND),dimension(mdet) :: aintf,concrain,concrainred,       &
                                    & optsedwt,optfixwt,optsolwt,ovf,sconc,    &
                                    & sumbasef,sumbasefl0, toputsol, toputsed,       &
                                    & sumbasefl1,sumbasefl2,sumintf,    &
                                    & sumovf,totb,toti,toto,dopt
      real(R8KIND),dimension(mpest) :: akoc,pestparama,pestparamb,      &
                                     & pestparamc
      real(R8KIND),dimension(10) :: amorecswt,rainwt
      real(R8KIND),dimension(mstate,mdet) :: basef
      real(R8KIND),dimension(mlscape) :: bovf
      real(R8KIND),dimension(mlscape,mlayer,mdet) :: cinit
      real(R8KIND),dimension(mpest,150,mlayer) :: cinitp
!      real(R8KIND),dimension(mlscape,mlayer) :: condsat,depth,origdepth,&
!                & wcfc,wcr,wcs
      real(R8KIND),dimension(mlayer) :: fluxtobase
      real(R8KIND),dimension(mtime,mdet) :: dds,totsumf
      real(R8KIND),dimension(mtime) :: evap_array
      real(R8KIND),dimension(mlscape,mlayer+1,mdet) :: fixconc
      integer,dimension(100,100) :: imcoord,irdcoord
      real(R8KIND),dimension(mdet) :: ovfdel,ovfl1,ovfl2,fvecdum
      real(R8KIND),dimension(20) :: param
      real(R8KIND),dimension(150,mtime) :: pe
      real(R8KIND),dimension(30) :: pite
      real(R8KIND),dimension(30) :: totlusev
		real(R8KIND),dimension(mlscape,mlscape) :: pls
      real(R8KIND),dimension(mhru,mlscape) :: plscape
      real(R8KIND),dimension(mstate,mstate) :: pnextstate
      real(R8KIND),dimension(mhru) :: propeff,sflchk
      real(R8KIND),dimension(mpsed,mdet) :: psed
      real(R8KIND),dimension(mpsed,mlscape,mdet) :: psedbase
      real(R8KIND),dimension(mtime) :: rain_array,temp
      real(R8KIND),dimension(mdet,mlayer+1) :: toputfix
      real(R8KIND),dimension(mlscape,mdet) :: rsedbase
!      real(R8KIND),dimension(mlayer,mdet) :: s,thiscinit,snew
      real(R8KIND),dimension(mpest,mlayer) :: sadspest
!      real(R8KIND),dimension(5,30) :: sdum
      real(R8KIND),dimension(mstate) :: sedmax
      real(R8KIND),dimension(mls,mlayer,mdet) :: so
      real(R8KIND),dimension(mtime,mlscape,mdet) :: sumf
		real(R8KIND),dimension(mtime,mlscape,mlayer) :: smd
      real(R8KIND),dimension(mpest,mlscape,mlayer) :: thiscinitpads
!      real(R8KIND),dimension(mlscape,mdet,mlayer+1) :: thisfix
      real(R8KIND),dimension(mtime,mdet,mchksite) :: totsum,totsumbasef,&
                & totsumintf,totsumovf
      character(20),dimension(100) :: vch
!
! COMMON /RSTATEUPD/
!
      character(60),dimension(mdet) :: anameshort
      character(60) :: andum
      character(20),dimension(mchksite) :: asite
      character(200) :: fname
      real(R8KIND),dimension(mdet) :: rainsource
!
! COMMON /TPERIOD/
!
      integer :: eday,ehr,emin,emon,esec,eyea,sday,shr,smina,smon,ssec, &
               & syea
!
! COMMON /XUNITS/
!
      character(2),dimension(30) :: aunits,bunits
		type lstype
			real(R8KIND)								:: area,aovf,bovf,ambalrain,sumae,abase,bbase, &
																& aksed,ambalstart,ambalend,ambalroff,sumrain,sumpe, &
																& rho
			real(R8KIND),dimension(mlayer)		:: condsat,depth,wcr,wcfc,origdepth,wcs
			real(R8KIND),dimension(mtime,mdet)	:: sumf,smd
			integer,dimension(mdet)					:: mindex	
			integer										:: peflag	
		end type lstype

		type cattype
			character(60)								:: name
			integer										:: nflchk
			integer,dimension(mhru)					:: hruflag 
			real(R8KIND)								:: area,sflchk
			real(R8KIND),dimension(mtime,mdet)  :: totsum,totsumbasef,&
																& totsumintf,totsumovf,valid,validr
			end type cattype
		
		type layertype
			character(60)								:: name
			real(R8KIND)								:: sav,smin,smax,ddepth,latflow,scfc,sinit
			real(R8KIND),dimension(mdet)			:: s,snew,sdum,thiscinit 
		end type layertype
		
		type dettype
			character(2)					:: aunits 
			character(60)									:: name,shortname,namel
			real(R8KIND)									:: optwt 
			real(R8KIND),dimension(mlscape)			:: thisacc 
!			real(R8KIND),dimension(mlayer)			:: rfix
			real(R8KIND),dimension(mlscape,mlayer)	:: thisfix
			integer											:: solflag,sedflag,fixflag,addlayer,vartype,todo,&
																   & solindex,sedindex,fixindex,indexpest	
		end type dettype

		type detredtype
			character(2)				:: bunits
			character(60)				:: bname
			real(R8KIND)									:: ovf,intf,basef,sumovf,sumintf,sumbasef
			real(R8KIND),dimension(mlscape)			:: thisacc 
!			real(R8KIND),dimension(mlayer)			:: rfix
			real(R8KIND),dimension(mlscape,mlayer)	:: thisfix
			integer											:: solflag,sedflag,fixflag,addlayer,vartype,todo,&
																   & solindex,sedindex,fixindex,indexpest	
		end type detredtype
		
		type raintype
			real(R8KIND)									:: depthm,depthmm,urbmult,urban
			real(R8KIND),dimension(mdet)				:: conc,concred
		end type raintype	  	
		
		type evaptype
			real(R8KIND)									::	potm,potmm,actm
		end type evaptype

		type(lstype) lscape(mlscape),lsc
		type(dettype) det(mdet)
		type(cattype) subcat(30),subc
		type(layertype) layer(mlayer)
		type(raintype) rain
		type(evaptype) evap
		type(detredtype) detred(mdet)
		end module i_common
		
MODULE S_STATEUPD
   INTERFACE
      subroutine stateupd
      use f77kinds
      use i_common
      implicit none
      END SUBROUTINE STATEUPD
   END INTERFACE
END MODULE S_STATEUPD
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_TRANSFORM
   INTERFACE
      subroutine transform
      use f77kinds
      use i_common
      implicit none
      END SUBROUTINE TRANSFORM
   END INTERFACE
END MODULE S_TRANSFORM
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_TRANSFORM2
   INTERFACE
      subroutine transform2 !unadjusted for variable time step
      use f77kinds
      use i_common
      implicit none
      END SUBROUTINE TRANSFORM2
   END INTERFACE
END MODULE S_TRANSFORM2
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_VERTFLOW
   INTERFACE
      subroutine vertflow
      use f77kinds
      use i_common
      implicit none
      END SUBROUTINE VERTFLOW
   END INTERFACE
END MODULE S_VERTFLOW
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_LATFLOW
   INTERFACE
      subroutine latflow
      use f77kinds
      use i_common
      implicit none
      END SUBROUTINE LATFLOW
   END INTERFACE
END MODULE S_LATFLOW
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_STRUCTURE
   INTERFACE
      subroutine structure
      use i_common
      implicit none
      END SUBROUTINE STRUCTURE
   END INTERFACE
END MODULE S_STRUCTURE
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_GETRAINEVAP
   INTERFACE
      subroutine getrainevap
      use i_common
      implicit none
      END SUBROUTINE GETRAINEVAP
   END INTERFACE
END MODULE S_GETRAINEVAP
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_READRAINEVAP
   INTERFACE
      subroutine readrainevap
      use f77kinds
      use i_common
      implicit none
      END SUBROUTINE READRAINEVAP
   END INTERFACE
END MODULE S_READRAINEVAP
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_SEDGENNEW
   INTERFACE
      subroutine sedgennew
      use i_common
      implicit none
      END SUBROUTINE SEDGENNEW
   END INTERFACE
END MODULE S_SEDGENNEW
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_READTEMP
   INTERFACE
      subroutine readtemp
      use i_common
      implicit none
      END SUBROUTINE READTEMP
   END INTERFACE
END MODULE S_READTEMP
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_LSQFUN
   INTERFACE
      subroutine lsqfun(mopt,nopt,paramopt,fveccopt)
      use f77kinds
      use i_common
      implicit none
      integer :: mopt,nopt
      real(R8KIND),dimension(mopt) :: fveccopt
      real(R8KIND),dimension(nopt) :: paramopt
      intent (in) mopt,nopt,paramopt
      intent (inout) fveccopt
      END SUBROUTINE LSQFUN
   END INTERFACE
END MODULE S_LSQFUN
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_LSQMON
   INTERFACE
      subroutine lsqmon(mopt,nopt,paramopt,fveccopt,fjaccopt,ljcoptt,   &
                      & sopt,igrade,niteropt,nfopt,iwopt,liwoptt,wopt,  &
                      & lwoptt)
      use f77kinds
      use i_common
      implicit none
      integer :: igrade,liwoptt,ljcoptt,lwoptt,mopt,nfopt,niteropt,nopt
      real(R8KIND),dimension(ljcoptt,nopt) :: fjaccopt
      real(R8KIND),dimension(mopt) :: fveccopt
      integer,dimension(liwoptt) :: iwopt
      real(R8KIND),dimension(nopt) :: paramopt,sopt
      real(R8KIND),dimension(lwoptt) :: wopt
      intent (in) liwoptt,ljcoptt,lwoptt,mopt,nfopt,nopt
      END SUBROUTINE LSQMON
   END INTERFACE
END MODULE S_LSQMON
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_DATEY
   INTERFACE
      subroutine datey(nn,ids,ims,iys,idf,imf,iyf,ifl)
      implicit none
      integer :: idf,ids,ifl,imf,ims,iyf,iys,nn
      intent (in) ids,ifl,ims,iys
      intent (inout) idf,imf,iyf,nn
      END SUBROUTINE DATEY
   END INTERFACE
END MODULE S_DATEY
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_JULDAY
   INTERFACE
      function julday(mm,id,iyyy)
      implicit none
!
! PARAMETER definitions
!
      integer,parameter :: igrega=15+31*(10+12*1582)
      integer :: id,iyyy,mm
      integer :: julday
      intent (in) id,iyyy,mm
      END FUNCTION JULDAY
   END INTERFACE
END MODULE S_JULDAY
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_CALDAT
   INTERFACE
      subroutine caldat(julian,mm,id,iyyy)
      implicit none
!
! PARAMETER definitions
!
      integer,parameter :: igreg=2299161
      integer :: id,iyyy,julian,mm
      intent (in) julian
      intent (out) id
      intent (inout) iyyy,mm
      END SUBROUTINE CALDAT
   END INTERFACE
END MODULE S_CALDAT
!*==intfaces.f90  created by SPAG 6.55Dc at 15:04 on  6 Mar 2006
MODULE S_INITIALISE
   INTERFACE
      subroutine initialise(startt,endt)
      use i_common
      implicit none
      character(14) :: endt,startt
      intent (inout) endt,startt
      END SUBROUTINE INITIALISE
   END INTERFACE
END MODULE S_INITIALISE

	program main	
		
		
	INCLUDE 'link_fnl_shared.h'
	use f77kinds
      use i_common
      use s_initialise
      use s_julday
      use s_lsqfun
      use s_lsqmon
      use s_readrainevap
      use s_readtemp
      use s_sedgennew
      implicit none
!*--********************************************************************
! calls       E04FCF   INITIALISE        JULDAY   len_trim   LSQFUN
!             LSQMON   READRAINEVAP      READTEMP SEDGENNEW
! called by   ** NOTHING **
! modifies    /FNAMES/ /IFLOW/  /NMCH/   /RFLOW/  /RSTATEUPD/
!             name nameL         ANAMESHORT        FNAME
!             IADDLAYER         IDETFLAG IDETINDEX         IDOPARAM
!             IFIXFLAG IFIXINDEX         IFLAGPE  ILEAP    INDHRU
!             IOPTWT   IREGCONSTR        ISEDFLAG ISEDINDEX
!             ISTOPPING         ITOTAREA IVARTYPE JDAY     KDO
!             NCALL    NDET     NDETMAX  NFDET    NLSCAPE  NPARAM
!             NSDET    OPTSEDWT OPTWT    PARAM
! passes arg  RMAIN    SUBCAT
! uses value  /FNAMES/ /IFLOW/  /RFLOW/  /RSTATEUPD/       name
!             FNAME    IADDLAYER         IDETFLAG IDETINDEX
!             IDOPARAM IFIXFLAG IFIXINDEX         INDHRU   IOPTWT
!             IOUTFBYHRU        ISEDFLAG ISEDINDEX         ISOPT
!             IVARTYPE KDO      NDET     NDETMAX  NFDET    NLSCAPE
!             NPARAM   NSDET    NSUBCAT  NTIME    PARAM    RMAIN
!             SUBCAT
! local vars  ADAY     AEEND    AIM      AKH      AKP      ALATM
!             ALPHA    AMAXEV   AREAD    AREND    AROOT    AV
!             BETA     BETA1    BETA2    BROOT    CHIREG   D1
!             DOPT     ENDT     ETA      FILENAME FINDACCUM
!             FINDEXNAME        FJAC     FLCHK    FNAME1   FOUT
!             FSUMSQ   FVEC     FVECDUM  GAMMA    GRASS    I
!             IACCTOT  IARD     IARTHIS  ICDAY    ICMON    ICNTD
!             ICNTTHIS ICO      ICYR     ID       ID1      ID2
!             ID3      ID4      ID5      ID6      IDAYV    IDD
!             IDEEPER  IDET     IDETTOT  IDF      IDH      IDUM
!             IDUM1    IDUM2    IFAIL    IFIXTOT  IFL      IFLAG
!             IFLAGOPT IFVMAX   IGF      IHR      IHRUAREA IHRUBRIV
!             IHRULAND IHRURRIV IHRV     IJ       ILSCAPE  ILUSE    IM
!             IMIN     IMINV    IMON     IMONTHV  INFILE   IORGD
!             IORGM    IORGY    IOSTAT   IPARAM   IPCOUNT  IPEST
!             IPR      IPRINT   IR       IRCAT    IRD      IRDAY
!             IRX      IRY      IRYEAR   ISEC     ISECV    ISEDTOT
!             ISIMS    ISITE    ISOIL    ico  ISVDATE  IT
!             ITE      ITOTSQNUM         ITYPE    IVD      IVM
!             IVY      IWOPT    IXM      IYM      IYR      J
!             JCO      JDET     JECO     JH       JHRU     JJ
!             JJDAY    JLSCAPE  JMON     JRYEAR   JTIME    JTYPE
!             JULRAINSTARTJSOIL K        KCH      KD       KDET
!             KDUM     KFDET    KH       KL       KSD      KSDET
!             KSUBCAT  L        LLAY     LUSE     MAXCAL   MCH
!             MDUM1    MDUM2    MM       MODFILE  MONTH    MOPT     MT
!             NCDET    NCHKSITE NDAYS    NDD      NEWMONTH NEWYEAR
!             NFLCHK   NFOPT    NITEROPT NOPT     PHOSTDOM PITEDUM
!             PLUSEDUM PLUSEVD  RESPE    RHOD     RRAIN    SAV
!             SAV1     SAV2     SDCONC   SDUMV    SFLCHK   SITE
!             SLAND    SMD      SNEW     SNEWPOT  SOM      SOPT
!             SOUT     SOVF     SSSSS    STARTT   STEPMX   SUM
!             SUM1     SUM2     SUMINS   SUMINSS  SUMPEST  SUMTILL
!             SYEAR    THETA    THISAER  THISPER  THISRAINR
!             TITLE1   TITLE2   TOTLUSEV TOTPN    TOTSUMF  TPRIME
!             TTOTSUM  URBAN    URBMULT  VAL      VERTM    VOPT
!             WOOD     WOPT     XOPT     XTOL
! uses PARAMs LIWOPT   LJOPT    LVOPT    LWOPT    MCHKSITE MDET
!             MLAYER   MLS      MLSCAPE  MPEST    MPSED    MREGION
!             MSTATE   MTIME
!*++********************************************************************
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
      real(R8KIND),intrinsic					:: dexp,dlog
      real(R8KIND),dimension(ljopt,lvopt) :: fjac
      real(R8KIND)								:: fsumsq,xtol
      real(R8KIND),dimension(ljopt)			:: fvec,fscale
      real(R8KIND),dimension(lvopt)			:: sopt,xopt,xscale,xout
      real(R8KIND),dimension(lvopt,lvopt) :: vopt
      real(R8KIND),dimension(lwopt)			:: wopt
		real(R8KIND),dimension(7)				:: rparam
      integer										:: ifail,iprint,maxcal,mopt,nfopt,niteropt,nopt
      integer,dimension(liwopt)				:: iwopt
		integer,dimension(6)						:: iparam2
 		character(14)								:: endt,startt

!*** End of declarations rewritten by SPAG
!
!*************************************************************************
!
! ----- CASCADE.DELIVERY.GENERAL
!
! ----- Delivery component computing HRU outputs
!       This version works with landscape units
!
! D.M.Cooper, B.Gannon and C-G.Wang
! CEH Wallingford
!
!      analytical solution of non-linear soil equations=
!
! f77 -w cascade.delivery.dmc2.f -o cascade.delivery.dmc2 -L$NAGDIR -lnag
!
!      This version uses only the dominant soil class in each HRU to compute
!      hydrology. It uses only 2 conceptual subsurface layers.
!
!      This combined version may be used for optimisation if requested,
!      and has options to halt optimisation and print the most recent simulations
!      Choice of weights by determinand needs to be coded for each
!      application
!
		do ico=1,30
			totlusev(ico)=0
		enddo
      chireg='      '
!
! Control file
!   filesource
!
! Input files dependent on catchment geometry:
!   <maincat>/hrubylscape
!   <maincat>/hruarea.data
!
! Input files dependent on landscape classification:
!   <maincat>/peweights.lscape
!   <maincat>/hrubylscape
!   <maincat>/lscape.soil.params
!   <maincat>/accum.lscape.<determ>
!   <maincat>/index.lscape.<determ>
!   <maincat>/fixconc.lscape.<determ>
!
! Input files dependent on determinands selected:
!   <maincat>/model.sim.dets
!   <maincat>/rain.conc
!
! Input files dependent on process description:
!   <maincat>/modelparam.data.general
!   <maincat>/accum.lscape.<determ>
!   <maincat>/fixconc.lscape.<determ>
!   <maincat>/rain.conc
!
! Input files for met. data
!   <maincat>/evap
!   <maincat>/modelrain.lscape.data
!   <maincat>/daily.temp
!
!
      call initialise(startt,endt)
!
! ----- open main output file
!
      open (40,file=rmain(1:len_trim(rmain))//'/modelout.chk.opt')
!
! ----- read determinand names
!
      open (70,file=rmain(1:len_trim(rmain))//'/fullnames')
      ico=1
      do
        read (70,'(a40,a10)',end=100) det(ico)%namel,                &
                                     & det(ico)%shortname
		  ico=ico+1
      enddo
 100  close (70)
!finished  reading names
!
! ----- Identify selected determinands with fixed internal ordering
!
      open (70,file=rmain(1:len_trim(rmain))//'/delivdets')
      do ico=1,12
         read (70,*)
      enddo
!
      ico=1
      do
         read (70,'(a10,8i5)',end=200) det(ico)%name,idum,idd,det(ico)%solflag,&
                                     & det(ico)%sedflag,det(ico)%fixflag,         &
                                     & det(ico)%addlayer,det(ico)%vartype,det(ico)%todo
         write (6,'(a10,8i5)') det(ico)%name,idum,idd,det(ico)%solflag,&
                                     & det(ico)%sedflag,det(ico)%fixflag,         &
                                     & det(ico)%addlayer,det(ico)%vartype,det(ico)%todo
         ico=ico+1
      enddo
!
 200  close (70)
      ndetmax=ico-1
!
! ----- read optimisation options
!
      open (70,file=rmain(1:len_trim(rmain))//'/optimisation')
      read (70,*) ioptwt
      write (6,*) ioptwt
      do ico=2,ndetmax
         read (70,*) idum,det(ico)%optwt
      enddo
      close (70)
!
! ----- initialise indexes
!
      do ico=2,ndetmax
         det(ico)%solindex=0
         det(ico)%sedindex=0
         det(ico)%fixindex=0
      enddo
!
      idettot=1
      isedtot=1
      ifixtot=0
      do ico=2,ndetmax
         if (det(ico)%solflag==1.and.det(ico)%todo==1) then
            idettot=idettot+1
            det(ico)%solindex=idettot
            write (6,'("ico,det(ico)%solindex",2i5)') ico,det(ico)%solindex
         endif
         if (det(ico)%sedflag==1.and.det(ico)%todo==1) then
            if (ico==11) then
               det(ico)%sedindex=1
            else
               isedtot=isedtot+1
               det(ico)%sedindex=isedtot
            endif
         endif
         if (det(ico)%fixflag==1.and.det(ico)%todo==1) then
            ifixtot=ifixtot+1
            det(ico)%fixindex=ifixtot
         endif
      enddo
      if (isedtot>1.and.det(11)%sedindex/=1) then
         write (*,'("Must have sediment included as determinand")')
         stop
      endif
      ncdet=0
      nsdet=0
      nfdet=0
      do ico=2,ndetmax
         if (det(ico)%solindex>=1) then
            ncdet=ncdet+1
				idum=det(ico)%solindex
            optsolwt(idum)=det(ico)%optwt
         endif
         if (det(ico)%sedindex>=1) then
            nsdet=nsdet+1
            idum=det(ico)%sedindex
				optsedwt(idum)=det(ico)%optwt
         endif
         if (det(ico)%fixindex>=1) then
            nfdet=nfdet+1
            idum=det(ico)%fixindex
				optfixwt(idum)=det(ico)%optwt
         endif
      enddo
      ndet=ncdet+1
! close(70) 
!finish reading <maincat>/order.params.lscape
!
! ----- identify HRUs to be modelled, using subcatchment information
!
      do ico=1,nsubcat
         fname=rmain(1:len_trim(rmain))//'/'//subcat(ico)%name             &
             & (1:len_trim(subcat(ico)%name))//'/subchk.data'
         open (70,file=fname(1:len_trim(fname)))
         jhru=0
         do
            jhru=jhru+1
            read (70,*,end=250) id1,id2,id3,subcat(ico)%hruflag(jhru)
			enddo
 250     close (70)
                   !finish reading <subcat>/subchk.data
      enddo
!
! ----- read main unfixed parameters
!
      fname=rmain(1:len_trim(rmain))//'/modelparam.data.general'
      fout=rout(1:len_trim(rout))//'/modelparam.data.general'
		open (60,file=fname(1:len_trim(fname)))
      open (61,file=fout(1:len_trim(fout))//'.forsplus')
      read (60,*) nparam
      write (61,*) nparam
      write (*,'("Parameter values and status this run")')
      do ico=1,nparam
         read (60,*) param(ico),idoparam(ico)
         write (61,*) param(ico),idoparam(ico)
         if (isopt==1.and.idoparam(ico)==0)                          &
            & write (*,'(i2,". ",f15.7,"      Fixed")') ico,         &
           & param(ico)
         if (isopt==0) write (*,'(i2,". ",f15.7," Fixed")') ico,     &
                            & param(ico)
         if (isopt==1.and.idoparam(ico)/=0)                          &
            & write (*,'(i2,". ",f15.7,"  Estimated")') ico,         &
           & param(ico)
         param(ico)=dlog(param(ico))
      enddo
      if (isopt==0.and.ioutfbyhru==1) then
         do ico=1,nsubcat
            fname1=rmain(1:len_trim(rmain))//'/'//subcat(ico)%name         &
                 & (1:len_trim(subcat(ico)%name))//'/output/Alatest.readme'
! write(6,*)fname1
            open (62,file=fname1)
            write (62,'("Parameter values for latest HRU output")')
            do iparam=1,nparam
               write (62,'(i2,". ",f15.7," Fixed")') iparam,            &
                    & dexp(param(iparam))
            enddo
            close (62)
         enddo
      endif
      close (60)
      close (61)
!
! ----- read series of driving variables
!
      write (6,'("Reading met data")')
      call readrainevap
      open (60,file=rmain(1:len_trim(rmain))//"/peweights.lscape")
      read (60,*) nlscape
      do ico=1,nlscape
         read (60,*) lscape(ico)%peflag
      enddo
      close (60)
      call readtemp
!
      write (6,'("Finished reading met data")')
!
! ----- define sediment concentration parameters by landscape
!
      call sedgennew
!
! ----- initialise date information
!
      do ico=1,12
         ileap(ico)=0
      enddo
!
! ----- setup soil and landuse parameters
!
!
! ----- Initialise
!
!
      do ico=1,nsubcat
         subcat(ico)%area=0
      enddo
 
! ----- prepare for optimisation
!
      iflag=0
      iprint=1
      maxcal=300
! eta=.5
      xtol=0.0D0
! stepmx=1.d2
      ifail=1
      mopt=ntime*ndet
      nopt=0
      ncall=0
      do ico=1,nparam
         if (idoparam(ico)>0) then
            xopt(idoparam(ico))=param(ico)
            nopt=nopt+1
         endif
      enddo
!
! ----- Either optimise or compute
!
!      write(6,'("Check")')
!      stop
      if (isopt==1) then
         istopping=0
         write (6,'("calling optimisation")')
!         call e04fcf(mopt,nopt,lsqfun,lsqmon,iprint,maxcal,eta,xtol,    &
!                   & stepmx,xopt,fsumsq,fvec,fjac,ljopt,sopt,vopt,lvopt,&
!                   & niteropt,nfopt,iwopt,liwopt,wopt,lwopt,ifail)
         call dunlsf(lsqfun,mopt,nopt,xopt,xscale,fscale,iparam2,rparam,xout,    &
			          &fvec,fjac,ljopt)
			
			write (*,'("Final value of ifail and iprint",2i8)') ifail,     &
              & iprint
         write (*,'("Final Parameter estimates")')
         istopping=1
         iflagopt=-1
         call lsqfun(mopt,nopt,xopt,fvec)
         do ico=1,nparam
            write (*,'(i2,". ",f20.12)') iparam,dexp(param(ico))
         enddo
         if (ioutfbyhru==1) then
            do ico=1,nsubcat
               fname1=rmain(1:len_trim(rmain))//'/'//subcat(ico)%name      &
                    & (1:len_trim(subcat(ico)%name))                       &
                     &//'/output/Alatest.readme'
               open (62,file=fname1)
               write (62,'("Parameter values for latest HRU output")')
               do jco=1,nparam
                  write (62,'(i2,". ",f20.12)') jco,                 &
                       & dexp(param(jco))
               enddo
               close (62)
            enddo
         endif
         stop
      elseif (isopt==0) then
         istopping=1
         iflagopt=-1
         write (6,'("calling lsqfun")')
         call lsqfun(mopt,nopt,xopt,fvec)
         write (6,'("out of lsqfun")')
      endif
!
      end
!*==STATEUPD.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!********************************************************************************
      subroutine stateupd
      use f77kinds
      use i_common
      use s_latflow
      use s_transform
      use s_transform2
      use s_vertflow
      implicit none
!*--STATEUPD374
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
      real(R8KIND),intrinsic :: dexp,dmax1,dmin1
!
!*** End of declarations rewritten by SPAG
!
!
! ----- For a given determinand, update states and compute fluxes this time step
!       Indices are:
!
!           j - lateral positional state within region
!           k - vertical positional state within region
!           l - interpenetrating state at given location
!
!     Assumes a unique lateral and vertical division, and time stepping
!
!
      if (ndet>23) then
         write (*,'("check array bounds in stateupd")')
         stop
      endif
!
! ----- modify accumulation rates according to prevailing conditions
!		  internal changes in concentration computed at start of time step 
!
      call transform
!
! ----- add proportional components of transformed variables to soil region
!       picking off only those of interest
!
      do jdet=2,ndetmax
! ----- determinands not associated with sediment
         kdet=det(jdet)%solindex
         if (kdet>=1) toputsol(kdet)=det(jdet)%thisacc(ilsno)
! ----- sediment determinands
         ksdet=det(jdet)%sedindex
         if (ksdet>1) toputsed(ksdet)=det(jdet)%thisacc(ilsno)
! ----- fixed concentration determinands
         kfdet=det(jdet)%fixindex
         if (kfdet>=1) then
            do llay=1,nlayer
               toputfix(kfdet,llay)=det(jdet)%thisfix(ilsno,llay)
            enddo
         endif
      enddo
!
! ----- add non-sediment accumulations to soil
!       toputsol is in units of kg/ha/time step
!       s is in units of kg per HRU
!
      do jdet=2,ndetmax
         kdet=det(jdet)%solindex
         if (kdet>=1) then
            llay=det(jdet)%addlayer
            if (toputsol(kdet)>0) then
               layer(llay)%s(kdet)=dmax1(0.D0,layer(llay)%s(kdet)+toputsol(kdet)/(100*100)&
                          & )                   !convert from kg/ha to kg/m**2
            else
               layer(llay)%s(kdet)=dmax1(0.D0,layer(llay)%s(kdet))
            endif
         endif
      enddo
!
!
! ----- do pesticide internal transformations 
!
      call transform2
		
      if (imbout==1) write (40,'(//"Date",3i5/"Landscape number",i6)')  &
                          & iday,imonth,iyear,ilsno
      do idet=1,ndet
         detred(idet)%sumovf=0.
      enddo
      totae=0.

	   do llay=1,nlayer
         do kdet=1,ndet
            layer(llay)%snew(kdet)=layer(llay)%s(kdet) !start to update soil contents
         enddo
      enddo

      do kdet=1,ndet
         detred(kdet)%ovf=0.
      enddo

!$$$$$$$$$$$$$$$$$$$$$$$$$ PARAM(1)  $$$$$$$$$$$$$$$$$$$$$$$$$
 
      rain%depthm=rain%depthmm/1000   !rescale to m**3/m**2 --- all calculations done per square metre
      evap%potm=evap%potmm/1000 !rescale to m**3/m**2
!
! ----- compute proportion of "urban" rainfall. This is in effect rapid response runoff
!
      rain%urbmult=10*dexp(param(1))/(1+dexp(param(1)))
      rain%urban=rain%depthm*rain%urbmult !compute "urban"rainfall
      rain%depthm=rain%depthm-rain%urban !compute remaining rainfall
!
! ----- compute maximum infiltration capacity
!
      capinf=lsc%condsat(1)*isubstep/(24*60.*60.) !m per time step from m/day
      ovlflow=rain%depthm-dmin1(rain%depthm,capinf)
      rain%depthm=rain%depthm-ovlflow
      layer(1)%snew(1)=layer(1)%s(1)+rain%depthm ! add sources which vary geographically
!
! ----- compute actual evaporation and remove from soil
!      
		call compute_actevap
!
! ----- generate overland flow as sum of "urban" and infiltration excess
!
      detred(1)%ovf=rain%urban+ovlflow
      do kdet=2,ndet
         detred(kdet)%ovf=(rain%urban+ovlflow)*rain%concred(kdet)/1000. !kg/m^2
      enddo
      if (layer(1)%snew(1)>layer(1)%smax) then
         do kdet=2,ndet
            detred(kdet)%ovf=detred(kdet)%ovf+(layer(1)%snew(1)-layer(1)%smax)*rain%concred(kdet)   &
                    & /1000.
            layer(1)%snew(kdet)=layer(1)%s(kdet)+(layer(1)%smax-layer(1)%s(1))*rain%concred(kdet)   &
                       & /1000.
         enddo
         detred(1)%ovf=detred(1)%ovf+layer(1)%snew(1)-layer(1)%smax
         layer(1)%snew(1)=layer(1)%smax
      else
         do kdet=2,ndet
            layer(1)%snew(kdet)=layer(1)%s(kdet)+rainsource(kdet)
         enddo
      endif
!
! ----- lag overland flow
!
      covf=1-lsc%aovf-lsc%bovf
      do kdet=1,ndet
         ovfin=amin1(1.,(itime-1.)/10.)*detred(kdet)%ovf 
         if (itime>2) ovfl2(kdet)=ovfl1(kdet)
         if (itime>1) ovfl1(kdet)=ovfdel(kdet)
         ovfdel(kdet)=lsc%aovf*ovfl1(kdet)+lsc%bovf*ovfl2(kdet)   &
                    & +covf*ovfin
      enddo
!
      do kdet=1,ndet
         detred(kdet)%sumovf=detred(idet)%sumovf+ovfdel(kdet)
      enddo
 
      do kdet=1,nfdet
         detred(ndet+nsdet+kdet)%sumovf=detred(1)%sumovf*toputfix(kdet,1)/1000 !!conc in kg/m**3
      enddo
!
      do llay=1,nlayer
         do kdet=1,ndet
            layer(llay)%s(kdet)=layer(llay)%snew(kdet)
			enddo
      enddo
	  	
      lsc%ambalrain=lsc%ambalrain+(rain%depthm+rain%urban-evap%actm)
                                                                !m^3 per m^2
      totae=totae+evap%actm
      lsc%sumae=lsc%sumae+totae*1000. !mm 
!
! ----- subsurface water transfers
!
      call vertflow
      call latflow
!
      if (ilsno==-5) write (6,                                          &
     &'("soil water content layers 1 and 2"                             &
     &                                              ,4f10.4)') layer(1)%smax, &
    & layer(1)%s(1),layer(2)%smax,layer(2)%s(2)
      end
!*==TRANSFORM.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!***************************************************************
      subroutine transform
      use f77kinds
      use i_common
      implicit none
!*--TRANSFORM616
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
      real(R8KIND),intrinsic :: dexp
!
!*** End of declarations rewritten by SPAG
!
!
! ----- generate accumulation rates this time step
!
 
!
! ----- indices are: landuse, time of year, determinand
!
! ----- pH
!       Not modelled
!
!		do ilsno=1,nlscape 
!			lsc=lscape(ilsno)
		   if (det(2)%solindex>=1) then
            mm=lsc%mindex(2)
            det(2)%thisacc(ilsno)=accum(mm,imonth,2)                       &
                             & *isubstep/(30*24*60*60.)
         endif
!
! ----- Cl
!       Accumulation rate zero; rainfall concentration 5mg/l
!
			if (det(3)%solindex>=1) then
	         mm=lsc%mindex(3)
            det(3)%thisacc(ilsno)=accum(mm,imonth,3)                       &
                             & *isubstep/(30*24*60*60.)
	      endif
!
! ----- Temperature
!       Not modelled
!
			if (det(4)%solindex>=1) then
            mm=lsc%mindex(4)
            det(4)%thisacc(ilsno)=accum(mm,imonth,4)                       &
                             & *isubstep/(30*24*60*60.)
         endif
!
! ----- Dissolved Oxygen
!       Not modelled
!
	      if (det(5)%solindex>=1) then
            mm=lsc%mindex(5)
            det(5)%thisacc(ilsno)=accum(mm,imonth,5)                       &
                             & *isubstep/(30*24*60*60.)
		   endif
!
! ----- Nitrate
!       Crop-dependent accumulation rate. Nominal zero concentration
!       in rainfall
!
	      if (det(6)%solindex>=1) then
            mm=lsc%mindex(6)
            det(6)%thisacc(ilsno)=accum(mm,imonth,6)                       &
                             & *isubstep/(30*24*60*60.)
! write(6,'("mm,ilsno,thisacc",2i6,f10.4)')mm,ilsno,
! &           thisacc(ilsno,6)
         endif
!
! ----- Ammonia
!       Not modelled
!
			if (det(7)%solindex>=1) then
            mm=lsc%mindex(7)
            det(7)%thisacc(ilsno)=accum(mm,imonth,7)                       &
                             & *isubstep/(30*24*60*60.)
         endif
!
! ----- Particulate organic N
!       .01 * sediment concentration except October
!       Follows particulate organic carbon
!
			if (det(8)%sedindex>=1) then
            mm=lsc%mindex(8)
            if (imonth/=10) then
               det(8)%thisacc(ilsno)=accum(mm,imonth,8)                    &
                                & *isubstep/(24*60*60.)
                         !base value is daily
            else
               det(8)%thisacc(ilsno)=.1*(30-.089*((iday+273)-289)**2)/100
               det(8)%thisacc(ilsno)=det(8)%thisacc(ilsno)                    &
                                & *isubstep/(24*60*60.)
            endif
         endif
!
! ----- Orthophosphate
!
!
		   if (det(9)%fixindex>=1) then
! write(6,'("Fixed concs for phosphate")')
            mm=lsc%mindex(9)
            do llay=1,nlayer+1
               det(9)%thisfix(ilsno,llay)=fixconc(mm,llay,9)
            enddo
! write(6,'(i4,3f10.4)')ilsno,(thisfix(ilsno,9,llay),llay=
! &           1,nlayer+1)
	      endif
!
! ----- Particulate P
!       .0025 * sediment concentration
!
			if (det(10)%sedindex>=1) then
            mm=lsc%mindex(10)
            det(10)%thisacc(ilsno)=accum(mm,imonth,10)
         endif
!
! ----- Suspended sediment
!       multiple of flow, depending on land use (see sedgennew)
!
			if (det(11)%sedindex>=1) then
            mm=lsc%mindex(11)
            det(11)%thisacc(ilsno)=accum(mm,imonth,11)
         endif
!
! ----- Nickel
!       .5 * sediment concentration
!
			if (det(12)%sedindex>=1) then
            mm=lsc%mindex(12)
            det(12)%thisacc(ilsno)=accum(mm,imonth,12)
         endif
!
! ----- Lead
!       .2 * sediment concentration
!
			if (det(13)%sedindex>=1) then
            mm=lsc%mindex(13)
            det(13)%thisacc(ilsno)=accum(mm,imonth,13)
         endif
!
! ----- SUP
!
			if (det(14)%fixindex>=1) then
            mm=lsc%mindex(14)
            do llay=1,nlayer+1
               det(14)%thisfix(ilsno,llay)=fixconc(mm,llay,14)
            enddo
         endif
!
! ----- Dissolved organic carbon
!       Accumulation dependent on soils as below
!
!     write(6,'("doing DOC")')
			if (det(15)%solindex>=1) then
            mm=lsc%mindex(15)
            sm=layer(1)%s(1)*100/layer(1)%smax
            if (sm<70.03) then
               cs=dexp(.036*sm)-1
            else
               cs=1747.8*dexp(-0.0718*sm)
            endif
! ----- cs*.298 is DOC production in g/day/m**3. Convert to kg/ha/day
            det(15)%thisacc(ilsno)=cs*.298*accum(mm,imonth,15)             &
                              & *(lscape(ilsno)%depth(1)+0.*lscape(ilsno)%depth(2)) &
                              & *.1
            det(15)%thisacc(ilsno)=det(15)%thisacc(ilsno)*isubstep/(24*60*60.)
                      !adjust for variable time stepn
! write(6,'("thisacc(mm,15)",f20.10)')thisacc(mm,15)
         endif
!
! ----- Particulate organic carbon
!       .1*sediment concentration, except October (below)
!
			if (det(16)%sedindex>=1) then
            mm=lsc%mindex(16)
            if (imonth/=10) then
               det(16)%thisacc(ilsno)=accum(mm,imonth,16)
            else
               det(16)%thisacc(ilsno)=(30-.089*((iday+273)-289)**2)/100
               det(16)%thisacc(ilsno)=det(16)%thisacc(ilsno)                  &
                                 & *isubstep/(24*60*60.)
            endif
			endif
!
! ----- Dissolved inorganic carbon
!       Accumulation depends on geology
!
			if (det(17)%solindex>=1) then
            mm=lsc%mindex(17)
            det(17)%thisacc(ilsno)=accum(mm,imonth,17)                     &
                              & *isubstep/(30*24*60*60.)
         endif
!
! ----- Phytoplankton
!       Not modelled
!
			if (det(18)%solindex>=1) then
            mm=lscape(ilsno)%mindex(18)
            det(18)%thisacc(ilsno)=accum(mm,imonth,18)                     &
                              & *isubstep/(30*24*60*60.)
         endif
!
!
! ----- Pesticide
!       Regional with processes
!
			do ipcount=19,ndetmax
				if (det(ipcount)%solindex>=1) then
               mm=lsc%mindex(ipcount)
               ipest=det(ipcount)%indexpest
               det(ipcount)%thisacc(ilsno)=iaccump(ipest,mm,imonth)        &
                                      & *isubstep/(30*24*60*60.)
            endif
			enddo
!		lscape(ilsno)=lsc
!		enddo
!
! write(6,'("leaving transform")')
      end
!*==TRANSFORM2.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!***************************************************************
      subroutine transform2 !unadjusted for variable time step
      use f77kinds
      use i_common
      implicit none
!*--TRANSFORM2877
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
      real(R8KIND),intrinsic :: dexp
!
!*** End of declarations rewritten by SPAG
!
!
! ----- carry out within-soil transformations
!
 
!
! ----- Pesticides
!
!       Pesticide is initially added to solution, and must be redistributed
      do ipcount=19,ndetmax
         kdet=det(ipcount)%solindex
         ipest=det(ipcount)%indexpest
         if (ipest>=1) then
            do llay=1,nlayer
               som=0
               if (llay==1) som=0.112*accum(ilsno,1,15)
               rhod=lsc%rho
               theta=layer(llay)%s(1)/lsc%depth(llay)
!
! ----- compute new value of total pesticide, to be partitioned using Kp=Cs/Cw
!       Note sumpest=Cw*theta*totvol+ Cs*rho*totvol
 
               sumpest=sadspest(ipest,llay)+layer(llay)%s(kdet)
                             ! units of sumpest are total mass of pesticide in the HRU
               tprime=.1204*pestparamc(ipest)*(20-thistemp)             &
                    & /(293*(thistemp+273))
               akh=pestparama(ipest)*(theta*100)**pestparamb(ipest)     &
                 & *dexp(tprime)
!
! ----- let it decay (nb -.693 = log .5)
!
!              if(llay.eq.1)write(6,'("sumpest pre",5f10.3)')sumpest,
!    &              theta,tprime,akh,thistemp
               sumpest=sumpest*dexp(-.693/akh)
! if(llay.eq.1)write(6,'("sumpest post",f20.10)')sumpest
               akp=akoc(ipest)*som
!
! ----- distribute between phases
!
               sadspest(ipest,llay)=sumpest*rhod*akp/(rhod*akp+theta)
               layer(llay)%s(kdet)=sumpest*theta/(rhod*akp+theta)
            enddo
! enddo
         endif
      enddo
!
      end
!*==VERTFLOW.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!***************************************************************
      subroutine vertflow
      use f77kinds
      use i_common
      implicit none
!*--VERTFLOW945
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
      real(R8KIND),intrinsic :: dexp,dmin1
!
!*** End of declarations rewritten by SPAG
!
!
! ----- move water and solute vertically. Linear in area, non-linear
!       with depth
!
!
      if (ndet>26) then
         write (*,'("check array bounds in vertflow")')
         stop
      endif
!
! ----- initialise
!
      sum1=0.
      sum2=0.
      sum1=sum1+detred(mbchk)%ovf
      sum2=sum2+detred(mbchk)%ovf
      do llay=1,nlayer
         do kdet=1,ndet
            layer(llay)%snew(kdet)=layer(llay)%s(kdet) !set layer content at beginning of timestep
         enddo
         sum1=sum1+layer(llay)%s(mbchk)
      enddo
! enddo
!
! ----- power law parameters
!
!$$$$$$$$$$$$$$$$$$$$$ PARAM 8 $$$$$$$$$$$$$$$$$$$$$$$$
      gamma=dexp(param(8))  
!$$$$$$$$$$$$$$$$$$$$$ PARAM 2 $$$$$$$$$$$$$$$$$$$$$$$$
      beta=dexp(param(2))
      sconc(1)=1.
      nlm1=nlayer-1
      do llay=nlm1,1,-1
         vertm=lsc%condsat(llay) !set saturated hydraulic conductivity
         if (llay==nlm1) then
            beta=beta
         else
            beta=gamma
         endif
!
         betavert=beta
         do kdet=1,ndet
            sconc(kdet)=layer(llay)%s(kdet)/layer(llay)%s(1) !set layer concentration at beginning of time step  
         enddo
!
! ----- If there is drainable water, estimate the volume drainable in a
!       time step. Move this down a layer assuming (i) There is space for it
!       (ii) The amount actually drained does not exceed the drainable
!       volume initially present
!
!
         if (layer(llay)%snew(1)-layer(llay)%smin>=0.) then ! there is some available water in the layer
!
! ----- 010799 Analytical solution, scaled by area
!
            available=layer(llay)%snew(1)-layer(llay)%smin ! available water

            if (beta/=0) potsnew=layer(llay)%smin                             & !snewpot is the potential new water content
                               & +(available**(-beta)+vertm*gamma*beta*isubstep      &
                               & /(24*60*60.))**(1.D0/(-beta))
            if (beta==0) potsnew=layer(llay)%smin                             &
                               & +available*dexp(-vertm*gamma*isubstep/              &
                               & (24*60*60.))
            potdrain=layer(llay)%snew(1)-potsnew ! potential drainage

            vertfl=dmin1(potdrain,layer(llay)%snew(1)-layer(llay)%smin, &
					&layer(llay+1)%smax-layer(llay+1)%snew(1)) !cannot drain more than available water, or more than can be accommodated in next layer

            do kdet=1,ndet
               layer(llay+1)%snew(kdet)=layer(llay+1)%snew(kdet)+vertfl*sconc(kdet) !move the water to next layer down
               if (llay+1==nlayer) fluxtobase(kdet)=vertfl*sconc(kdet)
               layer(llay)%snew(kdet)=layer(llay)%snew(kdet)-vertfl*sconc(kdet) !take the water from the present layer
            enddo
         endif
      enddo
!
      do llay=1,nlayer
         do kdet=1,ndet
            layer(llay)%s(kdet)=layer(llay)%snew(kdet)
         enddo
         sum2=sum2+layer(llay)%s(mbchk)
      enddo
! enddo
!
      if (imbout==1) write (40,                                         &
     &'("mass balance in vertflow"/                                     &
     &                                         "before/after"/,2f15.6)')&
    & sum1,sum2
      end
!*==LATFLOW.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!**************************************************************************
      subroutine latflow
      use f77kinds
      use i_common
      implicit none
!*--LATFLOW1067
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
      real,intrinsic :: amin1
      real(R8KIND),intrinsic :: dexp
!
!*** End of declarations rewritten by SPAG
!
!
! ----- move water and solute laterally. Note this is linear in area,
!       non-linear in depth
!
! ----- initialise updated state
!
      if (ndet>23) then
         write (*,'("check array bounds in latflow")')
         stop
      endif
      sum1=0.
      sum2=0.
      sovf=0.
      do llay=1,nlayer
         do kdet=1,ndet
            layer(llay)%sdum(kdet)=0.D0
         enddo
      enddo
      sum1=sum1+detred(mbchk)%ovf
      do llay=1,nlayer
         do kdet=1,ndet
            layer(llay)%snew(kdet)=layer(llay)%s(kdet)
         enddo
         sum1=sum1+layer(llay)%s(mbchk)
      enddo
! enddo
!
! ----- compute material flux from state j and to state jd, with
!       selected power law and saturation conductivity parameters
!
      sovf=sovf+detred(mbchk)%ovf
      do kdet=1,ndet
         detred(kdet)%intf=0.D0
			detred(kdet)%basef=0.d0
      enddo
      do llay=1,nlayer
!$$$$$$$$$$$$$$$$$$$$$$$$$ param 3 $$$$$$$$$$$$$$$$$$$$$$$$$$
         beta=dexp(param(3))
!$$$$$$$$$$$$$$$$$$$$$$$$$ param 4 $$$$$$$$$$$$$$$$$$$$$$$$$$
         if (llay==nlayer) beta=dexp(param(4))
!$$$$$$$$$$$$$$$$$$$$$$$$$ param 5 $$$$$$$$$$$$$$$$$$$$$$$$$$
         gamma=dexp(param(5))
!$$$$$$$$$$$$$$$$$$$$$$$$$ param 6 $$$$$$$$$$$$$$$$$$$$$$$$$
         if (llay==nlayer) gamma=dexp(param(6))
         alatm=lsc%condsat(llay)
         if (layer(llay)%s(1)-layer(llay)%smin>=0.) then !check water available to drain
			   if (llay/=nlayer) then
!
! ----- 010799 Analytical solution, scaled by area
!
               beta1=.2*betavert+beta
               available=layer(llay)%s(1)-layer(llay)%smin
               if (beta1/=0) potsnew=layer(llay)%smin                         &
                                   & +(available**(-beta1)+alatm*(gamma**beta1)&
                                   & *beta1*isubstep/(24*60*60.))       &
                                   & **(1.D0/(-beta1))
               if (beta1==0) potsnew=layer(llay)%smin                         &
                                   & +available*dexp(-gamma*isubstep/(24*60*60.&
                                   & ))
               layer(llay)%latflow=layer(llay)%s(1)-potsnew !potential lateral flow

            endif
            if (llay==nlayer) then
!
! ----- 010799 Analytical solution, scaled by area
!
               beta1=.2*betavert+beta
               available=layer(llay)%s(1)-layer(llay)%smin
               if (beta1/=0) potsnew=layer(llay)%smin                         &
                                   & +(available**(-beta1)+gamma*beta1*isubstep&
                                   & /(24*60*60.))**(1.D0/(-beta1))
               if (beta1==0) potsnew=layer(llay)%smin                        &
                                   & +available*dexp(-gamma*isubstep/(24*60*60)&
                                   & )
               layer(llay)%latflow=layer(llay)%s(1)-potsnew

            endif

            sconc(1)=1.
            do kdet=2,ndet
               sconc(kdet)=layer(llay)%s(kdet)/layer(llay)%s(1)
            enddo

            do kdet=1,ndet
               if (llay<nlayer) then
!
! ----- 010799 alter scaling
!
                  layer(llay)%sdum(kdet)=layer(llay)%sdum(kdet)+layer(llay)%latflow*sconc(kdet)
                                                              
                  detred(kdet)%intf=detred(kdet)%intf+layer(llay)%latflow*sconc(kdet)
                                                      
                  layer(llay)%s(kdet)=layer(llay)%s(kdet)-layer(llay)%latflow*sconc(kdet)
                                                        

               else
!
! ----- 010799 alter scaling
!
                  layer(llay)%sdum(kdet)=layer(llay)%sdum(kdet)+layer(llay)%latflow*sconc(kdet)
                                                              !6/2/98 !6/2/98
                  detred(kdet)%basef=detred(kdet)%basef+layer(llay)%latflow*sconc(kdet)
						
						layer(llay)%s(kdet)=layer(llay)%s(kdet)-layer(llay)%latflow*sconc(kdet)
                                                        !6/2/98 !6/2/98

               endif
            enddo
         endif
      enddo
!

      do llay=1,nlayer-1
         sum2=sum2+layer(llay)%s(mbchk)
      enddo
      sum2=sum2+layer(nlayer)%s(mbchk)
      cbase=1-lscape(ilsno)%abase-lscape(ilsno)%bbase

      do kdet=1,ndet
!         sumintf(kdet)=0.D0
!         do llay=1,nlayer-1
!            sumintf(kdet)=sumintf(kdet)+layer(llay)%sdum(kdet)
!         enddo
			detred(kdet)%sumintf=detred(kdet)%intf
!         sumbasefl0(kdet)=layer(nlayer)%sdum(kdet)
			sumbasefl0(kdet)=detred(kdet)%basef
         sumbasefl0(kdet)=amin1(1.,(itime-1.)/10.)*fluxtobase(kdet) !11/05/00
         if (itime>2) sumbasefl2(kdet)=sumbasefl1(kdet)
         if (itime>1) sumbasefl1(kdet)=detred(kdet)%sumbasef
         detred(kdet)%sumbasef=lscape(ilsno)%abase*sumbasefl1(kdet)+lscape(ilsno)%bbase     &
                      & *sumbasefl2(kdet)+cbase*sumbasefl0(kdet)
      enddo
      do kdet=1,nfdet
         detred(ndet+nsdet+kdet)%sumbasef=detred(1)%sumbasef*toputfix(kdet,3)/1000
         detred(ndet+nsdet+kdet)%sumintf=detred(1)%sumintf*toputfix(kdet,2)/1000
      enddo
      sum2=sum2+layer(nlayer)%sdum(mbchk)+detred(mbchk)%sumintf+sovf
      if (imbout==1) then
         write (40,                                                     &
     &'("mass balance in latflow "/                                     &
     &                                                             "befo&
     &re, after, baseflow, interflow, overland flow,baseflow2")')
         write (40,'(6f12.6)') sum1,sum2,layer(nlayer)%sdum(mbchk),            &
                             & detred(mbchk)%sumintf,sovf,detred(mbchk)%sumbasef
      endif
! if(itime.eq.40)write(6,'("Time,LS,Baseflow",2i5,2f20.10)')itime,
! &     ilsno,sumbasef(mbchk),sumbasef(5)
!
      end
!*==STRUCTURE.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!**********************************************************************
      subroutine structure
      use i_common
!      use s_len_trim
      implicit none
!*--STRUCTURE1274
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
!
! ----- assign rate and structural parameters to soil and land use classes
!       These are constant for any
!
!
! ----- soil parameters
!
!
! ----- indices are: landuse, time of year, determinand
!
! ----- Determinands are, in order following flow:
!       1 pH
!       2 chloride
!       3 temperature
!       4 dissolved oxygen
!       5 nitrate
!       6 ammonia
!       7 particulate N
!       8 orthophosphate
!       9 particulate P
!      10 sediment
!      11 Ni
!      12 Pb
!      13 Si
!      14 Dissolved organic carbon
!      15 Particulate organic carbon
!      16 Dissolved inorganic carbon
!      17 Phytoplankton
!      18 organochlorines
!      19 pesticides
!
!
! ----- read in landscape proportions
!
 
      open (20,file=rmain(1:len_trim(rmain))//"/hrubyls")
!
      ihru=1
      do
! read(20,*)
!!!!!!! proportions read in
         read (20,*,end=100) id,                                        &
                           & (plscape(ihru,jlscape),jlscape=1,nlscape)
         ihru=ihru+1
      enddo
 100  nhru=ihru-1
      close (20)
                !finsish reading <maincat>/hrubylscape
      write (6,'("Number of hrus",i5)') nhru
! ----- select determinands to model
!
!     write(6,'("nlscape = ",i4)')nlscape
      open (60,file=rmain(1:len_trim(rmain))//'/lscape.soil.params',      &
          & status="UNKNOWN")
      read (60,'(//)')
      read (60,*) nlayer
      do ilscape=1,nlscape
			lsc=lscape(ilscape)
         read (60,*)
         read (60,*) lsc%rho,aroot,broot,aovfroot,bovfroot
         lsc%abase=aroot+broot
         lsc%bbase=-(aroot*broot)
         lsc%aovf=aovfroot+bovfroot
         lsc%bovf=-(aovfroot*bovfroot)
! write(6,'("ilscape,bovf",i5,f8.3)')ilscape,bovf(ilscape)
         do llay=1,nlayer
            read (60,*) lsc%origdepth(llay),lsc%wcs(llay),      &
                      & lsc%wcfc(llay),lsc%wcr(llay), &
                      & lsc%condsat(llay)
				lsc%origdepth(llay)=lsc%origdepth(llay)/100 !convert to m
            lsc%wcs(llay)=lsc%wcs(llay)/100 !convert from % to proportion
            lsc%wcfc(llay)=lsc%wcfc(llay)/100 !convert from % to proportion
            lsc%wcr(llay)=lsc%wcr(llay)/100 !convert from % to proportion
            lsc%condsat(llay)=lsc%condsat(llay)/100 !convert to m/day 
		   enddo
			lscape(ilscape)=lsc
		enddo
      close (60)
!
! ----- read indexes governing variables to be modelled
!
      write (6,'("read indexes of governing variables")')
      ipest=0
      do jdet=2,ndetmax
         det(jdet)%indexpest=0
         if (det(jdet)%solindex>=1) then
! write(6,'("jdet = ",i6)')jdet
            kdet=det(jdet)%solindex
            detred(kdet)%bunits=det(jdet)%aunits
            detred(kdet)%bname=det(jdet)%namel
            itype=det(jdet)%vartype
            if (itype==4) then
               ipest=ipest+1
               det(jdet)%indexpest=ipest
               pestname(ipest)=det(jdet)%name
            endif
            if (itype<4) then   !pesticide treated separately
               findexname=rmain(1:len_trim(rmain))//'/index.lscape.'//    &
                        & det(jdet)%name
               open (60,file=findexname(1:len_trim(findexname)))
               do jco=1,nlscape
                  read (60,*) mdum1,mdum2
                  lscape(mdum2)%mindex(jdet)=mdum1
               enddo
               close (60)
               nacc(jdet)=-1
               do jco=1,nlscape
                  if (lscape(jco)%mindex(jdet)>nacc(jdet)) nacc(jdet)           &
                    & =lscape(jco)%mindex(jdet)
               enddo
! write(6,*)ndetmax,jdet,nacc(jdet)
            endif
            if (itype<4) then
               findaccum=rmain(1:len_trim(rmain))//'/accum.lscape.'//     &
                       & det(jdet)%name
               open (60,file=findaccum(1:len_trim(findaccum)))
               read (60,*)
               do ilscape=1,nacc(jdet)
                  read (60,*) (accum(ilscape,jmon,jdet),jmon=1,12)
               enddo
               read (60,*)
               do ilscape=1,nacc(jdet)
                  read (60,*) (cinit(ilscape,llay,jdet),llay=1,nlayer)
! write(6,'(2i5,2f10.4)')ilscape,kdet,(cinit(ilscape,
! &                 llay,jdet),llay=1,nlayer)
               enddo
            endif
            close (60)
         endif
         if (det(jdet)%sedindex>=1) then
! write(6,'("jdet = ",i6)')jdet
            kdet=det(jdet)%sedindex
            detred(kdet+ndet)%bname=det(jdet)%namel
            detred(kdet+ndet)%bunits=det(jdet)%aunits
            itype=det(jdet)%vartype
            findexname=rmain(1:len_trim(rmain))//'/index.lscape.'//       &
                     & det(jdet)%name
            open (60,file=findexname(1:len_trim(findexname)))
            do jco=1,nlscape
               read (60,*) mdum1,mdum2
               if (mdum2>1000) mdum2=mdum2-1000
               lscape(mdum2)%mindex(jdet)=mdum1
            enddo
            close (60)
            nacc(jdet)=-1
            do jco=1,nlscape
               if (lscape(jco)%mindex(jdet)>nacc(jdet)) nacc(jdet)              &
                 & =lscape(jco)%mindex(jdet)
            enddo
! write(6,*)ndetmax,jdet,nacc(jdet)
         endif
         if (det(jdet)%fixindex>=1) then
! write(6,'("jdet = ",i6)')jdet
            kdet=det(jdet)%fixindex
            detred(kdet+ndet+nsdet)%bname=det(jdet)%namel
            detred(kdet+ndet+nsdet)%bunits=det(jdet)%aunits
            itype=det(jdet)%vartype
            findexname=rmain(1:len_trim(rmain))//'/index.lscape.'//       &
                     & det(jdet)%name
            open (60,file=findexname(1:len_trim(findexname)))
            do jco=1,nlscape
               read (60,*) mdum1,mdum2
               if (mdum2>1000) mdum2=mdum2-1000
               lscape(mdum2)%mindex(jdet)=mdum1
            enddo
            close (60)
            nacc(jdet)=-1
            do jco=1,nlscape
               if (lscape(jco)%mindex(jdet)>nacc(jdet)) nacc(jdet)              &
                 & =lscape(jco)%mindex(jdet)
            enddo
! write(6,*)ndetmax,jdet,nacc(jdet)
            findaccum=rmain(1:len_trim(rmain))//'/fixconc.lscape.'//      &
                    & det(jdet)%name
            open (60,file=findaccum(1:len_trim(findaccum)))
            read (60,*)
            do ilscape=1,nacc(jdet)
               read (60,*) (fixconc(ilscape,llay,jdet),llay=1,nlayer+1)
            enddo
            read (60,*)
            do ilscape=1,nacc(jdet)
               read (60,*) (cinit(ilscape,llay,jdet),llay=1,nlayer)
! write(6,'(2i5,2f10.4)')ilscape,kdet,(cinit(ilscape,
! &                 llay,jdet),llay=1,nlayer)
            enddo
         endif
 
!
! ----- only accumulate sediment, not associated material
!
         if (det(jdet)%sedindex>0) then
            kdet=det(jdet)%sedindex
            nacc(jdet)=-1
            do jco=1,nlscape
               if (lscape(jco)%mindex(jdet)>nacc(jdet)) nacc(jdet)              &
                 & =lscape(jco)%mindex(jdet)
            enddo
! write(6,*)ndetmax,jdet,nacc(jdet)
            findaccum=rmain(1:len_trim(rmain))//'/accum.lscape.'//        &
                    & det(jdet)%name
            open (60,file=findaccum(1:len_trim(findaccum)))
            if (kdet/=1) then   !sediment-associated determinands
               read (60,*)
! write(6,*)nacc(jdet)
               do ilscape=1,nacc(jdet)
                  read (60,*) (accum(ilscape,jmon,jdet),jmon=1,12)
! write(6,*)(accum(ilscape,jmon,jdet),jmon=1,12)
               enddo
            endif
            close (60)
         endif
      enddo
!
      npest=ipest
!
!
! ----- set solute concentrations in rainfall
!
!      open(60,file=root(1:len_trim(root))//'/cascade/params/rain.conc')
      open (60,file=rmain(1:len_trim(rmain))//'/rain.conc')
      do kdet=2,ndetmax
         read (60,*) rain%conc(kdet)
      enddo
      close (60)
!
      do ipest=1,npest
         open (60,file=rmain(1:len_trim(rmain))//'/pestinit.'//           &
             & pestname(ipest)(1:len_trim(pestname(ipest))))
         do ilscape=1,nlscape
            read (60,*) (cinit(ipest,ilscape,llay),llay=1,nlayer)
            read (60,*) (iaccump(ipest,ilscape,jmon),jmon=1,12)
         enddo
         close (60)
         open (60,file=rmain(1:len_trim(rmain))//'/pestprop.'//           &
             & pestname(ipest)(1:len_trim(pestname(ipest))))
         read (60,*)
         do llay=1,nlayer
            read (60,*)
            read (60,*) (thiscinitpads(ipest,ilscape,llay),ilscape=1,   &
                      & nlscape)                                    ! read in as concentrations m/m
         enddo
         read (60,*)
         read (60,*) akoc(ipest),pestparama(ipest),pestparamb(ipest),   &
                   & pestparamc(ipest)
         pestparama(ipest)=pestparama(ipest)*33.**(-pestparamb(ipest))
! write(6,'("pestparamA for ipest = ",i4," is ",f15.7)')ipest,
! &        pestparamA(ipest)
         close (60)
      enddo
!
      write (6,'("leaving structure")')
      end
!*==GETRAINEVAP.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!***************************************************************
      subroutine getrainevap
      use i_common
      implicit none
!*--GETRAINEVAP1542
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
!
! get rainfall for this time step
!
!
      rain%depthmm=rain_array(itime+rainskip)
      evap%potmm=evap_array(itime+evapskip)
!
      end
!*==READRAINEVAP.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!********************************************************
      subroutine readrainevap
      use f77kinds
      use i_common
      use s_julday
!      use s_len_trim
      implicit none
!*--READRAINEVAP1652
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
      real(R8KIND),intrinsic :: dmax1
      integer,intrinsic :: nint
!
!*** End of declarations rewritten by SPAG
!
!
! ----- read in rainfall and evaporation data. The coordinates ixm and iym refer to a
!       10km grid
!
!
!
      write(6,'("Read rainfall and evaporation")')
		open (60,file=rmain(1:len_trim(rmain))//'/rain.data',iostat=iostat)
      open (61,file=rmain(1:len_trim(rmain))//'/evap.data',iostat=iostat)
      do jtype=1,2
         infile=60
         if (jtype==2) then
				infile=61
				read(infile,'(///////)')
			endif
         ico=0
         do
            ico=ico+1
            read (infile,'(2(i2,1x),i4,1x,2(i2,1x),i2,10f12.0)',end=50) &
                & imon,iday,iyr,ihr,imin,isec,val
            write (6,*) imon,iday,iyr,ihr,imin,isec,val
            if (jtype==1) then
               rain_array(ico)=dmax1(val,0.)
! write(6,*)val
               if (ico==1) arstart=julday(imon,iday,iyr)*1.D0+ihr/24.+  &
                                 & imin/(24*60.)+isec/(24*60*60.)
            else
               evap_array(ico)=dmax1(val,0.)
               if (ico==1) aestart=julday(imon,iday,iyr)*1.D0+ihr/24.+  &
                                 & imin/(24*60.)+isec/(24*60*60.)
            endif
         enddo
!
 50      if (jtype==1) then
            arend=julday(imon,iday,iyr)*1.D0+ihr/24.+imin/(24*60.)      &
                & +isec/(24*60*60.)
         else
            aeend=julday(imon,iday,iyr)*1.D0+ihr/24.+imin/(24*60.)      &
                & +isec/(24*60*60.)
         endif
      enddo
!
      rainskip=nint((astart-arstart)*24*60*60/isubstep)
      evapskip=nint((astart-aestart)*24*60*60/isubstep)
      if (rainskip<0) then
         write (6,                                                      &
     &'(                                                                &
     &                                                             "Term&
     &ination - no rainfall data",                                      &
     &                                                           " at st&
     &art of model period")')
         stop
      endif
      if (evapskip<0) then
         write (6,                                                      &
     &'(                                                                &
     &                                                             "Term&
     &ination - no evaporation",                                        &
     &                                                           " data &
     &at start of model period")')
         stop
      endif
      if (aend>arend) then
         write (6,                                                      &
     &'(                                                                &
     &                                                             "Term&
     &ination - no rainfall data at end of model period")')
         stop
      endif
      if (aend>aeend) then
         write (6,                                                      &
     &'(                                                                &
     &                                                             "Term&
     &ination - no evaporation",                                        &
     &                                                           " data &
     &at end of model period")')
         stop
      endif
!
      close (60)
      close (61)
!
      end
!*==SEDGENNEW.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!*****************************************************************
      subroutine sedgennew
      use i_common
!      use s_len_trim
      implicit none
!*--SEDGENNEW1811
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
!
! ----- generate sediment load using PSN + DMC paper (commented out is
!       old code
!
!
! compute sediment concentration in relationship sedconc=aksed*flow
!
      open (10,file=rmain(1:len_trim(rmain))//'/sedparams',iostat=iostat)
      read (10,*)
      do i=1,nlscape
         read (10,*) id,lscape(i)%aksed
         write (6,*) id,lscape(i)%aksed,nlscape
      enddo
      close (10)
      end
!*==READTEMP.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!*****************************************************************
      subroutine readtemp
      use i_common
!      use s_len_trim
      implicit none
!*--READTEMP1844
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
!
! ----- read temperature data
!
!
      open (60,file=rmain(1:len_trim(rmain))//'/dailytemp')
!
      i=0
      do
         i=i+1
         read (60,*,end=100) temp(i)
         if (temp(i)<-90.) temp(i)=temp(i-1)
      enddo
!
 100  close (60)
!
      end
!*==LSQFUN.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!*************************************************************************
!
! ----- optimisation routine
!
      subroutine lsqfun(mopt,nopt,paramopt,fveccopt)
      use f77kinds
      use i_common
      use s_caldat
      use s_datey
      use s_getrainevap
      use s_stateupd
      use s_structure
      implicit none
!*--LSQFUN1888
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      integer :: mopt,nopt
      real(R8KIND),dimension(mopt) :: fveccopt
      real(R8KIND),dimension(nopt) :: paramopt
      intent (in) mopt,nopt,paramopt
      intent (inout) fveccopt
!
! Local variables
!
      real(R8KIND),intrinsic :: dexp,dlog,dmax1,dsqrt
      real(R8KIND) :: fss,fvmax,tt1,tt2,tt3
      integer,intrinsic :: int
      integer :: itt1,itt2,itt3
      real,intrinsic :: sngl
      logical :: validok
      character(10) :: validdate
!
!*** End of declarations rewritten by SPAG
!
!
! write(6,'("Entered lsqfun")')
!
! ----- initialise
!
      ncall=ncall+1
      if (istopping==1) iflagopt=-1
      do iparam=1,nparam
         if (idoparam(iparam)>0.and.istopping==0) param(iparam)         &
           & =paramopt(idoparam(iparam))
      enddo
      if (isopt==1) then
         write (*,'("Parameter estimates iteration",i5)') ncall
         write (*,'(4f20.10)') (dexp(param(iparam)),iparam=1,nparam)
      endif
!
! write(6,'("Calling structure")')
      if (ncall==1) call structure
! write(6,'("Out of structure")')
!
! ----- adjust depths of soil by an adjustable parameter
!
      do ilscape=1,nlscape
         do llay=1,nlayer
            if (llay/=nlayer) lscape(ilscape)%depth(llay)                       &
              & =lscape(ilscape)%origdepth(llay)*dexp(param(7))
            if (llay==nlayer) lscape(ilscape)%depth(llay)                          &
              & =lscape(ilscape)%origdepth(llay)*dexp(param(9))
         enddo
      enddo
      do ico=1,nsubcat
         subcat(ico)%area=0
      enddo
      iflag=0
      do i=1,30
         totlusev(i)=0.
      enddo
      jsoil=1
      do jmon=1,12
         ileap(jmon)=0
      enddo
!
! ----- initialise mass balances
!
      do ilsno=1,nlscape
         lscape(ilsno)%ambalstart=0.
         lscape(ilsno)%ambalend=0.
         lscape(ilsno)%ambalrain=0.
         lscape(ilsno)%ambalroff=0.
      enddo
      do itime=1,ntime
         do jdet=1,ndet+nsdet+nfdet
            do ksubcat=1,nsubcat
               subcat(ksubcat)%totsum(itime,jdet)=0.
               subcat(ksubcat)%totsumovf(itime,jdet)=0
               subcat(ksubcat)%totsumintf(itime,jdet)=0
               subcat(ksubcat)%totsumbasef(itime,jdet)=0
            enddo
         enddo
      enddo
      do ilscape=1,nlscape
         lscape(ilscape)%sumrain=0.
         lscape(ilscape)%sumpe=0.
         lscape(ilscape)%sumae=0.
      enddo
      fvmax=0.
!
      do idet=1,ndet
         sumbasefl2(idet)=0.
         sumbasefl1(idet)=0.
      enddo
!
! ----- compute hru areas
!
      open (10,file=rmain(1:len_trim(rmain))//'/hruarea.data')
      ihru=1
      do
         read (10,*,end=100) ihruno,jhruarea(ihru)
         do ksubcat=1,nsubcat
            subcat(ksubcat)%area=subcat(ksubcat)%area+jhruarea(ihru)          &
                            & *subcat(ksubcat)%hruflag(ihru)
         enddo
         ihru=ihru+1
      enddo
 100  close (10)
 
! write(6,'("Starting loop on landscapes")')
!
! ----- LOOP ON LANDSCAPE TYPES ----------------------------------
 
      do ilsno=1,nlscape
			lsc=lscape(ilsno)
         write (6,'("ilsno = ",i4)') ilsno
!
!
 
! ----- constraints and initial values for non-sediment determinands
!
         do llay=1,nlayer
            do jdet=2,ndetmax
               kdet=det(jdet)%solindex
               ipest=det(jdet)%indexpest
               if (kdet>=1) then
                  mm=lsc%mindex(jdet)
                  layer(llay)%thiscinit(kdet)=cinit(mm,llay,jdet)
               endif
            enddo
         enddo
!
! ----- initialise soil capacities and contents for this landscape type
!
         do llay=1,nlayer
            layer(llay)%smax=lsc%depth(llay)*lsc%wcs(llay) !m**3 per m**2
            layer(llay)%scfc=lsc%depth(llay)*lsc%wcfc(llay)
            layer(llay)%smin=lsc%depth(llay)*lsc%wcr(llay) !m**3 per m**2
! sinit(llay)=smax(llay) ! initial water content set to saturation
         enddo
         do idet=1,ndet
            sumbasefl2(idet)=0.
            sumbasefl1(idet)=0.
            ovfl2(idet)=0.
            ovfl1(idet)=0.
			enddo
!
! ----- initialise water content by layer
!
         do llay=1,nlayer
! sinit(llay)=(scfc(llay)+smin(llay))/2
            layer(llay)%sinit=layer(llay)%smin
         enddo
!
         do llay=1,nlayer
            layer(llay)%s(1)=layer(llay)%sinit !m**3 per m**2 since sinit is m**3/gridsq
            lsc%ambalstart=lsc%ambalstart+layer(llay)%s(1)
            do ipest=1,npest

               sadspest(ipest,llay)=thiscinitpads(ipest,ilsno,llay)     &
                                  & *lsc%rho/1000.          !kg
            enddo
!
! ----- initialise mass of each chemical determinand by layer
!
            do kdet=2,ndet
               layer(llay)%s(kdet)=layer(llay)%thiscinit(kdet)*layer(llay)%sinit/1000 !kg
            enddo
         enddo
! enddo
! ----- LOOP ON TIME STEP* ******************************************************************LOOP ON TIME STEP
!
         aday=astart-isubstep/(24*60*60.)
         do itime=1,ntime
!
! ----- reset time counters
!
!				write(6,'("ilsno and itime ",2i5)')ilsno,itime			
            aday=aday+isubstep/(24*60*60.)
            jjday=int(aday)
            call caldat(jjday,imonth,iday,iyear)
            tt1=(aday-jjday)*24
            itt1=int(tt1)
            tt2=(tt1-itt1)*60
            itt2=int(tt2)
            tt3=(tt2-itt2)*60
            itt3=int(tt3)
            iyearv(itime)=iyear
            imonthv(itime)=imonth
            idayv(itime)=iday
            ihrv(itime)=itt1
            iminv(itime)=itt2
            isecv(itime)=itt3
!
! ----- get driving variables for this landscape class and time step
!
            call getrainevap
!
            thistemp=25
!
            lsc%sumrain=lsc%sumrain+rain%depthmm
            lsc%sumpe=lsc%sumpe+evap%potmm !mm x no of grid squares
 
! ----- find meteorological sources for this time step
!
            do jdet=2,ndetmax
               kdet=det(jdet)%solindex
               if (kdet>=1) then
                  rainsource(kdet)=rain%depthmm*rain%conc(jdet)/(1000.*1000.)
                                                                  !kg/m^2 from mg/l
                  rain%concred(kdet)=rain%conc(jdet)
               endif
            enddo
!
! ----- transfer water and determinands this time step
!

				call stateupd

!
! ----- store soil moisture deficit information by time step and landscape
!
            do llay=1,nlayer-1
               dsmd=1000*(layer(llay)%scfc-layer(llay)%s(1))
               lsc%smd(itime,llay)=dmax1(0.D0,dsmd)
            enddo
!
! ----- store fluxes by timestep and landscape
!
            dsmd=1000
            do jdet=1,ndet
               if (jdet==-1) write (6,'("itime,ilsno,jdet,sumovf,sumintf,sumbasef",3i5,3f10.4)') &
						&itime,ilsno,jdet,detred(jdet)%sumovf*dsmd,detred(jdet)%sumintf*dsmd,&
                & detred(jdet)%sumbasef*dsmd                      !mm per time step
               lsc%sumf(itime,jdet)=detred(jdet)%sumovf+detred(jdet)%sumintf        &
                                    & +detred(jdet)%sumbasef             !M**3/day output only
            enddo
            do kdet=1,nfdet
               jdet=ndet+nsdet+kdet
               lsc%sumf(itime,jdet)=detred(jdet)%sumovf+detred(jdet)%sumintf        &
                                    & +detred(jdet)%sumbasef             !M**3/day output only
            enddo
            lsc%ambalroff=lsc%ambalroff+lsc%sumf(itime,1)
!
! ----- sediment concentration (mg/l) = aksed * flow (mm/day)
!       divide 1000
!                              (kg/m**3) = aksed * flow (m/day)
            do jdet=1,nsdet
               if (jdet==1) then
 
                  sdconc=lsc%aksed*lsc%sumf(itime,1)
                                           !nb conc in kg/m**3
                  lsc%sumf(itime,jdet+ndet)=sdconc*lsc%sumf(itime,1)
               else
                  lsc%sumf(itime,jdet+ndet)=lsc%sumf(itime,1+ndet) &
                   & *toputsed(jdet)
               endif
            enddo
         enddo                        ! end of timestep loop
!
! ----- END OF TIME STEP LOOP**************************************************END OF TIME STEP LOOP
!
!			water content in each landscape at end of run         
			do llay=1,nlayer
            lsc%ambalend=lsc%ambalend+layer(llay)%s(1) !water content of soil at end of run 
         enddo
			lscape(ilsno)=lsc
!
      enddo
!
! ----- END OF LOOP ON LANDSCAPE TYPES ****************************************************END OF LOOP ON REGIONS
!
      write (6,'("Out of model loop")')
!
! ----- Distribute landscape type outputs around HRUs
!
      do itime=1,ntime
         do ico=1,nsubcat
            do kdet=1,ndetmax
               subcat(ksubcat)%totsum(itime,kdet)=0
            enddo
         enddo
      enddo
!
! ----- Compte flux out at each time step by summing contributions 
!			from each landscape in each HRU
!
      do ihru=1,nhru
         do ico=1,nsubcat
				subc=subcat(ico)
            if (subc%hruflag(ihru)==1) then
               do itime=1,ntime
                  do kdet=1,ndet+nsdet+nfdet
                     totsumf(itime,kdet)=0.
                     do jlscape=1,nlscape
                        totsumf(itime,kdet)=totsumf(itime,kdet)         &
                         & +plscape(ihru,jlscape)                       &
                         & *lscape(jlscape)%sumf(itime,kdet)
							enddo
                     totsumf(itime,kdet)=totsumf(itime,kdet)            &
                      & *jhruarea(ihru)*50*50  !total loss
                     subc%totsum(itime,kdet)                         &
                      & =subc%totsum(itime,kdet)+totsumf(itime,kdet)
                  enddo
               enddo
               do idh=1,4
                  if (ihru>=10**(idh-1).and.ihru<10**idh) kch=idh
               enddo
               write (chireg,'(i6)') ihru
               mch=6-kch+1
               fout=rmain(1:len_trim(rmain))//'/'//subc%name        &
                  & (1:len_trim(subc%name))//'/output/HRU'//        &
                  & chireg(mch:6)
! write(6,'(a80)')fout(1:len_trim(fout))
               open (41,file=fout(1:len_trim(fout)),status="UNKNOWN")
               fout=rmain(1:len_trim(rmain))//'/'//subc%name        &
                  & (1:len_trim(subcat(ico)%name))//'/smd.sims'
               open (42,file=fout(1:len_trim(fout)))
               write (41,                                               &
     &'("Start and end date",3x,i2.2,"/",i2.2                           &
     &                                                                ,"&
     &/",i4,1x,i2.2,":",i2.2,":",i2.2,3x,i2.2,"/"                       &
     &                                                              ,i2.&
     &2,"/",i4,1x,i2.2,":",i2.2,":",i2.2)') smon,sday,syea,shr,smina,   &
    & ssec,emon,eday,eyea,ehr,emin,esec
               write (41,'("Step length in seconds",i8)') isubstep
               write (41,'("Number of data",i8)') ntime
               write (41,                                               &
     &'("Number of determinands (including flow) ",                     &
     &                                                                i5&
     &)') nsdet+ndet+nfdet
               write (41,'("HRU number",i8)') ihru
               write (41,'(" 1. Discharge m**3/sec")')
               do jdet=2,ndet+nsdet+nfdet
                  write (41,'(i2,". ",a60)') jdet,detred(jdet)%bname
               enddo
               write (41,*)
               write (41,'("Mo Da Year Hr Mi Se",20(4x,i2,4x))')        &
                    & (jdet,jdet=1,ndet+nsdet+nfdet)
               write (41,*)
!
               do jdet=2,ndet+nsdet+nfdet
! write(41,'(a60)')bname(jdet)
                  do it=1,ntime
                     if (totsumf(it,1)/isubstep>1.D-4) then
                        dds(it,jdet)=totsumf(it,jdet)*1000/totsumf(it,1)
                        if (detred(jdet)%bunits=="ug") dds(it,jdet)            &
                          & =dds(it,jdet)*1000.
                     else
                        dds(it,jdet)=0.D0
                     endif
                  enddo
               enddo
               do it=1,ntime
                  write (41,                                            &
     &'(2(i2.2,"/"),i4,1x,i2.2,":",i2.2,":"                             &
     &                                                            ,i2.2,&
     &20f10.4)') imonthv(it),idayv(it),iyearv(it),ihrv(it),iminv(it),   &
               & isecv(it),totsumf(it,1)/isubstep,                      &
               & (dds(it,jdet),jdet=2,ndet+nsdet+nfdet)
               enddo
               write (42,'("HRU number",i5)') ihru
               close (41)
               close (42)
            endif
				subcat(ico)=subc
         enddo
      enddo
!
      ttotsum=0.
!
! ----- get actual flows for checking
!
!     write(6,'("nsubcat = ",i5)')nsubcat
!     stop
      write (6,'("Read validation data")')
      if (ncall==1) then
         do it=1,2922
            do jdet=1,23
               do ksubcat=1,nsubcat
                  subcat(ksubcat)%valid(it,jdet)=-1
               enddo
            enddo
         enddo
         do ksubcat=1,nsubcat
				subc=subcat(ksubcat)
            rsub=rmain(1:len_trim(rmain))//'/'//subc%name           &
               & (1:len_trim(subc%name))
            validok=.false.
            do while (.not.validok)
               open (80,file=rsub(1:len_trim(rsub))//'/model.valid',      &
                   & status="OLD",err=110)
 110           write (6,'("No validation data - stop program")')
               stop
 
            enddo
            stop
            read (80,'(10x,a8)') validdate
            read (validdate,'(i4,i2,i2)') ivy,ivm,ivd
            backspace (80)
            ifl=1
            call datey(ndays,sday,smon,syea,ivd,ivm,ivy,ifl)
            if (ndays<0) then
               do i=1,-ndays
                  read (80,*)
               enddo
               it=0
               subc%sflchk=0.
               do
                  it=it+1
                  read (80,*,end=120) id1,id2v(it),                     &
                                    & (subc%valid(it,j),j=1,23)
                  subc%valid(it,1)=subc%valid(it,1)-ps_correction
                  subc%sflchk=subc%sflchk+subc%valid(it,1)
                  if (it==1) isvdate=id1
               enddo
 120           subc%nflchk=it-1
            else
               it=ndays
               subc%sflchk=0.
               do while (validok)
                  it=it+1
                  read (80,*,end=130) id1,id2v(it),(subc%valid(it,j),j=1,23)
                  subc%valid(it,1)=subc%valid(it,1)-ps_correction
                  subc%sflchk=subc%sflchk+subc%valid(it,1)
                  if (it==ndays+1) isvdate=id1
               enddo
 130           subcat(ksubcat)%nflchk=it-1-ndays
            endif
            close (80)
            subc%sflchk=subc%sflchk/subc%nflchk
				subcat(ksubcat)=subc
         enddo
!
			do ksubcat=1,nsubcat
				subc=subcat(ksubcat)
				do it=1,subc%nflchk
					subc%validr(it,1)=subcat(ksubcat)%valid(it,1)
					do jdet=2,ndetmax
						kdet=det(jdet)%solindex
						if (kdet>=1) then
                     subc%validr(it,kdet)=subc%valid(it,jdet)
						endif
						kdet=det(jdet)%sedindex
						if (kdet>=1) then
                     subc%validr(it,kdet+ndet)=subc%valid(it,jdet)
						endif
						kdet=det(jdet)%fixindex
						if (kdet>=1) then
                     subc%validr(it,kdet+ndet+nsdet)=subc%valid(it,jdet)
		            endif
					enddo
				enddo
				subcat(ksubcat)=subc
         enddo
      endif
!
! ----- Output simulations for complete catchment as total flow
!       component and as surface, inter and baseflow separately
!
      if (istopping==1) then
         open (83,file=rmain(1:len_trim(rmain))//'/cdheader.file')
         do
            read (83,'(a)',end=150) ach
            jeco=0
            do kl=1,80
               if (ach(kl:kl)/=' ') jeco=kl
            enddo
            do ksubcat=1,nsubcat
               rsub=rmain(1:len_trim(rmain))//'/'//subcat(ksubcat)%name        &
                  & (1:len_trim(subcat(ksubcat)%name))
               open (82,file=rsub(1:len_trim(rsub))//'model3.sims.cd')
               write (82,'(a)') ach(1:jeco)
               close (82)
            enddo
         enddo
 150     close (83)
      endif
!
! ----- Compute optimisation criterion
!
      fss=0.D0
!
      iorgd=1
      iorgm=1
      iorgy=1985
      ifl=1
      call datey(ndd,iorgd,iorgm,iorgy,sday,smon,syea,ifl)
      igap=ndd
      do itime=1,ntime
         do kdet=1,ndet
            fveccopt(itime+(kdet-1)*ntime)=0.
         enddo
      enddo
      do ksubcat=1,nsubcat
			subc=subcat(ksubcat)
         do itime=1,ntime
            if (subc%totsum(itime,1)<1.E-8) then
               do jdet=ndet+1,ndet+nsdet
                  toto(jdet)=-1.
               enddo
            else
               do jdet=ndet+1,ndet+nsdet
                  toto(jdet)=subc%totsum(itime,jdet)                 &
                           & *1000/subc%totsum(itime,1)
               enddo
            endif
 
!
! ----- output
!
            if (istopping==1) then
               rsub=rmain(1:len_trim(rmain))//'/'//subc%name        &
                  & (1:len_trim(subc%name))
               open (80,file=rsub(1:len_trim(rsub))//'/model3.sims.full')
               open (91,file=rsub(1:len_trim(rsub))//'/model3.sims.short')
               open (82,file=rsub(1:len_trim(rsub))//'model3.sims.cd')
               write (80,'(2i8)') itime,itime+2446066+igap
               write (80,'(4f10.4)') subc%totsum(itime,1)/isubstep,  &
                                   & subc%totsumovf(itime,1)         &
                                   & /isubstep,                         &
                                   & subc%totsumintf(itime,1)        &
                                   & /isubstep,                         &
                                   & subc%totsumbasef(itime,1)       &
                                   & /isubstep
               do jdet=2,ndet
                  write (80,'(4f10.4)') subc%totsum(itime,jdet)       &
                                      & *1000/subc%totsum(itime,1),  &
                                      & toto(jdet),toti(jdet),totb(jdet)
               enddo
               do jdet=ndet+1,ndet+nsdet
                  write (80,'(2f10.4)')subc%totsum(itime,jdet)       &
                                      & *1000/subc%totsum(itime,1),  &
                                      & toto(jdet)
               enddo
 
            endif
!
! ----- Compute contribution to loss function
!
            if (isopt==1) then
               do kdet=1,ndet
                  fvecdum(kdet)=0.
               enddo
               if (ioptwt==1) then
                  if (itime>25.and.subc%validr(itime,1)>0) fvecdum(1) &
                    & =(1000*isubstep*subc%validr(itime,1)           &
                    & /(50*50*subc%area))                       &
                    & -(1000*subc%totsum(itime,1)                    &
                    & /(50*50*subc%area))                 !units are mm/day !6/2/98 !6/2/98
               elseif (ioptwt==-1) then
                  if (itime>25.and.subc%validr(itime,1)>0) fvecdum(1) &
                    & =dlog(1000*isubstep*subc%validr(itime,1)        &
                    & /(50*50*subc%area))                       &
                    & -dlog(1000*subc%totsum(itime,1)                &
                    & /(50*50*subc%area))                 !units are mm/day !6/2/98 !6/2/98
               else
                  write (*,'("Check value of ioptwt should be -1 or 1")')
                  stop
               endif
               do jdet=2,ndet
                  if (subc%validr(itime,jdet)>0.and.                  &
                    & subc%totsum(itime,1)>1.D-6) fvecdum(jdet)      &
                    & =optsolwt(jdet)                                      &
                    & *(subc%validr(itime,jdet)-subc%totsum(itime,jdet)       &
                                      & *1000/subc%totsum(itime,1))
               enddo
               fveccopt(itime)=fveccopt(itime)+fvecdum(1)*fvecdum(1)
               do kdet=2,ndet
                  fveccopt(itime+(kdet-1)*ntime)                        &
                   & =fveccopt(itime+(kdet-1)*ntime)+fvecdum(kdet)      &
                   & *fvecdum(kdet)
               enddo
 
               if (ksubcat==nsubcat) then
                  do kdet=1,ndet
                     fveccopt(itime+(kdet-1)*ntime)                     &
                      & =dsqrt(fveccopt(itime+(kdet-1)*ntime))
                     fss=fss+fveccopt(itime+(kdet-1)*ntime)             &
                       & *fveccopt(itime+(kdet-1)*ntime)
                  enddo
               endif
               if (fveccopt(itime)*fveccopt(itime)>fvmax) then
                  ifvmax=itime
                  fvmax=fveccopt(itime)*fveccopt(itime)
               endif
            endif
            if (istopping==1) then
               ipr=itime+2446066+igap
               call caldat(ipr,icmon,icday,icyr)
               write (91,'(i8,2i3,i5,20f12.5)') ipr,icday,icmon,icyr,   &
                    & subc%totsum(itime,1)/isubstep,                 &
                    & subc%validr(itime,1),                           &
                    & (subc%totsum(itime,j)       &
                    & *1000/subc%totsum(itime,1),subc%validr(itime,j),j=2,              &
                    & ndet+nsdet+nfdet)
            endif
!
			subcat(ksubcat)=subc
         enddo
         if (isopt==1.and.ksubcat==nsubcat) then
            write (*,'("Optimisation criterion this iteration")')
            write (*,'(f20.10)') fss
         endif
         close (91)
         close (80)
         close (82)
      enddo
 
      if (istopping==1) then
         open (80,file=rmain(1:len_trim(rmain))//'/mbal.chk')
         do ilscape=1,nlscape
				lsc=lscape(ilscape)
            write (80,'(//"**************************************")')
            write (80,'(/"Water budget for landscape type ",i4/)')      &
                 & ilscape
            write (80,'("Initial water content")')
            write (80,'(f30.5)') lsc%ambalstart
            write (80,'("Effective rainfall over period")')
            write (80,'(f30.5)') lsc%ambalrain
            write (80,'("Simulated runoff over period")')
            write (80,'(f30.5)') lsc%ambalroff
            write (80,                                                  &
                  &'("Initial water content+rainfall-simulated flow")')
            write (80,'(f30.5)') lsc%ambalstart+lsc%ambalrain &
                               & -lsc%ambalroff
            write (80,'("Final simulated water content")')
            write (80,'(f30.5)') lsc%ambalend
            write (80,'("Annual total rainfall (mm)",f10.5)')           &
                 & lsc%sumrain*365./ntime
            write (80,'("Annual total pe. (mm)",f10.5)') lsc%sumpe &
                 & *365./ntime
            write (80,'("Annual total ae. (mm)",f10.5)') lsc%sumae &
                 & *365./ntime
         enddo
         close (80)
      endif
!
      close (10)
      close (11)
      close (20)
      close (21)
      close (22)
      close (23)
!
      end
!*==LSQMON.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!
!*****************************************************************
!
      subroutine lsqmon(mopt,nopt,paramopt,fveccopt,fjaccopt,ljcoptt,   &
                      & sopt,igrade,niteropt,nfopt,iwopt,liwoptt,wopt,  &
                      & lwoptt)
      use f77kinds
      use i_common
      implicit none
!*--LSQMON2611
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      integer :: igrade,liwoptt,ljcoptt,lwoptt,mopt,nfopt,niteropt,nopt
      real(R8KIND),dimension(ljcoptt,nopt) :: fjaccopt
      real(R8KIND),dimension(mopt) :: fveccopt
      integer,dimension(liwoptt) :: iwopt
      real(R8KIND),dimension(nopt) :: paramopt,sopt
      real(R8KIND),dimension(lwoptt) :: wopt
      intent (in) liwoptt,ljcoptt,lwoptt,mopt,nfopt,nopt
!
! Local variables
!
      character(1) :: ag
!
!*** End of declarations rewritten by SPAG
!
!
      write (*,'("number of calls to lsqmon so far",2x,i5)') nfopt
      if (aintopt=='y') then
         write (*,'("continue optimization? y/n")')
         read (5,'(a1)') ag
      else
         ag='y'
      endif
      if (ag=='y') then
         write (*,'("Continuing optimisation")')
         istopping=0
      else
         istopping=1
         write (*,'("Halting optimisation and outputting results")')
      endif
!
      end
!*==DATEY.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
!*****************************************************************
 
      subroutine datey(nn,ids,ims,iys,idf,imf,iyf,ifl)
      implicit none
!*--DATEY2658
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      integer :: idf,ids,ifl,imf,ims,iyf,iys,nn
      intent (in) ids,ifl,ims,iys
      intent (inout) idf,imf,iyf,nn
!
! Local variables
!
      real :: dd
      integer :: i,k,m,n1,n2
!
!*** End of declarations rewritten by SPAG
!
!
! ----- returns number of days in date range (C-GW's external code)
!
      i=iys
      m=ims-3
      if (m<0) then
         m=m+12
         i=i-1
      endif
      k=i/100-16
      i=mod(i,100)
      n1=36524.*k+365.*i+(k/4+i/4+(153*m+2)/5+ids)-55714.
      if (ifl==2) then
         dd=n1+nn+445711./8.
         k=(4.*dd)/146097.
         dd=dd-36524.*k-k/4
         i=(4.*dd)/1461.
         dd=dd-365.*i-(i/4)
         i=100*(k+16)+i
         k=5.*dd-1.875
         m=k/153
         k=(k-153*m+5)/5
         if (m>9) then
            imf=m-9
            i=i+1
         else
            imf=m+3
         endif
         iyf=i
         idf=k
      else
         i=iyf
         m=imf-3
         if (m<0) then
            m=m+12
            i=i-1
         endif
         k=i/100-16
         i=mod(i,100)
         n2=36524.*k+365.*i+(k/4+i/4+(153*m+2)/5+idf)-55714.
         nn=n2-n1
         return
      endif
 
      end
!*==JULDAY.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
! ***************************************************************
!
      function julday(mm,id,iyyy)
      implicit none
!*--JULDAY2730
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
      integer,parameter :: igrega=15+31*(10+12*1582)
!
! Dummy arguments
!
      integer :: id,iyyy,mm
      integer :: julday
      intent (in) id,iyyy,mm
!
! Local variables
!
      integer,intrinsic :: int
      integer :: ja,jm,jy
!
!*** End of declarations rewritten by SPAG
!
!
! ----- converts month/day/year to Julian day
!
      jy=iyyy
      if (jy==0) pause 'julday: there is no year zero'
      if (jy<0) jy=jy+1
      if (mm>2) then
         jm=mm+1
      else
         jy=jy-1
         jm=mm+13
      endif
      julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
      if (id+31*(mm+12*iyyy)>=igrega) then
         ja=int(0.01*jy)
         julday=julday+2-ja+int(0.25*ja)
      endif
      end
!*==CALDAT.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
! (C) Copr. 1986-92 Numerical Recipes Software.
!
!********************************************************************************************
!
      subroutine caldat(julian,mm,id,iyyy)
      implicit none
!*--CALDAT2781
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
      integer,parameter :: igreg=2299161
!
! Dummy arguments
!
      integer :: id,iyyy,julian,mm
      intent (in) julian
      intent (out) id
      intent (inout) iyyy,mm
!
! Local variables
!
      integer,intrinsic :: int
      integer :: ja,jalpha,jb,jc,jd,je
!
!*** End of declarations rewritten by SPAG
!
!.....converts Julian day to month/day/year
      if (julian>=igreg) then
         jalpha=int(((julian-1867216)-0.25)/36524.25)
         ja=julian+1+jalpha-int(0.25*jalpha)
      else
         ja=julian
      endif
      jb=ja+1524
      jc=int(6680.+((jb-2439870)-122.1)/365.25)
      jd=365*jc+int(0.25*jc)
      je=int((jb-jd)/30.6001)
      id=jb-jd-int(30.6001*je)
      mm=je-1
      if (mm>12) mm=mm-12
      iyyy=jc-4715
      if (mm>2) iyyy=iyyy-1
      if (iyyy<=0) iyyy=iyyy-1
      end
!*==INITIALISE.s.f  processed by SPAG 6.55Dc at 15:04 on  6 Mar 2006
!*------------------   SPAG Configuration Options   --------------------
!*--0323,76 001101,-1 000000101011332110000002000020210201010,136 11 --
!*--1100000000012000000000000000,100,50,20,10 52,99000 1200000000000 --
!*--000000000000000000,72,72 73,42,38,33 30011012110000100000000     --
!*----------------------------------------------------------------------
! (C) Copr. 1986-92 Numerical Recipes Software )(.
!****************************************************************************************
!
      subroutine initialise(startt,endt)
      use i_common
      use s_julday
      implicit none
!*--INITIALISE2835
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      character(14) :: endt,startt
      intent (inout) endt,startt
!
! Local variables
!
      integer,intrinsic :: nint
!
!*** End of declarations rewritten by SPAG
!
!
      open (10,file="filesource")
! Read name of the main catchment directory
      read (10,'(a180)') rmain
      read (10,'(a180)') rout
		idailyevap=0
		write(6,*)'test1'
! Subcatchment names assumed held in <main catchment>/model.sites
      open (11,file=rmain(1:len_trim(rmain))//'/model.sites')
        write(6,*)'test2'!
! ----- read subcatchment names
!
      ico=1
      do
         read (11,'(a60)',end=100) subcat(ico)%name
         write (6,'(a60)') subcat(ico)%name
         ico=ico+1
      enddo
 100  nsubcat=ico-1
      close (11)
!
! ----- Optimisation parameters
!
!     Read point source correction (m**3/s)
      read (10,*) ps_correction
! Is this an optimisation run? y/n
      read (10,'(a1)') qopt
      isopt=0
      if (qopt=='y') then
         isopt=1
! Read maximum step length in optimisation
         read (10,*) stepmx
! Read eta value
         read (10,*) eta
! Interactive optimisation? y/n")')
         read (10,'(a1)') aintopt
      endif
!c ----- read information for this run
! ----- Changed by CWang. All information read in from keyboard (18/2/1999)
 
!     Enter the time period for this run:'
!     Read start date (yyyymmdd)
      read (10,'(a14)') startt
		write(6,'(a14)') startt
! Read end date (yyyymmdd)
      read (10,'(a14)') endt
		write(6,'(a14)') startt
      read (10,*) isubstep !seconds
!      write (6,'(a14,1x,a14)') startt,endt
		write (6,'(a14)') startt
!
      read (startt,'(i4,5i2)') syea,smon,sday,shr,smina,ssec
      write (6,'(i4,5(1x,i2))') syea,smon,sday,shr,smina,ssec
      astart=julday(smon,sday,syea)*1.D0+shr/24.+smina/(24*60.)         &
           & +ssec/(24*60*60.)
      read (endt,'(i4,5i2)') eyea,emon,eday,ehr,emin,esec
      write (6,'(i4,5(1x,i2))') eyea,emon,eday,ehr,emin,esec
      aend=julday(emon,eday,eyea)*1.D0+ehr/24.+emin/(24*60.)            &
         & +esec/(24*60*60.)
      ntime=nint((aend-astart)*24*60*60/isubstep)+1
!
      write (6,'("Number of time steps modelled: ",i8)') ntime
      write (6,*)
! Is HRU output to be generated (y/n)?")')
      read (10,'(a1)') answer
! read(5,'(a1)')answer
      if ((answer=='y').or.(answer=='Y')) then
         ioutfbyhru=1
      else
         ioutfbyhru=0
      endif
! Is a detailed mass balance check on 1 determinand required? y/n")')
      read (10,'(a1)') aimbout
      imbout=0
      if (aimbout=='y') then
         imbout=1
! Read index number of determinand
         read (10,*) mbchk
      endif
      if (imbout==0) then
         write (6,'("No detailed Mass balance check")')
      elseif (imbout==1) then
         write (6,'("Detailed mass balance check on determinand",i3)')  &
              & mbchk
         write (6,'("Output to <subcat>/modelout.chk")')
      endif
!
      close (10)
      write (6,'("Finished initialisation")')
!
      end
!
!
      subroutine compute_actevap
      use i_common
		implicit none

! ----- if net water loss (by evapotranspiration) exceeds available
!       water content of top layer, take from next layer. If rain
!       exceeds
!       available water capacity of top layer, generate overland flow
!
      do llay=1,nlayer
         layer(llay)%sav=dmin1(1.D0,(layer(llay)%snew(1)-layer(llay)%smin)/ &
			(layer(llay)%smax-layer(llay)%smin))
      enddo
!
! ----- evaporation from layer 1: sav1 is proportional saturation
!
      do llay=1,nlayer
         layer(llay)%ddepth=0.D0
      enddo
      dg=0.D0
      do llay=1,nlayer
         dg=dg+lsc%depth(llay)
         if (dg<1.D0) then
            layer(llay)%ddepth=lsc%depth(llay)
         else
            layer(llay)%ddepth=lsc%depth(llay)-(dg-1.D0)
            exit
         endif
      enddo
      evap%actm=0
      respe=evap%potm
      llay=0
      do
         llay=llay+1
         if (layer(llay)%ddepth==0.D0.or.llay>=nlayer) exit
         if (layer(llay)%sav>0.5) then
            amaxev=(.125+(layer(llay)%sav-.5))*(layer(llay)%snew(1)-layer(llay)%smin)
         else
            amaxev=.5*layer(llay)%sav*layer(llay)%sav*(layer(llay)%snew(1)-layer(llay)%smin)
         endif
         amaxev=amaxev*layer(llay)%ddepth/lsc%depth(llay)
         if (respe<=amaxev) then
            evap%actm=evap%actm+respe
            layer(llay)%snew(1)=layer(llay)%snew(1)-respe
            respe=0.D0
            ideeper=0
         else
            evap%actm=evap%actm+amaxev
            respe=respe-evap%actm
            layer(llay)%snew(1)=layer(llay)%snew(1)-amaxev
            ideeper=1
         endif
         if (ideeper/=1) exit
      enddo
		
		end