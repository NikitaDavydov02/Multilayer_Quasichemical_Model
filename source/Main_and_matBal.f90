implicit none			!   MAIN    PROGRAM
character*7 name_seg(150)
REAL*8 T,Pi,psi,k,e,Nav
COMMON /A/ T,Pi,psi,k,e,Nav
real*8 tot_energy,Xinput(150),r_max
integer s,I,M,NK,NG,NKK,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150),Nx,Nalf,Lmax,waterLmax,seg_ch
REAL*8 ZZ,ZR,JX,IJ,TL(150)
REAL*8 teta_chem(150,150),DPx(150),DPalf(150),EPS,EPF,Xbulk(150)
common /composition/Xbulk
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds  
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
common /forNEWT/Nx,Nalf,DPx,DPalf,EPS,EPF
REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid

T=298.15d0
Pi = 3.1415926d0
k = 1.380649d-16
e = 4.80286d-10
Nav = 6.022142d+23

Lmax=9                                !number of layers
r_solid=10.0d0																	
call indsurf(r_max) 		                   !parameters of system
TL=1.0d0 !r_max**0.333333d0 !1.0d0						!thickness of layers, may be different for every layer
TL(1)=1.0d0
!print *,r_max,TL(1),r_max**0.333333d0

call check_TL(TL,1.0d0,Lmax) !call check_TL(TL,r_max**0.333333d0,Lmax)			  !check minimal thickness, generate coordinates
! print *,r_max,TL(1)
!input data for phase equilubruim
NKK=Nmol						       !NKK=N(chains)+N(monomeric molecules)
	  DO I=1,NK
	  if (type_mol(I).eq.0) NKK=NKK+1 
	  enddo
if (NKK.eq.2) then
Xinput(1)=1.0d0
else
Xinput(2)=0.501d0 !1.0d0-Xinput(1) !-Xinput(2)			  !  initialize mixture composition
Xinput(3)=1.0d-19 !1.0d0-Xinput(1)-Xinput(2)
Xinput(1)=1.0d0-Xinput(3)-Xinput(2)
endif
Xinput(NKK)=0.0d0

!input parameters for NEWT
EPS=1.0D-10
EPF=40000.0d0
do I=1,NK
DPx(I)=1.0D-12            	 !!!!!-12  !step=DPx/Xloc						
enddo
DPx(NK)=1.0D-15
do s=1,NG-1
DPalf(NK-1+s)=1.0D-8					   !step=DPalf/Alfa_up
DPalf(NK+NG-2+Nbonds+s)=1.0D-8			   !step=DPalf/Alfa_dwn
enddo 
DO M=1,Nbonds
 DPalf(NK-1+NG-1+M)=1.0D-12  !1.0D-8				   !step=DPalf/chem_up
enddo

waterLmax=-1								 
!waterLmax=1 aqueous phase in the bulk (Lmax), =0 oil phase	in the bulk (Lmax), =-1 without calc LVE and Xbulk is composition for bulk

Xbulk(1)= 0.002698878 !0.30d0						  ! set arbitrary bulk composition
Xbulk(3)=0.9245223896 !0.50d0
Xbulk(2)=1.0d0-Xbulk(1)-Xbulk(3)

!Xbulk(1)=1.0d0-0.019d0 !1.0d0-0.019d0 !0.999d0
!Xbulk(2)=1.0d0-Xbulk(1) !-0.019d0
!Xbulk(3)=1.0d0-Xbulk(1)-Xbulk(2)
call micelle_energy('sph',TL,r_max,Lmax,Xinput,waterLmax,tot_energy)	
!(geometry(pln,sph or cyl),thickness of layer (in uinits of monomers),composition for LVE,final energy)		!TL=3A
end 
															 
!*************************************************************************************
SUBROUTINE init_profile(X2,X,Lmax)
      IMPLICIT NONE
	  character*7 name_seg(150)
	  real*8 pp,ZZ,ZR,JX,IJ,sigma,mediana(150),gauss,X2(150),X(150),teta_chem(150,150),PAR(150)
	  REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150)
	  integer L,Lmax,k,M,s,NK,NG,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150),phase,id(100),mainK,border
	  REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
	  real*8 Xold(150,150),AL_upold(10,150),AL_dwnold(10,150),AL_alold(10,150),ALlocold(10,150),chem_upold(150,150)
	  real*8 chem_up(150,150),chem_dwn(150,150),chem_al(150,150),ALchemold(150,150),Xchem(150),sumX(150)
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
	  COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	  common /chem_or/chem_up,ALchemold
	  common /chemdwn/chem_dwn,chem_al
	  common /old/Xold,AL_upold,AL_dwnold,AL_alold,ALlocold
      common /oldch/chem_upold
	  COMMON /C/T,TT,R,RS,Q,QS,vol_ref, q_ref
	  REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref

 border= Lmax/2
do L=1,Lmax
Xloc(nk,L)=0.0d0

do k=1,NK-1
!if (L.lt.Lmax-5) then
!Xloc(k,L)=(X(k)-X2(k))*(L-1)/(Lmax-1)+X2(k) !X2(k) +(X(k)-X2(k))*(Lmax-L)/Lmax
!else
!Xloc(k,L)=X(k)
!endif
 Xloc(k,L)=X2(k)

if (L.gt.border) Xloc(k,L)=X(k)
!if ((L.ge.16).and.(L.le.18)) Xloc(k,L)=X2(k)+0.1d0

!tanh profile
pp=L-(Lmax+1)/2 
!if (L.ne.Lmax) Xloc(k,L)= X2(k) +(X(k)-X2(k))*0.5d0*(1.0d0+dtanh(pp))

!gauss peaks
sigma=7.0d0																!thickness of peak (number of layers)
mediana(k)=border !15.0d0 !+(11.0d0- rank_seg(k))*0.005d0 !14.5d0 !rank_seg(k)+8.0d0
pp=L
!if (type_mol(k).eq.2) Xloc(k,L)=gauss(pp,mediana(k),1.06d-3,sigma/6.0d0,X2(k),X(k))  !(current layer, mediana,height,sigma,Xhom2,Xhom)
!if ((L.gt.14).and.(k.le.3)) Xloc(k,L)=(1.0d0 -X(NK)-10.0d0*Xloc(4,L))/3.0d0

Xloc(nk,L)=Xloc(nk,L)+Xloc(k,L)
enddo	 !k
Xloc(nk,L)=1.0d0-Xloc(nk,L)


Xloc(nk,L)=0.0d0
do k=1,NK-1
 !!!if ((L.gt.14).and.(k.le.3)) Xloc(k,L)=(1.0d0 -X(NK)-10.0d0*Xloc(4,L))/3.0d0

 Xloc(nk,L)=Xloc(nk,L)+Xloc(k,L)
enddo
 Xloc(nk,L)=1.0d0-Xloc(nk,L)
 enddo	!L

!read initial Xloc from inputX
!call readXinput(Lmax,NK,NG,Xloc)
 do L=1,Lmax
Xloc(nk,L)=0.0d0
do k=1,NK-1
Xloc(nk,L)=Xloc(nk,L)+Xloc(k,L)
enddo
Xloc(nk,L)=1.0d0-Xloc(nk,L)
 enddo

do k=1,NK
Xloc(k,Lmax+1)=Xloc(k,Lmax)
enddo



!guess initial orientations
call init_orient(Lmax,1)   !initial orientations without checking

do L=1,Lmax
do k=1,NK
PAR(k)=Xloc(k,L)
enddo
call check_borderX(PAR,L,Lmax,NK-1,id,mainK,1)

write (110,'(2X,I5,16(2X,E18.12))') L, Xloc(nk,L),(Xloc(k,L),k=1,NK-1)
  write (112,'(2X,I5,16(2X,E18.12))') L, Xloc(nk,L),(Xloc(k,L),k=1,NK-1)
enddo

call init_orient(Lmax,0)   !initial orientations with checking

!read initial orient from inputUP,DWN,AL
!call readALFinput(Lmax,NG,Nbonds,AL_up,AL_dwn,AL_al,chem_up,chem_dwn,chem_al)
 do k=1,NK
Xloc(k,Lmax+1)=Xloc(k,Lmax)
enddo

do L=1,Lmax
 do k=1,NK
Xold(k,L)=Xloc(k,L)
enddo
 
   write (112,*)'chem up,dwn,al'
   DO M=1,Nbonds
   chem_upold(M,L)=chem_up(M,L)
   ALchemold(M,L)=teta_chem(M,L)/ZZ
   write (112,'(I2,I2,16(4X,E18.12))')L,M,chem_up(M,L),chem_dwn(M,L),chem_al(M,L) !AL_up(4,L),AL_up(5,L)	
   enddo
   write (112,*)'AL up,dwn,al'	
   do s=1,NG
   AL_upold(s,L)=AL_up(s,L)
   AL_dwnold(s,L)=AL_dwn(s,L)
   ALlocold(s,L)=ALloc(s,L)
   write (112,'(I2,I2,16(4X,E18.12))')L,s,AL_up(s,L),AL_dwn(s,L),AL_al(s,L)			  
   enddo
 enddo !L

 return
 end
!**************************************************************************************************
 SUBROUTINE init_orient(Lmax,switch)
      IMPLICIT NONE
	  character*7 name_seg(150)
	  real*8 pp,ZZ,ZR,JX,IJ,teta_chem(150,150),XlocL(150),sumX(150),sumup,sumal,sumdwn,chem_min(150),chem_max(150)
	  REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
	  integer L,Lmax,k,s,NK,NG,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150),phase,seg_ch,KEY
	  integer Nx,Nalf,Nfix,CHANGE(150),No_or(150),kkkkk,switch
	  REAL*8 DPx(150),DPalf(150),EPS,EPF,PARalf(150)
	  REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
	  REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
	  real*8 chem_up(150,150),chem_dwn(150,150),chem_al(150,150),ALchemold(150,150),m_up,m_dwn,lat
	  COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
	  COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	   common /chem_or/chem_up,ALchemold
	  common /chemdwn/chem_dwn,chem_al
	   common /for_check_depend/KEY
	  common /ph/phase
	  common /forNEWT/Nx,Nalf,DPx,DPalf,EPS,EPF
	  common /without_or/ Nfix,No_or,CHANGE
	  	  REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
	  COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
	  kkkkk=0
	  pp=0.0d0
 Do L=1,Lmax+1
!35 write (110,'(2X,I5,16(2X,E18.12))') L, Xloc(nk,L),(Xloc(k,L),k=1,NK-1)
 ! write (112,'(2X,I5,16(2X,E18.12))') L, Xloc(nk,L),(Xloc(k,L),k=1,NK-1)
do k=1,NK
XlocL(k)=Xloc(k,L)
enddo
call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output

  DO s=1,NG
  ALloc(s,L)=0.0d0
  DO k=1,NK
  ALloc(s,L)=ALloc(s,L)+QS(k,s)*Xloc(k,L)/QT
  enddo 
 enddo
 !if (L.eq.1) print *,ALloc(1,L),ALloc(2,L),ALloc(3,L),QT 
  do s=1,NG
AL_up(s,L)=ALloc(s,L)
AL_dwn(s,L)=ALloc(s,L)
AL_al(s,L)=ALloc(s,L)			  !saving values in common to calc.XX on 1st iteration, then they will be updated		  
enddo
	
call local_chem(Xloc,L,QT) !,Nbonds,L,QT)		!writes teta_chem in common
enddo !L
!XlocL(1)=AL_up(1,1)+AL_up(2,1)+teta_chem(1,1)/ZZ+teta_chem(1,1)/ZZ+teta_chem(2,1)/ZZ+teta_chem(2,1)/ZZ
!XlocL(1)=1.0d0-XlocL(1)

!print *,XlocL(1)
Do L=1,Lmax
do k=1,NK
XlocL(k)=Xloc(k,L)
enddo
call calcALchem(XlocL,L,Xchem,sumX)		 !writes teta_chem in common
DO k=1,Nbonds
!if (teta_chem(k,L).gt.teta_chem(k,L+1)).and. then

if (teta_chem(k,L+1).ge.teta_chem(k,L)) then
 chem_up(k,L)=teta_chem(k,L)/ZZ
chem_dwn(k,L+1)=sumX(k)+chem_up(k,L)
else
!if (chem_dwn(k,L+1).le.0.0d0) GOTO 12
!!elseif ((teta_chem(k,L+1)/ZZ-sumX(k)).gt.0.0d0)	then
chem_dwn(k,L+1)=teta_chem(k,L+1)/ZZ
chem_up(k,L)=-sumX(k)+chem_dwn(k,L+1)
!if (chem_up(k,L).le.0.0d0) GOTO 11
!!else
!!chem_up(k,L)=sumX(k)*1.5d0
!!chem_dwn(k,L+1)=sumX(k)+chem_up(k,L)
endif
enddo
!XlocL(1)=AL_up(1,1)+AL_up(2,1)+teta_chem(1,1)/ZZ+teta_chem(1,1)/ZZ+teta_chem(2,1)/ZZ+teta_chem(2,1)/ZZ-1.0d0 !+chem_up(1,1)+chem_up(1,1)+chem_up(2,1)+chem_up(2,1)-1.0d0
!XlocL(1)=AL_up(1,1)+AL_up(2,1)+chem_up(1,1)+chem_up(2,1)+chem_dwn(1,2)+chem_dwn(2,2)-1.0d0
!print *,AL_up(1,1)+AL_up(2,1)+chem_up(1,L)+chem_up(2,L)chem_dwn(1,L+1)+chem_dwn(2,L+1)
!print *,1.0d0-AL_up(1,L)-AL_up(2,L)-chem_up(1,L)*2.0d0 -chem_up(2,L)*2.0d0 !-chem_up(1,L)-chem_up(2,L)-chem_dwn(1,L+1)-chem_dwn(2,L+1)
enddo !L


Do L=1,Lmax
do k=1,NK
XlocL(k)=Xloc(k,L)
enddo
11 call calcALchem(XlocL,L,Xchem,sumX)		 !writes teta_chem in common
 phase=0
if (L.eq.1) phase=-1
if (L.eq.Lmax) phase=1

 if ((L.eq.1).and.(r_solid.eq.0.0d0)) then	   
		DO k=1,Nbonds
		 chem_dwn(k,L)=0.0d0 
		  enddo
		  do s=1,NG
		AL_dwn(s,L)=0.0d0		  
		enddo    	
  endif	
  !!!!!!!!!!!!!!!!!!!!
 if (switch.lt.1) then
 DO k=1,Nbonds
 chem_min(k)=max(0.0d0,-sumX(k),(teta_chem(k,L+1)/ZZ/m_dwn(L+1)-sumX(k)-m_up(L+1)/m_dwn(L+1)-lat(L+1)/m_dwn(L+1)*0.5d0),(Xchem(k)/m_up(L)-m_dwn(L)*chem_dwn(k,L)/m_up(L)-lat(L)/m_up(L)*0.5d0))
chem_max(k)=min((Xchem(k)/m_up(L)-m_dwn(L)*chem_dwn(k,L)/m_up(L)),1.0d0,(1.0d0-sumX(k)),(teta_chem(k,L+1)/ZZ/m_dwn(L+1)-sumX(k)-m_up(L+1)/m_dwn(L+1)*chem_up(k,L+1)))	  
	if (chem_min(k).ge.chem_max(k))	then
	print *,'for initial chem_up min>max',k,l,chem_min(k), chem_max(k)
	print *,(teta_chem(k,L+1)/ZZ/m_dwn(L+1)-m_up(L+1)*chem_up(k,L+1)/m_dwn(L+1)-sumX(k))
	!Xloc(seg_ch(k),L)
	!GOTO 100
	stop
	endif
 !chem_min(D)=max(-sumX(D),(teta_chem(D,L+1)/ZZ/m_dwn(L+1)-m_up(L+1)/m_dwn(L+1)*chem_up(D,L+1)-sumX(D)-lat(L+1)/m_dwn(L+1)*0.5d0),(Allocsol(D)/m_up(L)-m_dwn(L)*chem_dwn(D,L)/m_up(L)-lat(L)/m_up(L)*0.5d0))
 !  chem_max(D)=min((1.0d0-sumX(D)),(Allocsol(D)/m_up(L)-m_dwn(L)*chem_dwn(D,L)/m_up(L)),(teta_chem(D,L+1)/ZZ/m_dwn(L+1)-sumX(D)-m_up(L+1)/m_dwn(L+1)*chem_up(D,L+1) ))	  	 

do while (chem_up(k,L).le.chem_min(k)) 
	chem_up(k,L)=dabs(chem_up(k,L))*1.001d0
		do while  (chem_up(k,L).ge.chem_max(k)) 
		chem_up(k,L)=chem_up(k,L)*0.999d0				 !
		enddo
	enddo
	do while (chem_up(k,L).ge.chem_max(k)) 
	chem_up(k,L)=chem_up(k,L)*0.999d0
		do while  (chem_up(k,L).le.chem_min(k)) 
		chem_up(k,L)=dabs(chem_up(k,L))*1.001d0				!
		enddo
	enddo
 chem_dwn(k,L+1)=sumX(k)+chem_up(k,L)

chem_al(k,L)=(teta_chem(k,L)/ZZ-m_up(L)*chem_up(k,L)-m_dwn(L)*chem_dwn(k,L))/lat(L) 
!endif

enddo	  !k
 	
call correct_ALF_forNG(L,AL_up,AL_dwn,ALloc,teta_chem,chem_up,chem_dwn,sumup,sumdwn,sumal,0,PARalf) !changes Al_up,AL_dwn

!endif

 call calcAL_al_and_NG(L)
 endif
  enddo	 !L
	
return
 end
 !**************************************************************************************************
SUBROUTINE readXinput(Lmax,NK,NG,Xloc)
IMPLICIT NONE
	  real*8  Xloc(150,150),PAR(150),Xchem(150),sumX(150)
	  integer Lmax,NK,jj,lay,k,s,NG
	  COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
	  REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
	  REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
	  COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	  
 OPEN (UNIT=1,FILE='inputX.txt')
 !read (1,'(A2,(4X,A14,4X),30(10X,I2,A5,5X))')
 do	lay=1,Lmax-1
 read (1,'(71X,50(4X,E18.12))')Xloc(nk,Lay),(Xloc(jj,Lay),jj=1,NK-1)

 Xloc(nk,Lay)=0.0d0
	do k=1,NK-1
	Xloc(nk,Lay)=Xloc(nk,Lay)+Xloc(k,Lay)
	enddo	!k
 Xloc(nk,Lay)=1.0d0-Xloc(nk,Lay)
	 do k=1,NK
	PAR(k)=Xloc(k,Lay)
	enddo

call total(PAR,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output

  DO s=1,NG
  ALloc(s,Lay)=0.0d0
  DO k=1,NK
  ALloc(s,Lay)=ALloc(s,Lay)+QS(k,s)*Xloc(k,Lay)/QT
  enddo 
  enddo
  call calcALchem(PAR,Lay,Xchem,sumX)		 !writes teta_chem in common

 enddo	   !lay
 return
 end
!*************************************************************************************
SUBROUTINE readALFinput(Lmax,NG,Nbonds,AL_up,AL_dwn,AL_al,chem_up,chem_dwn,chem_al)
IMPLICIT NONE
	 REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150)
	 real*8 chem_up(150,150),chem_dwn(150,150),chem_al(150,150)
	  integer Lmax,NG,Nbonds,s,k,lay

OPEN (UNIT=2,FILE='inputUP.txt')
OPEN (UNIT=3,FILE='inputDWN.txt')
!OPEN (UNIT=4,FILE='inputAL.txt')
!read (2,'(A2,30(5X,I2,A9,6X))')
!read (3,'(A2,30(5X,I2,A9,6X))')
!read (4,'(A2,30(5X,I2,A9,6X))')

do	lay=1,Lmax-1
read (2,'(2X,57(4X,E18.12))')(AL_up(s,Lay),s=1,NG),(chem_up(k,Lay),k=1,Nbonds) 
read (3,'(2X,57(4X,E18.12))')(AL_dwn(s,Lay),s=1,NG),(chem_dwn(k,Lay),k=1,Nbonds) 
!read (4,'(2X,17(4X,E18.12))')(AL_al(s,Lay),s=1,NG),(chem_al(k,Lay),k=1,Nbonds)
enddo
do	lay=1,Lmax
call calcAL_al_and_NG(Lay) 
enddo
 return
 end
!*************************************************************************************
SUBROUTINE	 check_TL(TL,tau_min,Lmax)
IMPLICIT NONE
	 REAL*8 TL(150),tau_min,Tcurrent
	 integer L,Lmax,iter,final,final1
	  common /cycles/iter,final,final1
	 OPEN (UNIT=22,FILE='coordinates of layers.txt')
	 write(22,'(2X,A5,9X,A5,14X,A12)')'layer','coord','tau of layer'
	 Tcurrent=0.0d0
	do L=1,Lmax
	if (TL(L).lt.tau_min) then
	print *,'thickness of layer ',L,'<min',TL(L),tau_min
!	TL(L)=tau_min
	endif
	Tcurrent=Tcurrent+TL(L)
	if (final /= 0) write(22,'(2X,I3,4X,E18.12,4X,E18.12)')L,Tcurrent,TL(L)
	write(22,'(2X,I3,4X,E18.12,4X,E18.12)')L,Tcurrent,TL(L)					  
	enddo
if (final /= 0)	write(22,'(A8,E18.12)')'tau min=',tau_min
if (final /= 0)	close (22)
	 return
 end
 
!*************************************************************************************
