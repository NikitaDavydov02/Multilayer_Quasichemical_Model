! ***************************  Numerical **************************************************
  REAL*8 FUNCTION FU_tau (X)
  IMPLICIT NONE
      REAL*8   PARXlocREM(150,150),PARTLmax,X(1),TL(150),TLmin,PAR(150),FU(150),ghom(150),FUtau,PARalf(150),DPalf(150),DDalf(150,150),EPS,EPF,PARinit(150),pN(150),Pnormal,vol1,Phi,vol2,volume ,dpT(150,150),PARTL,PARXloc(150,150),vol1max,vol2max
	  integer Lmax,L,Nx,Nalf,switch,k,Nbonds,KEY,iter,final,final1,id(100),mainK
	  common /layer/shape,TL,Lmax,TLmin
	   CHARACTER*3 shape
	   common /layer_f/Nx,Nalf,Nbonds,PARinit,L,ghom,DPalf,EPS,EPF,vol1,Phi,PARTL,PARTLmax,PARXloc
	   common /cycles/iter,final,final1
	   common /pressure/pN,dpT
	   real*8 ALlocnew(150),sumXX(150)
	 !  common/mechanical/vol1,Phi
	  ! print *,Nx,Nalf,PARinit(1:Nx),L,X(1),PARTL+X(1)*1.0d-15
	!    print *,'dzreal',X,FU_tau
	 PAR=PARinit
	 PARXlocREM=PARXloc
	!TL(L)=X(1)
	TL(L)=PARTL
	TL(Lmax)=PARTLmax
	vol1max=volume(Lmax)
	TL(L)=PARTL+X(1)*1.0d-5	
!	print *,vol1,Phi
	 ! 	Phi=JX
	!vol1=volume(L)
	!DPx(NK)=DDPx(NK)                                           
	!  TL(L)=X(NK)+DPx(NK)
	if (dabs(TL(L)).gt.1.0d10) then
	FU_tau=2.0d5
	else
	vol2=volume(L)
	do k=1,Nx
	PAR(k)=PARinit(k)/(1.0d0+(vol2-vol1)/vol1/Phi) !*(1.0d0-(vol2-vol1)/vol1/Phi)
	enddo
	 PAR(Nx)=PAR(Nx)+(vol2-vol1)/vol1/Phi/(1.0d0+(vol2-vol1)/vol1/Phi) 
!	print *,'!!!!',PAR(1:Nx),L,Lmax
	     call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR

	 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
	! print *,'!!!!',mainK,switch
	  if ((mainK.eq.0).and.(switch.eq.0))	   then
	!	TL(Lmax)=PARTLmax-X(1)*1.0d-5
	!!!	vol2max=volume(Lmax)
	!		print *,(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi
	!!!	do k=1,Nx
	 !!!   PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2max-vol1max)/vol1max/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	  !!!  enddo
	  !!!   call check_borderX(PAR,Lmax,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	!!!	 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,Lmax,switch)
		  call calcALloc(PAR,L,ALlocnew,sumXX)
	 
	  DO  k=Nx+1,Nx+Nalf+Nbonds				  !Nx=NK-1, Nalf=2*(NG-1)
	  PAR(k)=PARXlocREM(k,L)					!creating unified array PARth
	  enddo 
	!   if (L.eq.17) stop
	   print *,'dzreal'

 call check_borderALF(PAR,L,Nalf,ALlocnew,sumXX,id,mainK,0)	!0 = corrections   PARth
 print *,'dzreal'

	!	 print *,'!!!!',mainK,switch,vol2max-vol1max, PAR(1:Nx),PARXloc(1:Nx,Lmax)
		endif
	!  print *,TL(L),switch
	 if ((switch.eq.1).or.(mainK.eq.1)) then
	 FU_tau=2.0d5
	 !if (X(1).ge.1.0d0) FU_tau=2.0d5*dabs(X(1))
	! if (X(1).lt.1.0d0) FU_tau=2.0d5/dabs(X(1))
	 else
	 !if (final1.eq.0) then
		 do k=1,Nx
	PAR(k)=PARinit(k)/(1.0d0+(vol2-vol1)/vol1/Phi) !!!*(1.0d0-(vol2(L)-vol1(L))/vol1(L)/Phi)
	enddo
	PAR(Nx)=PAR(Nx)+(vol2-vol1)/vol1/Phi/(1.0d0+(vol2-vol1)/vol1/Phi)
!	print *,'dzreal'
	   call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	   print *,'dzreal'
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo
	  PARXloc(1:(Nx+Nalf+Nbonds),L)=PAR(1:(Nx+Nalf+Nbonds))
	!		print *,(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi
	!!!  	do k=1,Nx
	 !!   PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2max-vol1max)/vol1max/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	 !!!   enddo
		!  call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,Lmax,EPS,EPF,KEY,ghom)
	!	 DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      !PAR(K)=PARalf(K)
	 ! enddo  
	 ! PARXloc(1:(Nx+Nalf+Nbonds),Lmax)=PAR(1:(Nx+Nalf+Nbonds))   
	pN(L)=Pnormal(L,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
	pN(Lmax)=Pnormal(Lmax,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom) 
	PAR(1:(Nx+Nalf+Nbonds))=PARXloc(1:(Nx+Nalf+Nbonds),L) 
		CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)
      FU_tau=FUtau !dabs(FUtau) !**2
	  endif
	  endif
	  PARXloc=PARXlocREM
	 print *,'dzreal',X,FU_tau
	 TL(L)=PARTL
	 PAR=PARinit
	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
      RETURN
      END
 
!**************************
  REAL*8 FUNCTION FU_LVEp (X)
  IMPLICIT NONE
      REAL*8       X(1),Xinn(150),Xin(150),P,PP,sum,Xout2(150),Xout1(150),T
	  integer NK,I,NKK
	   common /equib/Xinn,NK,NKK,P,PP,T
	  Xin=Xinn
	  sum=0.0d0
	  if ((NKK-2).ge.2) then
	  do I=2,NKK-2
	  sum=sum+Xin(I)
	  enddo
	  endif
	  Xin(NKK-1)=X(1)
	 if ((X(1).le.1.0d-150).or.((sum+X(1)).ge.1.0d0)) then						 !.or.(X(1).gt.0.5d0)
	 FU_LVEp=1.0d20 !*X
	 else
	if ((NKK-1).ne.1) Xin(1)=1.0d0-sum-X(1)
	  call LVE(Xin,NK,Xout2,Xout1)

      FU_LVEp=P-PP
	  endif
	 print *,X,FU_LVEp,PP
      RETURN
      END
!**************************
 SUBROUTINE NEWTP(Nx,Nalf,PARx,PARalf,DDPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)	   !searching Xloc in layer L
	 IMPLICIT NONE
	  CHARACTER*3 shape
	  character*7 name_seg(150)
	   integer Lay,Nxx,Nalff,Nx,Nalf,KEY,L,I,J,KK,kkk,K,HH,chem(150),type_mol(150),phase,rank_seg(150),Lmax,D,iter,final,final1
	  REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150)
	  REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
      REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150),Falfeq,PARalfeq(150),Fxeq,PARxeq(150)
	  REAL*8 TL(150),EPS,EPF,A,DDalf(150,150),TLmin                                           
	  REAL*8 PPAR(150),PAR(150),PARR(150),FU(150),DD(150,150),X(150),F(150),FF(150),LAM,PARx(150),PARalf(150),DPx(150),DDPx(150),DPalf(150),ghom(150),PAR1(150),gghom(150)
	  integer NK,NG,Nmol,Nbonds,allK,KEYx(100),KEYxend(100),KEYxL(100),KEYid(100),id(100),mainK,KEY_step,switch 
	  REAL*8 ZZ,ZR,JX,IJ,FUtau,minn,m_up,m_dwn,lat,XMIN
	  common /layer/shape,TL,Lmax,TLmin
      COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
      COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  common /layer_one/PAR1
	  common /ph/phase
	   common /chem_or/chem_up,ALchemold
	   common /chemdwn/chem_dwn,chem_al	
	   common /cycles/iter,final,final1
	    common /eq/Falfeq,PARalfeq,Fxeq,PARxeq
	!	common /layer_f/Nxx,Nalff,PPAR,Lay,gghom
		external FU_tau
		Nxx=Nx
		Nalff=Nalf
		gghom=ghom
		Lay=L
	  ! TLmin=0.7d0 !(m_UP(L)+m_dwn(L))*TL(L)
	!	 common /start_of_NEWTX/iter_again
	!	 if (iter.gt.2) EPF=5.0d-2                                           
      KEY=0                                                             
      KK=0                                                              
      LAM=1.0d0 
	  HH=0
	  minn=100.0d0 
!	 print *,'start of NEWTx'	                                                    !Nx=NK-1, Nalf=2*(NG-1)
	  DO  I=1,Nx                                                       
      PAR(I)=PARx(I)
	  if (PARx(I).lt.minn) minn=PARx(I)
	  enddo
!	    if ((L.eq.17).and.(iter.gt.0)) then
!	print *,'NEWTx'
!	endif
	   call check_borderX(PAR,L,Lmax,Nx,id,mainK,0)	   !0= correction of initial PAR	
	

 	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	   

!	 write(*,*)'!!!!!!!!!!!!!! first NEWTalf ended'  
!	 stop 
	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(I)=PARalf(I)
	  enddo
	    PAR(Nx+Nalf+Nbonds+1)=TL(L)
	   do while ((PAR(NK+Nalf+Nbonds).lt.TLmin).or.(lat(L).le.0.0d0)	  )	
	   TL(L)=TL(L)*1.1d0
	   enddo
	  PARx(NK)=TL(L)
	   PAR(Nx+Nalf+Nbonds+1)=TL(L)
	!    print *,Xmin,L

	   CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)
	    !DO D=1,Nx
	  ! if (dabs(FU(D)).gt.1.0d-12) HH=1
	   if (dabs(FUtau).gt.1.0d-5) HH=1
	   !enddo
! iter_again=1 
!	    if ((L.le.18).and.(L.ge.17).and.(iter.lt.200)) HH=0	  !((L.le.1).and.(iter.lt.100)) 
	   if (L.eq.Lmax) HH=0
!	   if (L.eq.1) HH=0
!	  if (L.le.8) HH=0
PPAR=PAR
	! if (dabs(FUtau).gt.1.0d-5) then

	! call duvmgs(FU_tau, TLmin, 2000.0d0*PARx(NK), 1.0d-10, XMIN)
	! TL(L)=Xmin
	! print *,Xmin,L
	! endif
     if (HH /= 0) then                                                                                                                                                   
   5  DO  I=1,Nx                                                        
      F(I)=FU(I)                                            
      X(I)=PAR(I)
	  enddo
	  DO I=Nx+1,Nx+Nalf+Nbonds+1
	  X(I)=PAR(I)
	  enddo 
	  F(NK)=FUtau                                                     
      KK=KK+1
!	  print *,PAR(5:6),AL_up(2,L),AL_al(2,L),chem_al(1,L),chem_al(2,L)                                                          
    !  PAR=X
 
 !tau
  DPx(NK)=DDPx(NK) !*PAR(I)                                                    
      PAR(NK+Nalf+Nbonds)=X(NK+Nalf+Nbonds)+DPx(NK)
	   TL(L)=PAR(NK+Nalf+Nbonds)
	   	  call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
	  do while (switch.ne.0)
	  !do while (PAR(NK+Nalf+Nbonds).lt.max(TLmin,(m_up(L)+m_dwn(L))*TL(L)))
	  !do while ((PAR(NK+Nalf+Nbonds).lt.TLmin).or.((PAR(NK+Nalf+Nbonds)-PARx(NK)).gt.PARx(NK)).or.(lat(L).le.0.0d0).or.(switch.ne.0)	  )					!.or.((PAR(NK+Nalf+Nbonds)-PARx(NK)).gt.PARx(NK))
	  if (dabs(DPx(NK)).lt.1.0d-30) then				 !-16
		DPx(NK)=-DDPx(NK)
		endif
	  if ((PAR(NK+Nalf+Nbonds).lt.TLmin).or.(lat(L).le.0.0d0)	  ) then
	  	if (DPx(NK).gt.0.0d0) DPx(NK)=DPx(NK)*1.1d0
		if (DPx(NK).lt.0.0d0) DPx(NK)=DPx(NK)*0.9d0
	!	elseif   (lat(L).le.0.0d0) then
		!DPx(NK)=DPx(NK)*0.9d0
		!print *,DPX(NK)
		else
		if (DPx(NK).gt.0.0d0) DPx(NK)=DPx(NK)*0.9d0
		if (DPx(NK).lt.0.0d0) DPx(NK)=DPx(NK)*1.1d0
	 endif
	   PAR(NK+Nalf+Nbonds)=X(NK+Nalf+Nbonds)+DPx(NK)
	    TL(L)=PAR(NK+Nalf+Nbonds)
	  call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)

	  enddo
	  TL(L)=PAR(NK+Nalf+Nbonds)
	   call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)

	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	!  print *,I,DPx(I),PARalf(K)
	  enddo

 CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)
                                            
  !    DO  J=1,Nx                   
!	  DD(J,NK)=(FU(J)-F(J))/DPx(NK) !*X(I) !(dlog(PAR(I)/X(I)) )  !DPx(I) 		! DD(J,I) changed!!!
!	  print *,J,DD(J,NK),DPX(NK)
!	  enddo                  !cycle J 
	  DD(1,1)=(FUtau-F(NK))/DPx(NK)
	  print *,NK,DD(1,1) !,(Al_up(1,L)+Al_dwn(1,L)-2.0d0*Alloc(1,L)),(Al_up(2,L)+Al_dwn(2,L)-2.0d0*Alloc(2,L)),(Al_up(3,L)+Al_dwn(3,L)-2.0d0*Alloc(3,L)) !DPX(NK)                                                           
      PAR(NK+Nalf+Nbonds)=PAR(NK+Nalf+Nbonds)-DPx(NK) !*X(I)
 !tau
	  do I=1,NK
	  PARR(I)=0.0d0 !PAR(I)
	  enddo
	  FF(1)=-F(NK)
!	   DO  J=1,NK                   
!	   FF(J)=-F(J)
!	   enddo                  !cycle J 
  	 !    if ((L.eq.2).and.(iter.eq.13)) then 
!	 write(*,*)'!!!!!!!!!!!!!! NEWTalf after dx ended',DPx(1:3) !,DD(1:3,1:3)
!	 write(*,*) ((DD(i,j), J=1,3),I=1,3) 
	 !endif  																					 
!	   print *,((DD(i,j),j=1,Nx),i=1,Nx)
!		  stop
	! if (DD(NK,NK).eq.0.0d0) then
	!  CALL DLSLRG(Nx, DD, 150, FF, 1, PARR)
	!  else
	if (DD(1,1).ne.0.0d0) then   
	  CALL DLSLRG(1, DD, 150, FF, 1, PARR)
	  	  PAR(NK+Nalf+Nbonds)=X(NK+Nalf+Nbonds)+PARR(1)
		  TL(L)=PAR(NK+Nalf+Nbonds)
		  	   	  call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
		  print *,parr(1)
		  do while (switch.ne.0)
	! do while (PAR(NK+Nalf+Nbonds).lt.max(TLmin,(m_up(L)+m_dwn(L))*TL(L)))
	  ! do while ((PAR(NK+Nalf+Nbonds).lt.TLmin).or.((PAR(NK+Nalf+Nbonds)-PARx(NK)).gt.PARx(NK)).or.(lat(L).le.0.0d0).or.(switch.ne.0)	  )		  !.or.((PAR(NK+Nalf+Nbonds)-PARx(NK)).gt.PARx(NK))
	 ! if (PAR(NK+Nalf+Nbonds).lt.TLmin) then
	  	PARR(1)=PARR(1)*0.8d0
	  !endif
	   PAR(NK+Nalf+Nbonds)=X(NK+Nalf+Nbonds)+PARR(1)
	   TL(L)=PAR(NK+Nalf+Nbonds)
	  	   	  call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
	  enddo
	  endif
	 ! TL(L)=PAR(NK+Nalf+Nbonds)
	! print *,L,TL(L),PARR(1)    

    A=DABS(PARR(NK)) !0.0d0 
	
if (kkk.lt.150)   call ALFA_step (Nx,Nalf,X,PARalf,PARR,0,L)
if (kkk.ge.150)	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)	  
	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(I)=PARalf(I)
	  enddo         
	  ! write(*,*)'x=',kkk
	  CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom) 
 
	  write(*,*)'!!!Layer = ',L
	                                            
     ! A=A/Nx  
	  IF(KK.GT.1000) write(*,*)'>1000 iterations in NEWTx',A
	 ! if (iter.ne.1) then 
	 ! do k=1,Nx
	!  if (PARR(k).gt.(X(k)*1.0d-5)) GOTO 5
	!  enddo
	!  else                                                          
      IF(A.GT.EPS)  GOTO 5
	!  endif                                             
      A=0.0d0
!	  if (Nfix) then
 !  do I=1,Nfix
 !   A=A+DABS(FU(Nx+No_or(I)))+DABS(FU(Nx+Nalf/2+Nbonds+No_or(I)))
 !  enddo 
!   endif                                                               
      DO 6 I=1,Nx 
!	  IF(KK.GT.180) write(*,*)'200',DABS(FU(I))                                                     
   6  A=A+DABS(FU(I))
 
   write(*,*)'A=sum of FU =',A                                          
      KEY=1                                                             
      IF(A.GT.EPF)  KEY=-1.0d0
	  
	  if (KEY.lt.0.0d0) GOTO 5
	!  	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  
!	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
   !   PAR(I)=PARalf(I)
!	  enddo         
!	  CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom) 

  endif          !HH  	
99	  DO  I=1,Nx                                                      
      PARx(I)=PAR(I)
	  enddo                                       
   7  CONTINUE                    
                           
      RETURN                                                            
      END                                         
!********************************************************************************************************************
REAL*8 function Pnormal(L,Xloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
implicit none
real*8 vol1(150),vol2(150),volume,TL(150),DP,PARx,TLmin,EPS,EPF,ghom(150),DPalf(150),Xloc(150,150),Phi,PAR(150),PARalf(150),FU(150),DDalf(150,150),free_energy,ZZ,ZR,JX,IJ,ALlocnew(150),sumXX(150),dpT(150,150),pN(150) 
integer L,LL,Lmax,Nx,Nalf,k,kkk,KEYid(100),id(100),mainK,KEY,NK,NG,Nmol,Nbonds,switch,phase
real*8 loc_tens(150),tension,grand_potential,delta,loc_tensold
character*3 shape
common /layer/shape,TL,Lmax,TLmin
COMMON /D/ ZZ,ZR,NK,NG,Jx,IJ,Nmol,Nbonds
	common /pressure/pN,dpT
	 do k=1,Nx
	PAR(k)=Xloc(k,L)
	enddo
	call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      Xloc(K,L)=PARalf(K)
	  enddo
	Pnormal=0.0d0
	do LL=1,Lmax
Pnormal=Pnormal-free_energy(LL)
!print *,'free',L,free_energy(LL),LL !,JX
enddo
!print *,'free1',L,Pnormal,volume(L)
Phi=JX
	vol1(L)=volume(L)
	PARx=TL(L)
	DP=1.0d-7
	if (L.ne.Lmax) then
		do LL=L+1,Lmax
	!	CALL grandpot_pressure(LL,ghom,loc_tens,tension,grand_potential,delta)
	!	dPt(LL,L)=-loc_tens(LL)*volume(LL)
		vol1(LL)=volume(LL)
	!		if (LL.eq.(L+1))	print *,'!1',LL,L,dPt(LL,L),loc_tens(LL)
		enddo
	endif	
	  TL(L)=PARx+DP !*PARX

	!PAR=XlocL
1	vol2(L)=volume(L)
LL=L
	!	do LL=L,Lmax
		vol2(LL)=volume(LL)
	!	print *,vol2(LL),vol1(LL),LL
	do k=1,Nx
	PAR(k)=Xloc(k,LL)/(1.0d0+(vol2(LL)-vol1(LL))/vol1(LL)/Phi) !!!Xloc(k,LL)*(1.0d0-(vol2(LL)-vol1(LL))/vol1(LL)/Phi)
	enddo
	 PAR(Nx)=PAR(Nx)+(vol2(LL)-vol1(LL))/vol1(LL)/Phi/(1.0d0+(vol2(LL)-vol1(LL))/vol1(LL)/Phi) 

	     call check_borderX(PAR,LL,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
		  call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,LL,switch)
	 kkk=0
	 do while ((mainK.ne.0).or.(switch.ne.0))
	 print *,DP,PAR(1:Nx),LL,L
	 kkk=kkk+1
	 DP=DP*0.9d0
	 	if (dabs(DP).lt.1.0d-20) then	!-16
		DP=-DP*1.0d13 !DP*1.0d10
		endif 
	  TL(L)=PARx+DP
	  !	print *,LL,PAR(1:Nx)
	  if (LL.ne.L) goto 1
	!  print *,TL(L)
	vol2(LL)=volume(LL)
	do k=1,Nx
	PAR(k)=Xloc(k,LL)/(1.0d0+(vol2(LL)-vol1(LL))/vol1(LL)/Phi) !!!Xloc(k,LL)*(1.0d0-(vol2(LL)-vol1(LL))/vol1(LL)/Phi)
	enddo
	 PAR(Nx)=PAR(Nx)+(vol2(LL)-vol1(LL))/vol1(LL)/Phi/(1.0d0+(vol2(LL)-vol1(LL))/vol1(LL)/Phi) 
	  call check_borderX(PAR,LL,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	  if (LL.eq.L)  call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,LL,switch)
	  if (kkk.gt.1000) then
	  Pnormal=0.0d0
	  goto 11
	  endif
	 enddo		!while
	! print *,TL(L),DP
!	 call calcALloc(PAR,L,ALlocnew,sumXX)
!	call check_borderALF(PAR,L,Nalf,ALlocnew,sumXX,id,mainK,0)	!0 = corrections   PARth	 	
!	CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)
!	print *,PAR(1),Xloc(1,LL),LL
	  !!! call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,LL,EPS,EPF,KEY,ghom)
!	   enddo	!LL
	   do k=1,Nx
	PAR(k)=Xloc(k,L)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!Xloc(k,LL)*(1.0d0-(vol2(LL)-vol1(LL))/vol1(LL)/Phi)
	enddo
	PAR(Nx)=PAR(Nx)+(vol2(L)-vol1(L))/vol1(L)/Phi/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi)
	   call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo
!	  print *,'!!'
	vol2(L)=volume(L)
	do LL=1,Lmax
 phase=0
if (LL.eq.1) phase=-1
if (LL.eq.Lmax) phase=1
call IKS_or(phase,LL)		   !solve quasichemical equations,calc XX_or(L)
enddo
	do LL=1,Lmax
	Pnormal=Pnormal+free_energy(LL)
!	print *,'free',L,free_energy(LL),LL,vol2(L)-vol1(L)
	enddo
	!print *,'free',L,Pnormal !,free_energy(Lmax-1),free_energy(Lmax)
!		CALL grandpot_pressure(L,ghom,loc_tens,tension,grand_potential,delta)
!	Pnormal=Pnormal-loc_tens(L)*vol2	
!	print *,vol2-vol1,free_energy(L)
	Pnormal=-Pnormal/(vol2(L)-vol1(L))
!	print *,Pnormal
	   if (L.ne.Lmax) then
		do LL=L+1,Lmax
	!	CALL grandpot_pressure(LL,ghom,loc_tens,tension,grand_potential,delta)
	!	dPt(LL,L)=	dPt(LL,L)+loc_tens(LL)*volume(LL)
	!	dPt(LL,L)=	-dPt(LL,L)/(vol2(L)-vol1(L))
	!if (LL.eq.(L+1))	print *,'!',LL,L,dPt(LL,L),loc_tens(LL)
		enddo 
		endif
11	TL(L)=PARx
	!!!	do LL=L,Lmax
    do k=1,Nx
	PAR(k)=Xloc(k,L)
	enddo
	call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      Xloc(K,L)=PARalf(K)
	  enddo
	  	do LL=1,Lmax
 phase=0
if (LL.eq.1) phase=-1
if (LL.eq.Lmax) phase=1
call IKS_or(phase,LL)		   !solve quasichemical equations,calc XX_or(L)
enddo

	!!!  enddo
	   !PAR(1:(Nx+Nalf+Nbonds))
	end
!*******************************************************************************************************************
	  SUBROUTINE NEWT_mechanical(Nx,Nalf,PARx,PARalf,DDPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)	   !searching thickness of layer L
	  IMPLICIT NONE
	  CHARACTER*3 shape
	  character*7 name_seg(150)
	   integer Nx,Nalf,KEY,L,I,J,KK,kkk,K,HH,chem(150),type_mol(150),phase,rank_seg(150),Lmax,D,iter,final,final1,switch,info(1),Nxx,Nalff,NNbonds,Lay
	  REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150),free_energy,FUt,vol1(150),vol2(150),volume
	  REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
      REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150),Falfeq,PARalfeq(150),Fxeq,PARxeq(150),SUMXX(150),ALLOCNEW(150)
	  REAL*8 TL(150),EPS,EPF,A,DDalf(150,150),dpT(150,150)                                           
	  REAL*8 PAR(150),PARR(150),FU(150),DD(150,150),X(150),F(150),FF(150),LAM,PARx(150),PARalf(150),DPx(150),DDPx(150),DPalf(150),ghom(150),PAR1(150),loc_tens_old(150),PARXloc(150,150),Xlmax(150),TLmax
	  integer NK,NG,Nmol,Nbonds,allK,KEYx(100),KEYxend(100),KEYxL(100),KEYid(100),id(100),mainK,KEY_step 
	  REAL*8 ZZ,ZR,JX,IJ,TLmin,minn,FUtau,ALlocREM(10,150),teta_chemREM(150,150),QTrem,QT,RT,XlocL(150),Phi,pN(150),Pnormal,X1(1),Xmin(1),FU_tau,PPAR(150),gghom(150),DDPalf(150),EEPS,EEPF,TLx,voll1 ,PPARXloc(150,150)
	  common /layer/shape,TL,Lmax,TLmin
      COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
      COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  common /layer_one/PAR1
	  common /ph/phase
	   common /chem_or/chem_up,ALchemold
	   common /chemdwn/chem_dwn,chem_al	
	   common /cycles/iter,final,final1
	    common /eq/Falfeq,PARalfeq,Fxeq,PARxeq
	   common /for_steps/ALlocREM,teta_chemREM,QTrem
	common /pressure/pN,dpT
	external FU_tau
	common /layer_f/Nxx,Nalff,NNbonds,PPAR,Lay,gghom,DDPalf,EEPS,EEPF,voll1,Phi,TLx,TLmax,PPARXloc
		Nxx=Nx
		Nalff=Nalf
		gghom=ghom
		Lay=L
		DDPalf=DPalf
		EEps=EPS
		EEPF=EPF
		NNbonds=Nbonds
		PPAR=PARx
		
	  PARXloc=Xloc
	  PAR=PARx
	pN(L)=Pnormal(L,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
	pN(Lmax)=Pnormal(Lmax,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
	PARx(NK)=TL(L)
	TLx=TL(L)
	TLmax=TL(Lmax)
	PAR(1:(Nx+Nalf+Nbonds))=PARXloc(1:(Nx+Nalf+Nbonds),L)
	 PPARXloc=PARXloc
	CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)
!	print *,'Xloc=',Xloc(1,L),L 
	kk=0
	if (dabs(FUtau).gt.1.0d-12) HH=1
	if (pN(L).eq.0.0d0) HH=0
!	HH=0
		X1(1)=1.0d-7
	Phi=JX
	vol1(L)=volume(L)
	vol1(Lmax)=volume(Lmax)
	voll1=vol1(L)
!	print *,X1(1)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (HH.eq.1) then
	CALL DZREAL (FU_tau, 1.0d-5, 1.0d-2, 1.0d-2, 1.0d-2, 1, 1000,X1, Xmin, INFO)		   !min value FUtau,  relative change of two successive approximations
	 TL(L)=PARx(NK)+Xmin(1)*1.0d-5
	 print *,Xmin(1),L
!	 endif
	  if (TL(L).lt.TLmin) TL(L)=TLmin
	  vol2(L)=volume(L)
	do k=1,Nx
	PAR(k)=PARx(k)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!PARx(k)*(1.0d0-(vol2-vol1)/vol1/Phi)
	enddo
	 PAR(Nx)=PAR(Nx)+(vol2(L)-vol1(L))/vol1(L)/Phi/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) 

		    call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
		 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
	    call calcALloc(PAR,L,ALlocnew,sumXX)
	 
	  DO  k=Nx+1,Nx+Nalf+Nbonds				  !Nx=NK-1, Nalf=2*(NG-1)
	  PAR(k)=PARXloc(k,L)					!creating unified array PARth
	 enddo 
	!   if (L.eq.17) stop

 call check_borderALF(PAR,L,Nalf,ALlocnew,sumXX,id,mainK,0)	!0 = corrections   PARth
	  if (mainK.ne.0) then
	   TL(L)=PARx(NK)
	   PAR=PARx
	   endif
	  !!	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  
	!!  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
    !!  PAR(I)=PARalf(I)
	!!  enddo 
	!!  if ((mainK.eq.0).and.(switch.eq.0))	   then
!		TL(Lmax)=TLmax-Xmin(1)*1.0d-5
	!	vol2(Lmax)=volume(Lmax)
	!		print *,(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi
	!	do k=1,Nx
	!    PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	 !   enddo
	   !!  call check_borderX(PAR,Lmax,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	!!	 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,Lmax,switch)
	!!	endif
!	 enddo
	 do k=1,Nx
	PAR(k)=PARx(k)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!*(1.0d0-(vol2(L)-vol1(L))/vol1(L)/Phi)
	enddo
		 PAR(Nx)=PAR(Nx)+(vol2(L)-vol1(L))/vol1(L)/Phi/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) 
	   call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo
	  PARXloc(1:(Nx+Nalf+Nbonds),L)=PAR(1:(Nx+Nalf+Nbonds))
	!		print *,(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi
	  !	do k=1,Nx
	  !  PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	 !   enddo
	!	  call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,Lmax,EPS,EPF,KEY,ghom)
	!	 DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
    !  PAR(K)=PARalf(K)
	!  enddo  
	!  PARXloc(1:(Nx+Nalf+Nbonds),Lmax)=PAR(1:(Nx+Nalf+Nbonds))
	   	print *,'Xloc='
	pN(L)=Pnormal(L,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
			print *,'Xloc='
	pN(Lmax)=Pnormal(Lmax,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom) 
	PAR(1:(Nx+Nalf+Nbonds))=PARXloc(1:(Nx+Nalf+Nbonds),L)                                                 
      CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)                  
	 ! pN(L)=Pnormal(L,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
	 ! CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)
	  HH=0
	  endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!	  	if (dabs(FUtau).gt.1.0d-12) HH=1   
 !Print *, 'pN=',pN(L),TL(L),pN(Lmax)
	if (HH.eq.1) then
5	kk=kk+1
Print *, 'pN=',pN(L),TL(L) ,pN(Lmax)
	F(1)=FUtau
	X=PAR
	Xlmax(1:Nx)=PARXloc(1:Nx,Lmax)
	X(NK)=TL(L)
	Xlmax(NK)=TL(Lmax)
	Phi=JX
	vol1(L)=volume(L)
	vol1(Lmax)=volume(Lmax)
	DPx(NK)=DDPx(NK)                                           
	  TL(L)=X(NK)+X(NK)*DPx(NK)
	vol2(L)=volume(L)

	do k=1,Nx
	PAR(k)=X(k)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	enddo
	     call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
		 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
		if ((mainK.eq.0).and.(switch.eq.0))	   then
	!	TL(Lmax)=Xlmax(NK)-DPx(NK)*X(NK)
		vol2(Lmax)=volume(Lmax)
	!		print *,(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi
		do k=1,Nx
	    PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	    enddo
	     call check_borderX(PAR,Lmax,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
		 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,Lmax,switch)
		endif
	 kkk=0
	 do while ((mainK.ne.0).or.(switch.ne.0))
	 kkk=kkk+1
	 DPx(NK)=DPx(NK)*0.9d0
	 	if (dabs(DPx(NK)).lt.1.0d-12) then	!-16
		DPx(NK)=-DPx(NK)*1.0d10
		endif 
	  TL(L)=X(NK)+X(NK)*DPx(NK)
	vol2(L)=volume(L)
	do k=1,Nx
	PAR(k)=X(k)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!*(1.0d0-(vol2(L)-vol1(L))/vol1(L)/Phi)
	enddo
	  call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	   call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
	   if ((mainK.eq.0).and.(switch.eq.0))	   then
	!	TL(Lmax)=Xlmax(NK)-DPx(NK)*X(NK)
		vol2(Lmax)=volume(Lmax)
	!		print *,(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi
		do k=1,Nx
	    PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	    enddo
	     call check_borderX(PAR,Lmax,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
		 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,Lmax,switch)
		endif
	 enddo
	 do k=1,Nx
	PAR(k)=X(k)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!*(1.0d0-(vol2(L)-vol1(L))/vol1(L)/Phi)
	enddo
	   call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo
	  PARXloc(1:(Nx+Nalf+Nbonds),L)=PAR(1:(Nx+Nalf+Nbonds))
	!		print *,(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi
	  	do k=1,Nx
	    PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	    enddo
		  call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,Lmax,EPS,EPF,KEY,ghom)
		 DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo  
	  PARXloc(1:(Nx+Nalf+Nbonds),Lmax)=PAR(1:(Nx+Nalf+Nbonds))   
	pN(L)=Pnormal(L,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
	pN(Lmax)=Pnormal(Lmax,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom) 
	PAR(1:(Nx+Nalf+Nbonds))=PARXloc(1:(Nx+Nalf+Nbonds),L)                                                 
      CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)                                    
    !  print *,(FUtau-F(1)),DPx(NK)             
	  DD(1,1)=(FUtau-F(1))/(TL(L)-X(NK)) !DPx(NK)  
	PAR=X
	PARXloc(1:Nx,Lmax)=Xlmax(1:Nx)                                                           
     TL(L)=X(NK) 
	 TL(Lmax)=Xlmax(NK)
	  PARR(1)=0.0d0                    
	   FF(1)=-F(1)
	   if (dabs(FUtau-F(1)).le.1.0d-100) goto 7
!	 write(*,*)'!!!!!!!!!!!!!! NEWTalf after dx ended',DPx(1:3) !,DD(1:3,1:3)
!	 write(*,*) ((DD(i,j), J=1,3),I=1,3) 																					 
!	   print *,((DD(i,j),j=1,Nx),i=1,Nx)
	  CALL DLSLRG(1, DD, 150, FF, 1, PARR)
 TL(L)=X(NK)+PARR(1)
vol2=volume(L)
	do k=1,Nx
	PAR(k)=X(k)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	enddo
	     call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
	 if ((mainK.eq.0).and.(switch.eq.0))	   then
	!	TL(Lmax)=Xlmax(NK)-PARR(1)
		vol2(Lmax)=volume(Lmax)
		do k=1,Nx
	    PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	    enddo
	     call check_borderX(PAR,Lmax,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
		 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,Lmax,switch)
		endif
	 kkk=0
	 do while ((mainK.ne.0).or.(switch.ne.0))
	 kkk=kkk+1
	 PARR(1)=PARR(1)*0.8d0	  
	  TL(L)=X(NK)+PARR(1)
	vol2(L)=volume(L)
	do k=1,Nx
	PAR(k)=X(k)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	enddo
	  call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	   call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,L,switch)
	   if ((mainK.eq.0).and.(switch.eq.0))	   then
	!	TL(Lmax)=Xlmax(NK)-PARR(1)
		vol2(Lmax)=volume(Lmax)
		do k=1,Nx
	    PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	    enddo
	     call check_borderX(PAR,Lmax,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
		 call minmaxALF_fromTL(PAR,Nx,Nalf,Nalf,Lmax,switch)
		endif
	 enddo
	 	do k=1,Nx
	PAR(k)=X(k)/(1.0d0+(vol2(L)-vol1(L))/vol1(L)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	enddo
	   call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo
	  PARXloc(1:(Nx+Nalf+Nbonds),L)=PAR(1:(Nx+Nalf+Nbonds))
	  	do k=1,Nx
	    PAR(k)=PARXloc(k,Lmax)/(1.0d0+(vol2(Lmax)-vol1(Lmax))/vol1(Lmax)/Phi) !!!*(1.0d0-(vol2-vol1)/vol1/Phi)
	    enddo
	    call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,Lmax,EPS,EPF,KEY,ghom)
		DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo
	   PARXloc(1:(Nx+Nalf+Nbonds),Lmax)=PAR(1:(Nx+Nalf+Nbonds))
	pN(L)=Pnormal(L,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
	pN(Lmax)=Pnormal(Lmax,PARXloc,Nx,Nalf,DPalf,EPS,EPF,ghom)
	PAR(1:(Nx+Nalf+Nbonds))=PARXloc(1:(Nx+Nalf+Nbonds),L)    
    CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)  
	  A=DABS(PARR(1))
	!  print *,TL(L),PARR(1),FUtau
	  IF(KK.GT.1000) write(*,*)'>1000 iterations in NEWT_mechanical',A                                                     
      IF(A.GT.EPS*1.0d0)  GOTO 5                                                                                                                                                              
	A=DABS(FUtau)
 
   write(*,*)'A=sum of FUmech =',A
   endif                                          
      KEY=1                                                             
   !   IF(A.GT.1.0d-4)  KEY=-1.0d0
	  
!	  if (KEY.lt.0.0d0) GOTO 5
                                     
   7  CONTINUE                    
                           
      RETURN                                                            
      END                                 
!********************************************************************************************************************
     SUBROUTINE NEWTx(Nx,Nalf,PARx,PARalf,DDPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)	   !searching Xloc in layer L
	 IMPLICIT NONE
	  CHARACTER*3 shape
	  character*7 name_seg(150)
	   integer Nx,Nalf,KEY,L,I,J,KK,kkk,K,HH,chem(150),type_mol(150),phase,rank_seg(150),Lmax,D,iter,final,final1
	  REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150)
	  REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
      REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150),Falfeq,PARalfeq(150),Fxeq,PARxeq(150)
	  REAL*8 TL(150),EPS,EPF,A,DDalf(150,150)                                           
	  REAL*8 PAR(150),PARR(150),FU(150),DD(150,150),X(150),F(150),FF(150),LAM,PARx(150),PARalf(150),DPx(150),DDPx(150),DPalf(150),ghom(150),PAR1(150)
	  integer NK,NG,Nmol,Nbonds,allK,KEYx(100),KEYxend(100),KEYxL(100),KEYid(100),id(100),mainK,KEY_step 
	  REAL*8 ZZ,ZR,JX,IJ,TLmin,minn,FUtau,ALlocREM(10,150),teta_chemREM(150,150),QTrem,QT,RT,XlocL(150)
	  common /layer/shape,TL,Lmax,TLmin
      COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
      COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  common /layer_one/PAR1
	  common /ph/phase
	   common /chem_or/chem_up,ALchemold
	   common /chemdwn/chem_dwn,chem_al	
	   common /cycles/iter,final,final1
	    common /eq/Falfeq,PARalfeq,Fxeq,PARxeq
	   common /for_steps/ALlocREM,teta_chemREM,QTrem
	!	 common /start_of_NEWTX/iter_again
	!	 if (iter.gt.2) EPF=5.0d-2                                           
      KEY=0                                                             
      KK=0                                                              
      LAM=1.0d0 
	  HH=0
	  minn=100.0d0 
!	 print *,'start of NEWTx'	                                                    !Nx=NK-1, Nalf=2*(NG-1)
	  DO  I=1,Nx                                                       
      PAR(I)=PARx(I)
	  if (PARx(I).lt.minn) minn=PARx(I)
	  enddo
	
!	    if ((L.eq.17).and.(iter.gt.0)) then
!	print *,'NEWTx'
!	endif
	   call check_borderX(PAR,L,Lmax,Nx,id,mainK,0)	   !0= correction of initial PAR	

	 	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	 write(*,*)'!!!!!!!!!!!!!! first NEWTalf ended'  
!	 stop 
	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(I)=PARalf(I)
	  enddo
	   CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)
	    DO D=1,Nx
	   if (dabs(FU(D)).gt.1.0d-12) HH=1
	   enddo
! iter_again=1 
!	    if ((L.le.18).and.(L.ge.17).and.(iter.lt.200)) HH=0	  !((L.le.1).and.(iter.lt.100)) 
	   if (L.eq.Lmax) HH=0
!	   if (L.eq.1) HH=0
!	  if (L.le.8) HH=0

     if (HH /= 0) then                                                                                                                                                   
   5  DO  I=1,Nx                                                        
      F(I)=FU(I)                                            
      X(I)=PAR(I)
	  enddo
	  DO I=Nx+1,Nx+Nalf+Nbonds
	  X(I)=PAR(I)
	  enddo
	  DO  I=1,NK 
	  XlocL(I)=Xloc(I,L)
	  enddo 
	  call total(XlocL,QT,RT)
	  ALlocREM=ALloc
	  teta_chemREM=teta_chem
	  QTrem=QT                                                     
      KK=KK+1
!	  print *,PAR(5:6),AL_up(2,L),AL_al(2,L),chem_al(1,L),chem_al(2,L)                                                          
      DO I=1,Nx
	  DPx(I)=DDPx(I) !*PAR(I)                                                    
      PAR(I)=PAR(I)+DPx(I) !*X(I)							!Xloc(I,L)
	 do J=Nx+1,Nx+Nalf+Nbonds
	  PAR(J)=X(J)
	 enddo
	  call minmaxX(PAR,L,Lmax,Nx,allK,KEYx,KEYxL,KEYxend,KEYid,-I)
	!  print *,'!'
	if (KEYid(I).eq.0) call check_ALFA_step(PAR,L,Nx,Nalf,X,DPx,I,KEY_step)
	if (KEYid(I).eq.0) KEYid(I)=KEY_step
	 ! if (I.eq.2) print *,KEYx(I),KEYxL(I),KEYxend(I),KEYid(I)
	  kkk=0
	 do while (KEYid(I).ne.0)
	 kkk=kkk+1
			if (kkk.gt.1000) then
			print *,'many steps of DP in NEWTx',DPx(I),I,'keyid=',KEYid(I)
			do J=1,Nx
			PAR(J)=X(J)
			enddo
			goto 99
			!stop
			endif
	 PAR(I)=PAR(I)-DPx(I) !*X(I)
		if (dabs(DPx(I)).lt.1.0d-12*X(I)) then	!-16
		DPx(I)=-DPx(I)*1.0d10
		endif 

	if  (KEYid(I).eq.1) then
		if (DPx(I).gt.0.0d0) DPx(I)=DPx(I)*1.1d0
		if (DPx(I).lt.0.0d0) DPx(I)=DPx(I)*0.9d0
	else if  (KEYid(I).eq.(-1)) then
		if (DPx(I).gt.0.0d0) DPx(I)=DPx(I)*0.9d0
		if (DPx(I).lt.0.0d0) DPx(I)=DPx(I)*1.1d0
	endif
	
	PAR(I)=PAR(I)+DPx(I) !*X(I)
	call minmaxX(PAR,L,Lmax,Nx,allK,KEYx,KEYxL,KEYxend,KEYid,-I)
	!!!!!!!
	if (kkk.lt.150) then
	if (KEYid(I).eq.0)	call check_ALFA_step(PAR,L,Nx,Nalf,X,DPx,I,KEY_step)
	if (KEYid(I).eq.0) KEYid(I)=KEY_step
	endif
!		  if (I.eq.2) print *,KEYx(I),KEYxL(I),KEYxend(I),KEYid(I)
	enddo	!while
!	 print *,'after DPx',I,PAR(I),L
!	  print *,I !,PAR(5:6),AL_up(2,L),AL_al(2,L),chem_al(1,L),chem_al(2,L) 
 if (iter.gt.4) then 
	if (kkk.ge.150)	 	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	  	  print *,I,PAR(5:6),AL_up(2,L),AL_al(2,L),chem_al(1,L),chem_al(2,L)
	if (kkk.lt.150) call ALFA_step (Nx,Nalf,X,PARalf,DPx,I,L)
	else
 	   call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	endif
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	!  print *,I,DPx(I),PARalf(K)
	  enddo
                                                       
      CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)
                                            
      DO  J=1,Nx                   
	  DD(J,I)=(FU(J)-F(J))/DPx(I) !*X(I) !(dlog(PAR(I)/X(I)) )  !DPx(I) 		! DD(J,I) changed!!!
	  enddo                  !cycle J                                                             
      PAR(I)=PAR(I)-DPx(I) !*X(I)
	  enddo  				   !cycle I 
 
	  do I=1,Nx
	  PARR(I)=0.0d0 !PAR(I)
	  enddo
	   DO  J=1,Nx                   
	   FF(J)=-F(J)
	   enddo                  !cycle J 
  	 !    if ((L.eq.2).and.(iter.eq.13)) then 
!	 write(*,*)'!!!!!!!!!!!!!! NEWTalf after dx ended',DPx(1:3) !,DD(1:3,1:3)
!	 write(*,*) ((DD(i,j), J=1,3),I=1,3) 
	 !endif  																					 
!	   print *,((DD(i,j),j=1,Nx),i=1,Nx)
!		  stop

	  CALL DLSLRG(Nx, DD, 150, FF, 1, PARR)                                          
	  ! !X(1)before=',X(1),X(2),X(3),'X(1)after=',X(1)+PARR(1),X(2)+PARR(2),X(3)+PARR(3)                                             
!	if ((L.eq.17).and.(iter.gt.0)) then
!	print *,'NEWTx'
!	endif

	 DO I=1,Nx
	  	PARR(I)=PARR(I) !*X(I)		   
		 PAR(I)=X(I)+PARR(I)
	 enddo
	   call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	 if (mainK.eq.0)	call check_ALFA_step(PAR,L,Nx,Nalf,X,PARR,0,KEY_step)
	if (mainK.eq.0) mainK=KEY_step
	  kkk=0
	 do while (mainK.ne.0)
	 kkk=kkk+1
	 DO J=1,Nx
	! if (id(J).ne.0) then
	!	if ((PARR(J)*id(J)).gt.0.0d0) then 
	!	PARR(J)=PARR(J)*1.1d0
	!	 else
		PARR(J)=PARR(J)*0.8d0	
	!	endif
	 PAR(J)=X(J)+PARR(J)
			if (kkk.gt.10000) then
			print *,'many steps after iter in NEWTx',id,PAR
			!stop
				do I=1,Nx
			PAR(I)=X(I)
			enddo
			goto 99
			endif
	! endif
	enddo	!J
	 call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)           !-1 = checking   PAR
	 if (kkk.lt.150) then
	if (mainK.eq.0)	call check_ALFA_step(PAR,L,Nx,Nalf,X,PARR,0,KEY_step)
	if (mainK.eq.0) mainK=KEY_step
	endif

	enddo		 !while
!	 write(*,*)'!!!after dlslrgX'
 if (L.eq.10) then
 A=0.0d0
 endif
    A=0.0d0 
	  DO I=1,Nx 
	  A=A+DABS(PARR(I))
      PAR(I)=X(I)+PARR(I)
	  enddo
	  if (iter.gt.4) then
if (kkk.lt.150)   call ALFA_step (Nx,Nalf,X,PARalf,PARR,0,L)
if (kkk.ge.150)	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
else
call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
endif	  
	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(I)=PARalf(I)
	  enddo         
	!   write(*,*)'x=',kkk
	  CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom) 
 
	  write(*,*)'!!!Layer = ',L
	                                            
      A=A/Nx  
	  IF(KK.GT.1000) write(*,*)'>1000 iterations in NEWTx',A
	 ! if (iter.ne.1) then 
	 ! do k=1,Nx
	!  if (PARR(k).gt.(X(k)*1.0d-5)) GOTO 5
	!  enddo
	!  else                                                          
      IF(A.GT.EPS)  GOTO 5
	!  endif                                             
      A=0.0d0
!	  if (Nfix) then
 !  do I=1,Nfix
 !   A=A+DABS(FU(Nx+No_or(I)))+DABS(FU(Nx+Nalf/2+Nbonds+No_or(I)))
 !  enddo 
!   endif                                                               
      DO 6 I=1,Nx 
!	  IF(KK.GT.180) write(*,*)'200',DABS(FU(I))                                                     
   6  A=A+DABS(FU(I))
 
   write(*,*)'A=sum of FU =',A                                          
      KEY=1                                                             
      IF(A.GT.EPF)  KEY=-1.0d0
	  
	  if (KEY.lt.0.0d0) GOTO 5
	  	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	  
	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(I)=PARalf(I)
	  enddo         
	  CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom) 

  endif          !HH  	
99	  DO  I=1,Nx                                                      
      PARx(I)=PAR(I)
	  enddo                                       
   7  CONTINUE                    
                           
      RETURN                                                            
      END                                        
!********************************************************************************************************************
 SUBROUTINE ALFA_step (Nx,Nalf,PARth,PARalf,DPx,k,L)
 IMPLICIT NONE
 REAL*8	 m_up,PARth(150),PARalf(150),DPx(150)
 REAL*8 QT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
 REAL*8 ZZ,ZR,JX,IJ
 integer s,k,Nx,Nalf,NK,NG,Nmol,Nbonds,L,seg_ch,M,I, seg_max_ch,chem(150),type_mol(150),rank_seg(150)
 COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
 COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
  common /for_steps/ALloc,teta_chem,QTrem
REAL*8 ALloc(10,150), teta_chem(150,150),QTrem
 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid 
QT=QTrem
 if (k.ne.0) then
QT=QT+(Q(k)-Q(NK))*DPx(k)
!print *,QT,R(NK),NK	  
		 do s=1,Nalf/2
	   	 PARalf(Nx+s)=PARth(Nx+s)+(QS(k,s)-QS(NK,s)*R(k)/R(NK))/QT*DPx(k)-ALloc(s,L)*(Q(k)-Q(NK)*R(k)/R(NK))/QT*DPx(k)					
		if (L.ne.1)	 PARalf(Nx+Nalf/2+Nbonds+s)=PARth(Nx+Nalf/2+Nbonds+s)+(QS(k,s)-QS(NK,s)*R(k)/R(NK))/QT*DPx(k)-ALloc(s,L)*(Q(k)-Q(NK)*R(k)/R(NK))/QT*DPx(k)		
		if ((L.eq.1).and.(r_solid.lt.1.0d-50))  PARalf(Nx+Nalf/2+Nbonds+s)=0.0d0
		 enddo	
		 do M=1,Nbonds
		  PARalf(Nx+Nalf/2+M)=PARth(Nx+Nalf/2+M)
	if (k.eq.seg_ch(M))	 PARalf(Nx+Nalf/2+M)=PARth(Nx+Nalf/2+M)+DPx(k)/ZZ/QT-teta_chem(M,L)*(Q(k)-Q(NK)*R(k)/R(NK))/QT/ZZ*DPx(k)
	if (k.eq.(seg_max_ch(M)))	 PARalf(Nx+Nalf/2+M)=PARth(Nx+Nalf/2+M)+(1.0d0-teta_chem(M,L)/PARth(seg_ch(M))*PARth(k)*(Q(k)-Q(NK)*R(k)/R(NK))-1.0d0/m_up(L))*DPx(k)/ZZ/QT					
		 enddo 
 else
	   do s=1,Nalf/2
	   PARalf(Nx+s)=PARth(Nx+s)					
	   if (L.ne.1)	 PARalf(Nx+Nalf/2+Nbonds+s)=PARth(Nx+Nalf/2+Nbonds+s)
	   if ((L.eq.1).and.(r_solid.lt.1.0d-50))  PARalf(Nx+Nalf/2+Nbonds+s)=0.0d0
	   enddo	
		 do M=1,Nbonds
		 PARalf(Nx+Nalf/2+M)=PARth(Nx+Nalf/2+M)					
		 enddo
do I=1,Nx
QT=QT+(Q(I)-Q(NK))*DPx(I)
enddo 
	do I=1,Nx
	 do s=1,Nalf/2
	   	 PARalf(Nx+s)=PARalf(Nx+s)+(QS(I,s)-QS(NK,s)*R(I)/R(NK))/QT*DPx(I)-ALloc(s,L)*(Q(I)-Q(NK)*R(I)/R(NK))/QT*DPx(I)					
		if (L.ne.1)	 PARalf(Nx+Nalf/2+Nbonds+s)=PARalf(Nx+Nalf/2+Nbonds+s)+(QS(I,s)-QS(NK,s)*R(I)/R(NK))/QT*DPx(I)-ALloc(s,L)*(Q(I)-Q(NK)*R(I)/R(NK))/QT*DPx(I)		
		if ((L.eq.1).and.(r_solid.lt.1.0d-50))  PARalf(Nx+Nalf/2+Nbonds+s)=0.0d0
		 enddo	
		 do M=1,Nbonds
		if (I.eq.seg_ch(M))	 PARalf(Nx+Nalf/2+M)=PARalf(Nx+Nalf/2+M)+DPx(I)/ZZ/QT-teta_chem(M,L)*(Q(I)-Q(NK)*R(I)/R(NK))/QT/ZZ*DPx(I)
	if (I.eq.(seg_max_ch(M)))	 PARalf(Nx+Nalf/2+M)=PARalf(Nx+Nalf/2+M)+(1.0d0-teta_chem(M,L)/PARth(seg_ch(M))*PARth(I)*(Q(I)-Q(NK)*R(I)/R(NK))-1.0d0/m_up(L))*DPx(I)/ZZ/QT				
		 enddo 
	enddo	 !I
 endif
       RETURN                                                            
      END                                        
!********************************************************************************************************************
SUBROUTINE check_ALFA_step (PARth,L,Nx,Nalf,X,DPx,I,KEYid)
 IMPLICIT NONE
 REAL*8	 m_up,PARth(150),PARalf(150),DPx(150),ALlocnew(150),sumXX(150),X(150)
 REAL*8 T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
 REAL*8 ZZ,ZR,JX,IJ
 integer s,k,Nx,Nalf,NK,NG,Nmol,Nbonds,L,seg_ch,M,I,KEYid,id(100),mainK
 COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
! COMMON /C/ QT,RT,T,TT,R,RS,Q,QS,vol_ref, q_ref
 COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
 ! print *,R(1:NK)

  call calcALloc(PARth,L,ALlocnew,sumXX)
  !print *,R(NK)

	call ALFA_step (Nx,Nalf,X,PARalf,DPx,I,L)
	 DO  s=Nx+1,Nx+Nalf+Nbonds                                                       
     PARth(s)=PARalf(s)
	 enddo
	  do s=1,Nbonds
	!	 PARalf(Nx+Nalf/2+s)=sumXX(s)/(dexp(-2.0d0*Wprev(s))-1.0d0)
	  enddo 
	call check_borderALF(PARth,L,Nalf,ALlocnew,sumXX,id,mainK,-1)
	if (mainK.ne.0) then
		if (I.eq.0) then
		KEYid=1
		else
	if (DPx(I).gt.0.0d0)  KEYid=-1
	if (DPx(I).lt.0.0d0)  KEYid=1
		endif
	else
	KEYid=0
	endif
       RETURN                                                            
      END                                        
!********************************************************************************************************************
      SUBROUTINE NEWTalf_chem(Nx,Nalf,PARth,PARalf,DPalfth,FUth,DDalf,L,EPS,EPF,KEY,ghom)
	  IMPLICIT NONE
	  character*7 name_seg(150)
	  integer KEYid(100),id(100),KEYch(100),mainK,sup,chemup,sdwn,allK
	  integer M,D,S,Nx,Nalf,KEY,L,I,J,KK,kkk,chem(150),type_mol(150),CHANGE(150),Nsolve,rank_seg(150),HH
	  REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150)
	  REAL*8 ALchemold(150,150),chem_upold(150,150),FUtau
      REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALlocnew(150),sumXX(150)
	  REAL*8 EPS,EPF,A,ZZ,ZR,JX,IJ,Falfeq,PARalfeq(150),Fxeq,PARxeq(150)                                         
      REAL*8 PAR(150),PARth(150),PARR(150),FU(150),FUth(150),DDalf(150,150),X(150),F(150),FF(150),LAM,PARalf(150),DPalf(150),DPalfth(150),ghom(150),alf_max(150),alf_min(150)   
      COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /without_or/ Nfix,No_or,CHANGE
	  common /chem_or/chem_up,ALchemold
	  common /chemdwn/chem_dwn,chem_al
	  common /oldch/chem_upold
	  common /cycles/iter,final,final1
	  common /eq/Falfeq,PARalfeq,Fxeq,PARxeq
	  common /old/Xold,AL_upold,AL_dwnold,AL_alold,ALlocold
	  common /min_maxALF/alf_max,alf_min
	 ! common /start_of_NEWTX/iter_again
	  REAL*8 Xold(150,150),AL_upold(10,150),AL_dwnold(10,150),AL_alold(10,150),ALlocold(10,150)	
	  integer NK,NG,Nmol,Nfix,Nbonds,iter,final,final1,No_or(150)
	   REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
		COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid       
	 ! print *,'start of NEWTalf'
	!  if (iter.gt.20) EPF=5.0d-4 
	  HH=0
	  KK=0                         !number of iterations in NEWTalf_chem                                                            
      LAM=1.0d0 
	  Nsolve=Nalf-Nfix*2  		   !number of alf_up+alf_dwn to find  (without chem)
	   
	  call calcALloc(PARth,L,ALlocnew,sumXX)
	  do S=1,NG 
	  ALloc(s,L)=ALlocnew(s)
	  enddo 
	    do s=1,Nalf/2
	   	 PARalf(Nx+s)=AL_up(s,L)					
		if (L.ne.1)	 PARalf(Nx+Nalf/2+Nbonds+s)=AL_dwn(s,L)
		if ((L.eq.1).and.(r_solid.lt.1.0d-50))  PARalf(Nx+Nalf/2+Nbonds+s)=0.0d0
		 enddo	
		 do M=1,Nbonds
		 teta_chem(M,L)=ZZ*ALlocnew(M+NG)
		 PARalf(Nx+Nalf/2+M)=chem_up(M,L)
	!	 if (teta_chem(M,L).gt.(2.5d0*teta_chem(M,L+1))) PARalf(Nx+Nalf/2+M)=chem_dwn(M,L+1)-sumXX(M) !ALlocnew(M+NG)*chem_upold(M,L)/ALchemold(M,L)
	!	 if (teta_chem(M,L).gt.teta_chem(M,L+1)) then
	!		if ((teta_chem(M,L+1)/ZZ-sumXX(M)).gt.0.0d0) PARalf(Nx+Nalf/2+M)=(teta_chem(M,L+1)/ZZ-sumXX(M))*chem_upold(M,L)/ALchemold(M,L)
	!		if ((teta_chem(M,L+1)/ZZ-sumXX(M)).le.0.0d0) PARalf(Nx+Nalf/2+M)=sumXX(M)*0.1d0*chem_upold(M,L)/ALchemold(M,L)
	!	 endif					
		 enddo     
	   
	  DO  I=Nx+1,Nx+Nalf+Nbonds				  !Nx=NK-1, Nalf=2*(NG-1)
	  PARth(I)=PARalf(I)					!creating unified array PARth
	  enddo 
	!   if (L.eq.17) stop

 call check_borderALF(PARth,L,Nsolve,ALlocnew,sumXX,id,mainK,0)	!0 = corrections   PARth	 	

!	 if (L.eq.19) stop
	  DO D=1,Nsolve+Nbonds
	  PAR(D)=PARth(CHANGE(D))
	  DPalf(D)=DPalfth(CHANGE(D))
	  enddo
	  !	 print *,'ffff'

	 CALL inhomsol(Nx,Nalf,PARth,FUth,L,FUtau,ghom)			 !calc FU for initial point
	 !	 print *,'ffff'

	 ! if (L.eq.12) print *,'12!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'	    
  	  DO D=1,Nsolve+Nbonds
	  FU(D)=FUth(CHANGE(D))
	  if (dabs(FU(D)).gt.1.0d-12) HH=1
	HH=1
	  enddo

 if (HH /= 0) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                         
   5  DO  I=1,Nsolve+Nbonds                                                        
      F(I)=FU(I)                                                                                                      
      X(I)=PAR(I) !/PARalf(CHANGE(I))
	  enddo
	  KK=KK+1                                                                                                               
      DO I=1,Nsolve+Nbonds
	!	if (L.eq.12)  print *,'!',I,Par(I),DPalf(I)
	   if  ((I.le.(Nsolve/2+Nbonds)).and.(PAR(I).eq.(alf_max(I)-1.0d-150)).and.(DPalf(I).gt.0.0d0)) DPalf(I)=-DPalf(I)
		  if  ((I.gt.(Nsolve/2+Nbonds)).and.(PAR(I).eq.(alf_max(I)-1.0d-150)).and.(DPalf(I).gt.0.0d0)) DPalf(I)=-DPalf(I)
	      if  ((PAR(I).le.(alf_min(I)+1.0d-150)).and.(DPalf(I).lt.0.0d0)) DPalf(I)=-DPalf(I)
		  if  ((PAR(I).ge.(alf_max(I)-1.0d-150)).and.(DPalf(I).gt.0.0d0)) DPalf(I)=-DPalf(I)
		  if ( (dabs(DPalf(I))).gt.(alf_max(I)-alf_min(I)) )  DPalf(I)=(alf_max(I)-alf_min(I))*0.1d0 !/PARalf(CHANGE(I))
	    PAR(I)=PAR(I)+DPalf(I) !*PARalf(CHANGE(I))						 !alfa up and down
		PARth(CHANGE(I))=PAR(I)
		call minmaxALF(PARth,Nx,Nalf,Nsolve,ALlocnew,sumXX,L,KEYid,KEYch,sup,chemup,sdwn,allK,-I)
	kkk=0
	do while (KEYid(I).ne.0)
	kkk=kkk+1
	!	if (L.eq.14)  print *,'!',I,Par(I),DPalf(I)

			if (kkk.gt.1000) then
			print *,'many steps of DP in NEWTalf',DPalf(I),I
			DO J=1,Nsolve+Nbonds
			PARth(CHANGE(J))=X(J)
			enddo
			goto 88
			!stop
			endif

   	PAR(I)=PAR(I)-DPalf(I) !*PARalf(CHANGE(I))
		if (dabs(DPalf(I)).lt.1.0d-12*X(I)) then				 !-16
		DPalf(I)=-DPalf(I)*1.0d10
		endif 
	if  (KEYid(I).eq.1) then
		if (DPalf(I).gt.0.0d0) DPalf(I)=DPalf(I)*1.1d0
		if (DPalf(I).lt.0.0d0) DPalf(I)=DPalf(I)*0.9d0
	else if  (KEYid(I).eq.(-1)) then
		if (DPalf(I).gt.0.0d0) DPalf(I)=DPalf(I)*0.9d0
		if (DPalf(I).lt.0.0d0) DPalf(I)=DPalf(I)*1.1d0
	endif
	PAR(I)=PAR(I)+DPalf(I) !*PARalf(CHANGE(I))
	PARth(CHANGE(I))=PAR(I)
	call minmaxALF(PARth,Nx,Nalf,Nsolve,ALlocnew,sumXX,L,KEYid,KEYch,sup,chemup,sdwn,allK,-I)
	enddo	!while
	!print *,'ggg'                             
      CALL inhomsol(Nx,Nalf,PARth,FUth,L,FUtau,ghom) 
	  
	  DO D=1,Nsolve+Nbonds
	  FU(D)=FUth(CHANGE(D))
	  enddo
                                          
       DO  J=1,Nsolve+Nbonds                   
	   DDalf(J,I)=(FU(J)-F(J))/DPalf(I) !*X(I) !*PARalf(CHANGE(I)) !(dlog(PAR(I)/X(I)) ) !DPalf(I) 
	   enddo                  !cycle J                     

	  PAR(I)=X(I) !PAR(I)-DPalf(I) !*PARalf(CHANGE(I))
	  PARth(CHANGE(I))=PAR(I)
	  if ((L.eq.1).and.(I.eq.(Nsolve/2+Nbonds)).and.(r_solid.lt.1.0d-50))  GOTO 150
	enddo                    !cycle I  
	150  do I=1,Nsolve+Nbonds
	  PARR(I)=0.0d0 
	  enddo
	   DO  J=1,Nsolve+Nbonds                   
	   FF(J)=-F(J)*LAM
	   enddo                  !cycle J    
	   	 
if ((L.eq.1).and.(r_solid.lt.1.0d-50)) then
	 CALL DLSLRG(Nsolve/2+Nbonds, DDalf, 150, FF, 1, PARR)  
else
  CALL DLSLRG(Nsolve+Nbonds, DDalf, 150, FF, 1, PARR)           !PARR - steps after solution of system's equations                                 
  endif                            
      DO I=1,Nsolve+Nbonds
	  if ((L.eq.1).and.(I.gt.(Nsolve/2+Nbonds)).and.(r_solid.lt.1.0d-50)) then 			   
		 PARR(I)=0.0d0
		 endif
		! X(I)=X(I)*PARalf(CHANGE(I))
		 PARR(I)=PARR(I) !*X(I) !*PARalf(CHANGE(I))
		  PARth(CHANGE(I))=X(I)+PARR(I)
	   enddo !I
	  
	   call check_borderALF(PARth,L,Nsolve,ALlocnew,sumXX,id,mainK,-1)	!-1 = checking   PARth; id and mainK - output
	  
	 kkk=0
	 do while (mainK.ne.0)
	 kkk=kkk+1
	 DO J=1,Nsolve+Nbonds
!	 if (PARth(CHANGE(J)).le.(alf_min(J)+1.0d-150)) then
!	 PARR(I)=PARR(J)+0.05d0*(alf_max(J)-alf_min(J))
!	 elseif (PARth(CHANGE(J)).ge.(alf_max(J)-1.0d-150)) then
!	 PARR(I)=PARR(J)-0.05d0*(alf_max(J)-alf_min(J))
!	 else
	! if (id(J).ne.0) then
	!	if ((PARR(J)*id(J)).gt.0.0d0) then 
	!	PARR(J)=PARR(J)*1.1d0
	!	 else
		PARR(J)=PARR(J)*0.99d0	
!		endif
	PARth(CHANGE(J))=X(J)+PARR(J)
			if (kkk.gt.10000) then
			print *,'many steps after iter in NEWTalf',id,PARth
			!stop
			DO I=1,Nsolve+Nbonds
			PARth(CHANGE(I))=X(I)
			enddo
			goto 88
			endif
	! endif
	enddo	!J
	call check_borderALF(PARth,L,Nsolve,ALlocnew,sumXX,id,mainK,-1)	!-1 = checking   PARth
	enddo		 !while
                       
      CALL inhomsol(Nx,Nalf,PARth,FUth,L,FUtau,ghom) 
	  DO D=1,Nsolve+Nbonds
	  FU(D)=FUth(CHANGE(D))
	  !if (dabs(FU(d)).gt.dabs(F(d))) PARR(D)=0.0d0
	  enddo

      A=0.0d0 
	  DO I=1,Nsolve+Nbonds
	  A=A+DABS(PARR(I))
      PAR(I)=X(I)+PARR(I)
	  enddo                                 
      A=A/(Nsolve+Nbonds) 
	  !print *,A,EPS !A,FU(1:(Nsolve+Nbonds))   
	   IF(KK.GT.1000) write(*,*)'>1000 iterations in NEWTalf',A
	 !  if (iter.ne.1) then 
	 ! do I=1,Nsolve+Nbonds
	  !if (PARR(I).gt.(X(I)*1.0d-4)) GOTO 5
	 ! enddo
	  !else                                                          
      IF(A.GT.EPS)  GOTO 5
	 ! endif                                                      
      A=0.0d0                                                              
      DO 6 I=1,Nsolve+Nbonds 
!	  IF(KK.GT.180) write(*,*)'200',DABS(FU(I))                                                      
   6  A=A+DABS(FU(I))
 !  if (Nfix) then
 !  do I=1,Nfix
 !   A=A+DABS(FUth(Nx+No_or(I)))+DABS(FUth(Nx+Nalf/2+Nbonds+No_or(I)))
  ! enddo 
  ! endif  
   	!  print *,A,EPF
	 KEY=1                                                             
      IF(A.GT.EPF)  KEY=-1

	  if (KEY.lt.0.0d0) GOTO 5
	  endif !!!!!!!!!!!!!!!!!!!!!!!!!!!!HH
	88  DO  I=Nx+1,Nx+Nalf+Nbonds
	  PARalf(I)=PARth(I)					!output values as input for the next iteration in micelle_energy 
	  enddo
                                       
   7  CONTINUE       
  ! write(*,*)'NEWTalf finished'                                              
      RETURN                                                            
      END                                                                         
  !********************************************************************************************************************
     SUBROUTINE LVE(X,NK,X1,X2)  
	 IMPLICIT NONE                                                                                         
      INTEGER I   
	  character*7 name_seg(150)                                                      
      REAL*8 N2,SL,SS,WG,WK,PRES,X1R(150),X2R(150),RE,N1,Xinn(150),T
	  REAL*8 ZZ,ZR,JX,IJ                  
      INTEGER NK,NG,NT,Nmol,NKKK,NKK,Nbonds,rank_seg(150),NNK                                               
      REAL*8 X(150),X1(150),X2(150),Y(150),GAM1(150),GAM2(150),P,PP
	  COMMON /S/ Xtohom,X1e,X2e,XXe,Ye,GAM1e,GAM2e 
	  COMMON /D/ ZZ,ZR,NKKK,NG,JX,IJ,Nmol,Nbonds    
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
      real*8 teta_chem(150,150)
	  integer chem(150),type_mol(150)                                       
	  REAL*8 Xtohom(150),X1e(150),X2e(150),XXe(150),Ye(150),GAM1e(150),GAM2e(150),ghom(150)
     common /equib/Xinn,NNk,NKK,P,PP,T
	                                         
      NT=Nmol
	  DO I=1,NK
	  if (type_mol(I).eq.0) NT=NT+1 
	  enddo
	  NKK=NT
	  NT=NT-1 
	 !print *, X(1:NKK) 
	  X1(NKK)=0.04d0         !X1 = 1,..,Nmol,..others                                              
!     X1(NK)=0.08                                                       
      RE=0.0d0 	  
	                                                             
  !!!    NT=NK-1                                                           
      DO 1 I=1,NT                                                       
      X1(I)=X(I)*(1.0d0-X1(NKK))                                              
      X2(I)=0.4D-05                                                     
!     X2(I)=0.1                                                         
  1   RE=RE+X2(I)   
 ! X1(1)=1.0D-150 
 ! X1(2)=1.0D0-X1(1)                                                   
      X2(NKK)=1.0d0-RE                                                   
      PRES=1.0D-08                                                      
      PRES=1.0D-10 !1.0D-07                                                      
      IJ=0                                                            
      SS=0.0d0                                                            
   	 ! print *, X1(1:NKK) 
	  call XcomptoXseg(X1,Xtohom)	   !X1 - molecular composition, Xtohom - segments' composition
	  call homogsol(ghom,GAM1e)
	   DO I=1,NKK 
	  GAM1(I)=dexp(GAM1e(I))/X1(I)
	  enddo
	! print *, X2(1:NKK) 
	 call XcomptoXseg(X2,Xtohom)	   !X1 - molecular composition, Xtohom - segments' composition
	  call homogsol(ghom,GAM2e)
	  DO I=1,NKK  
	  GAM2(I)=dexp(GAM2e(I))/X2(I)
	  enddo
	  !print *,GAM1e(1),GAM1e(2),X1(1),X1(2)
	  !print *,GAM2e(1),GAM2e(2),X2(1),X2(2)

	                                                
      DO 18 I=1,NKK                                                      
  18  SS=SS+DABS(X1(I)*GAM1(I)-X2(I)*GAM2(I))                           
  19  IJ=IJ+1                                                           
      RE=0.0d0                                                             
      DO  20 I=1,NT                                                     
  20  RE=RE+GAM1(I)*X1(I)/GAM2(I)                                       
      RE=DABS(1-RE)                                                     
      X2(NKK)=RE                                                         
      N2=RE                                                             
      RE=RE*GAM2(NKK)/GAM1(NKK)                                           
      N1=RE                                                             
      DO 21 I=1,NT                                                      
      X1(I)=X(I)*DABS(1.0D0-RE)                                         
      N1=N1+X1(I)                                                       
      X2(I)=X1(I)*GAM1(I)/GAM2(I)                                       
  21  N2=N2+X2(I)                                                       
      X1(NKK)=RE/N1                                                      
      X2(NKK)=X2(NKK)/N2                                                  
      DO 22 I=1,NT                                                      
      X1(I)=X(I)*(1.0d0-X1(NKK))                                             
  22  X2(I)=X2(I)/N2                                                    
      IF (IJ-4) 28,23,23                                                
  23  IF (SS-SL) 24,24,28                                               
  24  IF (SS-0.2d0) 25,25,28                                              
  25  IF(DABS(SS-SL).LT.1.0D-10) GO TO 28            !1.0D-07                   
      WG=SS/(SL-SS)                                                     
      WK=1+WG                                                           
      N1=0.0d0                                                             
      N2=0.0d0                                                             
      DO 30 I=1,NKK                                                      
      X1(I)=WK*X1(I)-WG*X1R(I)                                          
      X2(I)=WK*X2(I)-WG*X2R(I)                                          
      X1(I)=DABS(X1(I))                                                 
      X2(I)=DABS(X2(I))                                                 
      N1=N1+X1(I)                                                       
  30  N2=N2+X2(I)                                                       
      X1(NKK)=X1(NKK)/N1                                                  
      X2(NKK)=X2(NKK)/N2                                                  
      DO  31 I=1,NT                                                     
      X1(I)=X(I)*(1.0d0-X1(NKK))                                             
  31  X2(I)=X2(I)/N2                                                    
  28 call XcomptoXseg(X1,Xtohom)	   !X1 - molecular composition, Xtohom - segments' composition
	  call homogsol(ghom,GAM1e) 
	   DO I=1,NKK 
	  GAM1(I)=dexp(GAM1e(I))/X1(I)
	  enddo

	 call XcomptoXseg(X2,Xtohom)	   !X1 - molecular composition, Xtohom - segments' composition 
	  call homogsol(ghom,GAM2e)
	  DO I=1,NKK  
	  GAM2(I)=dexp(GAM2e(I))/X2(I)
	  enddo
      print *,X1(1:2),X2(1:2)                                        
      SL=SS                                                             
      SS=0.0d0                                                             
      DO 34  I=1,NKK                                                     
  34  SS=SS+DABS(X1(I)*GAM1(I)-X2(I)*GAM2(I))                           
!     PRINT 1234 ,IJ,SS,PRES                                            
!     WRITE  (*,1234) IJ,SS,SL                                          
 1234 FORMAT('I',I3,':',2E16.8)                                         
      IF (SS-PRES) 40,35,35                                             
  35  DO  38  I=1,NKK                                                    
      X1R(I)=X1(I)                                                      
  38  X2R(I)=X2(I)                                                      
      IF(IJ-300) 150,150,51                                               
  150  GO  TO  19                                                        
  51  PRINT 52                                                          
  52  FORMAT (2X,'NUMBER OF ITER MORE THAN 300')                        
  40  RE=0.0d0                                                             
      DO  41  I=1,NT                                                    
  41  RE=RE+X2(I)                                                       
      DO  42  I=1,NT                                                    
  42  Y(I)=X2(I)/RE                                                     
      PP=-8.314462d0*T*GAM2e(NKK-1)/14.244D-06
!-8.3143d0*T*DLOG(GAM1(NKK-1)*X1(NKK-1))/14.244D-06
!	  print *, X1(1:NKK),GAM2e(NKK-1),T
!	  print *, X2(1:NKK)
   IF (IJ.gt.300) PP=1.0d20 !stop                                                                               
      RETURN                                                            
      END                                                               
	  !******************************************************************************************************************** 
	      SUBROUTINE NEWTx_canonic(Nx,Nalf,PARx,PARalf,DDPx,DPalf,FU,DD,L,EPS,EPF,KEY,ghom)
	 IMPLICIT NONE
	  CHARACTER*3 shape
	  character*7 name_seg(150)
	   integer Nx,Nalf,KEY,L,I,J,KK,K,HH,chem(150),type_mol(150),num_seg(150),phase,rank_seg(150),Lmax,D,M,KKEY,iter,final,s,final1
	  REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150)
	  REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150),sumXX(150),ALlocnew(150)
      REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150),Falfeq,PARalfeq(150),Fxeq,PARxeq(150)
	  REAL*8 TL(150),EPS,EPF,A,DDalf(150,150),m_dwn,m_up,lat                                            
	  REAL*8 PAR(150),PARR(150),FU(150),DD(150,150),X(150),F(150),FF(150),LAM,PARx(150),PARalf(150),DPx(150),DDPx(150),DPalf(150),ghom(150),PAR1(150),Wprev
	  integer NK,NG,Nmol,Nbonds 
	  REAL*8 ZZ,ZR,JX,IJ,summa,volume,TLmin,QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref,QTpr,RTpr,QTnext,RTnext,RTx,QTcur,RTcur,Philoc(150,150)
	  real*8 LL,geom(150),DPmin(150),DPmax(150),sumL,sumLnext,sumLprev,PARnext(150),PARalfnext(150),PARprev(150),PARalfprev(150),PARold(150)
	   integer allK,KEYx(100),KEYxend(100),KEYxL(100),KEYid(100),id(100),mainK,KEY_step,kkk,Lmaxx 
	   COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
	  common /layer/shape,TL,Lmaxx,TLmin
      COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
      COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	 common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  common /layer_one/PAR1
	  common /ph/phase
	   common /chem_or/chem_up,ALchemold
	   common /chemdwn/chem_dwn,chem_al	
	   common /cycles/iter,final,final1
	    common /eq/Falfeq,PARalfeq,Fxeq,PARxeq     
      EXTERNAL inhomsol 
!	  EPF=1.0D-5                                                   
      KEY=0                                                             
      KK=0                                                              
      LAM=1.0 
	  HH=1
	  Lmax=Lmaxx+1
!	  if (iter.eq.1) then
!	   DO  I=1,Nx                                                       
  !    PAR(I)=PARx(I)
!	  enddo
!	  	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
  !    PAR(I)=PARalf(I)
!	  enddo
!	   CALL inhomsol(Nx,Nalf,PAR,FU,L,Wprev,ghom)
!	   HH=0
!endif
	 ! if (final1.ne.1) then
	!  if ((L.eq.1).or.(L.eq.2)) HH=0
	  !if (L.le.4) HH=0
	! endif
	!   if (L.eq.Lmax) HH=0
	!   endif
	!  if ((L.eq.2).and.(iter.le.2)) HH=0
	!   if ((L.eq.3).and.(iter.le.10)) HH=0
	if (L.eq.(Lmax)) then
	PAR(NK)=1.0d0
	 DO  I=1,Nx                                                       
      PAR(I)=Xloc(I,L-1)
	  PAR(NK)=PAR(NK)-PAR(I)
	  enddo
	  call total(PAR,QTpr,RTpr) !,NK,Q,R,QT,RT)		!QT, RT - output
	  endif
	  	PAR(NK)=1.0d0
   	 DO  I=1,Nx                                                       
      PAR(I)=Xloc(I,L+1)
	  PAR(NK)=PAR(NK)-PAR(I)
	  enddo
	  call total(PAR,QTnext,RTnext) !,NK,Q,R,QT,RT)		!QT, RT - output
	 ! N=Nx+Nalf                         !Nx=NK-1, Nalf=2*(NG-1)
	 	PAR(NK)=1.0d0
	  DO  I=1,Nx                                                       
      PAR(I)=PARx(I)
	  PAR(NK)=PAR(NK)-PAR(I)
	  enddo
	  !  DO D=1,Nx
	  ! if (dabs(FU(D)).gt.1.0d-12) HH=1
	  ! enddo
!	    if ((L.le.18).and.(L.ge.17).and.(iter.lt.200)) HH=0	  !((L.le.1).and.(iter.lt.100)) 
!	   if (L.eq.Lmax) HH=0
    call total(PAR,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
	RTx=RT
do i=1,Lmax
geom(i)=volume(i)
!LL=i
!if (shape.eq.'sph') then
! geom(i)=!(3.0d0*LL**2-3.0d0*LL+1.0d0) !/(1.0d0*LLmax**3)
! elseif (shape.eq.'cyl') then
! geom(i)=(2.0d0*LL-1.0d0) !/LLmax**2
 !else
 !geom(i)=1.0d0 !/LLmax
! endif
enddo 
print *,geom(Lmax)
sumL=0.0d0
	do I=1,Nx
	Philoc(I,L)=R(I)*Xloc(I,L)/RT
	!Phi0=Philoc(I,L)
	sumL=sumL+Philoc(I,L) !Xloc(I,L)
	enddo
sumLnext=0.0d0
	do I=1,Nx
	Philoc(I,L+1)=R(I)*Xloc(I,L+1)/RTnext
	sumLnext=sumLnext+Philoc(I,L+1) !Xloc(I,L+1)
	enddo
	if (L.eq.(Lmax)) then
sumLprev=0.0d0
	do I=1,Nx
	Philoc(I,L-1)=R(I)*Xloc(I,L-1)/RTpr
	sumLprev=sumLprev+Philoc(I,L-1) !Xloc(I,L-1)
	enddo
	endif
do D=1,Nx
if (L.lt.(Lmax)) then
PARold(D)=Philoc(D,L+1)
DPmin(D)=max(-Philoc(D,L),(sumLnext-1.0d0)*geom(L+1)/geom(L))
DPmax(D)=min((1.0d0-sumL),Philoc(D,L+1)*geom(L+1)/geom(L))
elseif (L.eq.(Lmax)) then
PARold(D)=Philoc(D,L-1)
DPmin(D)=max(-Philoc(D,L),(sumLprev-1.0d0)*geom(L-1)/geom(L))
DPmax(D)=min((1.0d0-sumL),Philoc(D,L-1)*geom(L-1)/geom(L))
endif

if (DPmin(D).ge.DPmax(D)) then
print *,'DPmin>DPmax in canonic newtX',D,DPmin(D),DPmax(D),L
stop
endif
enddo
	   call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !0= correction of initial PAR
	   	 do while (mainK.ne.0)
333	 kkk=kkk+1
 
	 DO J=1,Nx
	! if (id(J).ne.0) then
	!	if ((PARR(J)*id(J)).gt.0.0d0) then 
	!	PARR(J)=PARR(J)*1.1d0
	!	 else
	!!	PARR(J)=PARR(J)*0.8d0	
	!	endif
	 PAR(J)=PARx(J)+(PAR(J)-PARx(J))*0.9d0 !X(J)+PARR(J)
	!!		if (kkk.gt.100000) then
	!!		print *,'many steps after iter in NEWTx',id,PAR
	!!		stop
	!!		endif
	! endif
	enddo	!J
	 call check_borderX(PAR,L,Lmax,Nx,id,mainK,0)           !-1 = checking   PAR
	 PAR(NK)=1.0d0
	 	 DO J=1,Nx
	 	PAR(NK)=PAR(NK)-PAR(J)
	enddo	!J
	  call total(PAR,QTcur,RTcur) !,NK,Q,R,QT,RT)		!QT, RT - output

	maink=0
	 summa=0.0d0
      DO K=1,Nx
	  	summa=summa+PAR(K)*R(k)/RTcur-PARx(k)*R(k)/RT
      enddo
	if (L.lt.(Lmax)) then
	if (summa.le.max((sumLnext-1.0d0)*geom(L+1)/geom(L),-sumL)) 	mainK=1
	if (summa.ge.min(sumLnext*geom(L+1)/geom(L),(1.0d0-sumL)))	mainK=1
elseif (L.eq.(Lmax)) then
if (summa.le.max((sumLprev-1.0d0)*geom(L-1)/geom(L),-sumL)) 	mainK=1
if (summa.ge.min(sumLprev*geom(L-1)/geom(L),(1.0d0-sumL)))	mainK=1
endif

	 DO J=1,Nx
	if (((PAR(J)*R(J)/RTcur-PARx(J)*R(J)/RT).le.DPmin(J)).or.((PAR(J)*R(J)/RTcur-PARx(J)*R(J)/RT).ge.DPmax(J))	 ) mainK=1
	enddo
	enddo		 !while
   	 PAR(NK)=1.0d0
	 	 DO J=1,Nx
	 	PAR(NK)=PAR(NK)-PAR(J)
	enddo	!J
	  call total(PAR,QTcur,RTcur) !,NK,Q,R,QT,RT)		!QT, RT - output

!!!!change neighbor layer

if (L.lt.(Lmax)) then
RTnext=(1.0d0-sumLnext)/R(NK)
DO I=1,Nx
Philoc(I,L+1)=PARold(I)-(PAR(I)*R(I)/RTcur-PARx(I)*R(I)/RT)*geom(L)/geom(L+1)
RTnext=RTnext+Philoc(I,L+1)/R(I)-(Philoc(I,L+1)-PARold(I))/R(NK)
enddo
RTnext=1.0d0/RTnext
	do k=1,Nx
	PARnext(k)=Philoc(k,L+1)/R(k)*RTnext !Xloc(k,L+1)
	enddo
  call check_borderX(PARnext,L+1,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	 if (maink.ne.0) then
	 goto 333
	 endif
!	call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
 !  	 write(*,*)'!!!!!!!!!!!!!! first NEWTalf ended' 
 if (L.ne.Lmaxx) then
call NEWTalf_chem(Nx,Nalf,PARnext,PARalfnext,DPalf,FU,DDalf,L+1,EPS,EPF,KEY,ghom)
   	DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
     PARnext(K)=PARalfnext(K)
	 enddo
CALL inhomsol(Nx,Nalf,PARnext,FU,L+1,Wprev,ghom)
 endif
elseif (L.eq.(Lmax)) then
RTpr=(1.0d0-sumLprev)/R(NK)
DO I=1,Nx
  Philoc(I,L-1)=PARold(I)-(PAR(I)*R(I)/RTcur-PARx(I)*R(I)/RT)*geom(L)/geom(L-1)
  RTpr=RTpr+Philoc(I,L-1)/R(I)-(Philoc(I,L-1)-PARold(I))/R(NK)
  enddo
  RTpr=1.0d0/RTpr
	do k=1,Nx
	PARprev(k)=Philoc(k,L-1)/R(k)*RTpr !Xloc(k,L-1)
	enddo
	 call check_borderX(PARprev,L-1,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	 if (maink.ne.0) then
	 goto 333
	 endif
 ! 	call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
 !  	 write(*,*)'!!!!!!!!!!!!!! first NEWTalf ended' 
  call NEWTalf_chem(Nx,Nalf,PARprev,PARalfprev,DPalf,FU,DDalf,L-1,EPS,EPF,KEY,ghom)
  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
     PARprev(K)=PARalfprev(K)
	 enddo
CALL inhomsol(Nx,Nalf,PARprev,FU,L-1,Wprev,ghom)
endif
!enddo
!!!! 
	
	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	 write(*,*)'!!!!!!!!!!!!!! first NEWTalf ended'  
!	 stop 
	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(I)=PARalf(I)
	  enddo
	   CALL inhomsol(Nx,Nalf,PAR,FU,L,Wprev,ghom)
	   PARx=PAR
	   RTx=RTcur
	   RT=RTcur
     if (HH /= 0) then                                              
                                              
  !    NT2=Nx+1 !Nalf+1                                                           
   5  DO 1 I=1,Nx                                                        
      F(I)=FU(I) !+FU1(I)                                                        
   !   DD(I,NT2)=-F(I)*LAM                                              
   1  X(I)=PAR(I)                                                      
      KK=KK+1                                                           
      DO I=1,Nx
	  DPx(I)=DDPx(I)
	   !if (I.le.Nx) then 
	  !  if  ((PAR(I).eq.1.0d0).and.(DPx(I).gt.0.0d0)) DPx(I)=-DPx(I)
	  !   if  ((PAR(I).le.1.0d-150).and.(DPx(I).lt.0.0d0)) DPx(I)=-DPx(I)
	  Philoc(I,L)=PAR(I)*R(I)/RT+DPx(I)
	  RTcur=(1.0d0-sumL-DPx(I))/R(NK)
      DO k=1,Nx
	RTcur=RTcur+Philoc(k,L)/R(k)
	enddo
	RTcur=1.0d0/RTcur
	DO k=1,Nx
	PAR(k)=	Philoc(k,L)*RTcur/R(k)
	enddo                                              
      !!! PAR(I)=PAR(I)+DPx(I)							!mole fract.
   !    call minmaxX(PAR,L,Lmax,Nx,allK,KEYx,KEYxL,KEYxend,KEYid,-I)
    call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
!	if (mainK.eq.0) call check_ALFA_step(PAR,L,Nx,Nalf,X,DPx,I,KEY_step)
!	if (mainK.eq.0) KEYid(I)=KEY_step
	 ! if (I.eq.2) print *,KEYx(I),KEYxL(I),KEYxend(I),KEYid(I)
	  kkk=0
	 do while ((mainK.ne.0).or.(DPx(I).le.DPmin(I)).or.(DPx(I).ge.DPmax(I))	 )
	111 kkk=kkk+1
			if (kkk.gt.1000) then
			print *,'many steps of DP in NEWTx',DPx(I),I,'keyid=',KEYid(I)
				Philoc(I,L)=Philoc(I,L)-DPx(I)
			  	DO k=1,Nx
				PAR(k)=	Philoc(k,L)*RT/R(k)
				enddo
			call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
			print *,'many'
			if (L.lt.(Lmax)) then
		RTnext=(1.0d0-sumLnext)/R(NK)
		DO k=1,Nx
		Philoc(k,L+1)=PARold(k)-(PAR(k)*R(k)/RT-PARx(k)*R(k)/RTx)*geom(L)/geom(L+1)
		RTnext=RTnext+Philoc(k,L+1)/R(k)-(Philoc(k,L+1)-PARold(k))/R(NK)
				enddo
		RTnext=1.0d0/RTnext
	do k=1,Nx
	PARnext(k)=Philoc(k,L+1)/R(k)*RTnext !Xloc(k,L+1)
	enddo
   	call NEWTalf_chem(Nx,Nalf,PARnext,PARalf,DPalf,FU,DDalf,L+1,EPS,EPF,KEY,ghom)
endif
			goto 99
		!	stop
			endif
		Philoc(I,L)=Philoc(I,L)-DPx(I)
	!!! PAR(I)=PAR(I)-DPx(I) !*X(I)
		if (dabs(DPx(I)).lt.1.0d-16) then
		DPx(I)=-DPx(I)*1.0d10
		endif 

	if  ((KEYid(I).eq.1).or.((PAR(I)*R(I)/RTcur-PARx(I)*R(I)/RTx).le.DPmin(I)) ) then
		if (DPx(I).gt.0.0d0) DPx(I)=DPx(I)*1.1d0
		if (DPx(I).lt.0.0d0) DPx(I)=DPx(I)*0.9d0
	else if  ((KEYid(I).eq.(-1)).or.((PAR(I)*R(I)/RTcur-PARx(I)*R(I)/RTx).ge.DPmax(I))	 ) then
		if (DPx(I).gt.0.0d0) DPx(I)=DPx(I)*0.9d0
		if (DPx(I).lt.0.0d0) DPx(I)=DPx(I)*1.1d0
		else if (mainK.ne.0) then
		DPx(I)=DPx(I)*0.9d0
	endif
		DO k=1,Nx
	PAR(k)=	Philoc(k,L)*RT/R(k)
	enddo
  Philoc(I,L)=PAR(I)*R(I)/RT+DPx(I)
	  RTcur=(1.0d0-sumL-DPx(I))/R(NK)
      DO k=1,Nx
	RTcur=RTcur+Philoc(k,L)/R(k)
	enddo
	RTcur=1.0d0/RTcur
	DO k=1,Nx
	PAR(k)=	Philoc(k,L)*RTcur/R(k)
	enddo
	
!!!	PAR(I)=PAR(I)+DPx(I) !*X(I)
!	call minmaxX(PAR,L,Lmax,Nx,allK,KEYx,KEYxL,KEYxend,KEYid,-I)
	 call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	!!!!!!!
	if ((kkk.lt.150).and.(L.ne.Lmax)) then
!	if (KEYid(I).eq.0)	call check_ALFA_step(PAR,L,Nx,Nalf,X,DPx,I,KEY_step)
!	if (KEYid(I).eq.0) KEYid(I)=KEY_step
	endif
!		  if (I.eq.2) print *,KEYx(I),KEYxL(I),KEYxend(I),KEYid(I)
	enddo	!while
!	 print *,'after DPx',I,PAR(1:Nx),L,PARx(1:Nx)
!	  print *,I !,PAR(5:6),AL_up(2,L),AL_al(2,L),chem_al(1,L),chem_al(2,L)
call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
  
!	if ((kkk.ge.150).or.(L.eq.Lmax))	 	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	  	  print *,I,PAR(5:6),AL_up(2,L),AL_al(2,L),chem_al(1,L),chem_al(2,L)
!	if ((kkk.lt.150).and.(L.ne.Lmax))  call ALFA_step (Nx,Nalf,X,PARalf,DPx,I,L)
    DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo

!!!	  call calcALloc(PAR,L,ALlocnew,sumXX)
!!!	  	call check_chemdwn(PAR,L,KKEY)

	   summa=0.0d0
      DO D=1,Nx
	  	summa=summa+PAR(D)
      enddo
call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	 
!	 write(*,*)'!!!!!!!!!!!!!! NEWTalf after dx ended'   
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo
                                                       
      CALL inhomsol(Nx,Nalf,PAR,FU,L,Wprev,ghom)
                                            
       DO  J=1,Nx                                                       
	 DD(J,I)=(FU(J)-F(J))/DPx(I) !(dlog(PAR(I)/X(I)) )  !DPx(I) 		! DD(J,I) changed!!!
	   enddo                  !cycle J                     
                                            
    !!!  PAR(I)=PAR(I)-DPx(I)
		Philoc(I,L)=Philoc(I,L)-DPx(I)
		DO  J=1,Nx  
		PAR(J)=	Philoc(J,L)*RT/R(J)
		enddo
	  enddo                    !cycle I 

	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	 write(*,*)'!!!!!!!!!!!!!!'   
	  DO  K=Nx+1,Nx+Nalf+Nbonds                                                      
      PAR(K)=PARalf(K)
	  enddo
                                                       
      CALL inhomsol(Nx,Nalf,PAR,FU,L,Wprev,ghom) 
 
	  do I=1,Nx
	  PARR(I)=0.0d0 !PAR(I)
	  enddo
	   DO  J=1,Nx                   
	   FF(J)=-F(J)
	   enddo                  !cycle J 
	!   print *,'x=',X
	!   print *,'DD=',DD   
	  CALL DLSLRG(Nx, DD, 150, FF, 1, PARR)                                          
!      CALL GA(DD,Nx,NT2,PAR) 	!PAR is deltaX after GA
	 ! write(*,*)'!!!after dlslrg!!X(1)before=',X(1),X(2),X(3),'X(1)after=',X(1)+PARR(1),X(2)+PARR(2),X(3)+PARR(3)                                             
    
	 DO I=1,Nx 			   
		!!! PAR(I)=X(I)+PARR(I)
		Philoc(I,L)=PAR(I)*R(I)/RT+PARR(I)
	 enddo						
	  RTcur=(1.0d0-sumL)/R(NK)
DO I=1,Nx
RTcur=RTcur+Philoc(I,L)/R(I)-PARR(I)/R(NK)
enddo
RTcur=1.0d0/RTcur
DO I=1,Nx
PAR(I)=Philoc(I,L)*RTcur/R(I)
enddo
 	   call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
!	 if (mainK.eq.0)	call check_ALFA_step(PAR,L,Nx,Nalf,X,PARR,0,KEY_step)
!	if (mainK.eq.0) mainK=KEY_step
	  kkk=0
	 do while (mainK.ne.0)
222	 kkk=kkk+1
		 DO I=1,Nx 			   
		Philoc(I,L)=Philoc(I,L)-PARR(I)
	 enddo
	 DO J=1,Nx
	! if (id(J).ne.0) then
	!	if ((PARR(J)*id(J)).gt.0.0d0) then 
	!	PARR(J)=PARR(J)*1.1d0
	!	 else
		PARR(J)=PARR(J)*0.95d0	
	!	endif
!!!	 PAR(J)=X(J)+PARR(J)
		
	! endif

	enddo	!J
		if (kkk.gt.1000) then
			print *,'many steps after iter in NEWTx',id(1:Nx),PAR(1:Nx)
		!	PAR(J)=X(J)
			PARR=0.0d0
			PAR=PARx
			goto 100
			!stop
			endif
						
!	  RTcur=(1.0d0-sumL)/R(NK)
!DO I=1,Nx
!RTcur=RTcur+Philoc(I,L)/R(I)-PARR(I)/R(NK)
!enddo
!RTcur=1.0d0/RTcur
DO I=1,Nx
PAR(I)=Philoc(I,L)*RT/R(I)
enddo

		 DO I=1,Nx 			   
		!!! PAR(I)=X(I)+PARR(I)
		Philoc(I,L)=PAR(I)*R(I)/RT+PARR(I)
	 enddo						
	  RTcur=(1.0d0-sumL)/R(NK)
DO I=1,Nx
RTcur=RTcur+Philoc(I,L)/R(I)-PARR(I)/R(NK)
enddo
RTcur=1.0d0/RTcur
DO I=1,Nx
PAR(I)=Philoc(I,L)*RTcur/R(I)
enddo
	 call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)           !-1 = checking   PAR
	 if ((kkk.lt.150).and.(L.ne.Lmax))  then
!	if (mainK.eq.0)	call check_ALFA_step(PAR,L,Nx,Nalf,X,PARR,0,KEY_step)
!	if (mainK.eq.0) mainK=KEY_step
	endif
	 summa=0.0d0
      DO K=1,Nx
	  	summa=summa+Philoc(k,L)-PARx(k)*R(k)/RTx !PARR(K)
      enddo
	if (L.lt.(Lmax)) then
	if (summa.le.max((sumLnext-1.0d0)*geom(L+1)/geom(L),-sumL)) 	mainK=1
	if (summa.ge.min(sumLnext*geom(L+1)/geom(L),(1.0d0-sumL)))	mainK=1
elseif (L.eq.(Lmax)) then
if (summa.le.max((sumLprev-1.0d0)*geom(L-1)/geom(L),-sumL)) 	mainK=1
if (summa.ge.min(sumLprev*geom(L-1)/geom(L),(1.0d0-sumL)))	mainK=1
endif

	 DO J=1,Nx
	if (((Philoc(J,L)-PARx(J)*R(J)/RTx).le.DPmin(J)).or.((Philoc(J,L)-PARx(J)*R(J)/RTx).ge.DPmax(J))	 ) mainK=1
	enddo
	enddo		 !while
 PRINT *,'.'
 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)

 print *,'L'
  !	if ((kkk.ge.150).or.(L.eq.Lmax))	 	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	  	  print *,I,PAR(5:6),AL_up(2,L),AL_al(2,L),chem_al(1,L),chem_al(2,L)
!	if ((kkk.lt.150).and.(L.ne.Lmax))  call ALFA_step (Nx,Nalf,X,PARalf,DPx,I,L)
    DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo
 100	  A=0.0d0   
 DO I=1,Nx 
	  A=A+DABS(PARR(I))
    !  PAR(I)=PAR(I)+PARR(I)
	!  if (KK.gt.180) write(*,*)'x=',DABS(PAR(I))
	  enddo
	!!!!change neighbor layer
!DO I=1,Nx
if (L.lt.(Lmax)) then
RTnext=(1.0d0-sumLnext)/R(NK)
DO I=1,Nx
Philoc(I,L+1)=PARold(I)-(PAR(I)*R(I)/RTcur-PARx(I)*R(I)/RTx)*geom(L)/geom(L+1)
RTnext=RTnext+Philoc(I,L+1)/R(I)-(Philoc(I,L+1)-PARold(I))/R(NK)
enddo
RTnext=1.0d0/RTnext
	do k=1,Nx
	PARnext(k)=Philoc(k,L+1)/R(k)*RTnext !Xloc(k,L+1)
	enddo
!!!DO I=1,Nx
!!!Xloc(I,L+1)=PARold(I)-(PAR(I)-PARx(I))*geom(L)/geom(L+1)
!!!enddo
!!!	do k=1,Nx
!!!	PARnext(k)=Xloc(k,L+1)
!!!	enddo
  call check_borderX(PARnext,L+1,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
  PRINT *,'!'
	 if (maink.ne.0) then
	 RTnext=(1.0d0-sumLnext)/R(NK)
DO I=1,Nx
Philoc(I,L+1)=PARold(I)
RTnext=RTnext+Philoc(I,L+1)/R(I)
enddo
RTnext=1.0d0/RTnext
	 do k=1,Nx
	PARnext(k)=PARold(k)/R(k)*RTnext !Xloc(k,L+1)
	enddo
	call check_borderX(PARnext,L+1,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR 	
	 goto 222
	 endif
  if (L.ne.Lmaxx) then
  PRINT *,'!!'
call NEWTalf_chem(Nx,Nalf,PARnext,PARalfnext,DPalf,FU,DDalf,L+1,EPS,EPF,KEY,ghom)
print *,'L+1'
   	DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
     PARnext(K)=PARalfnext(K)
	 enddo
CALL inhomsol(Nx,Nalf,PARnext,FU,L+1,Wprev,ghom)
  endif
elseif (L.eq.(Lmax)) then
RTpr=(1.0d0-sumLprev)/R(NK)
DO I=1,Nx
  Philoc(I,L-1)=PARold(I)-(PAR(I)*R(I)/RTcur-PARx(I)*R(I)/RTx)*geom(L)/geom(L-1)
  RTpr=RTpr+Philoc(I,L-1)/R(I)-(Philoc(I,L-1)-PARold(I))/R(NK)
  enddo
  RTpr=1.0d0/RTpr
	do k=1,Nx
	PARprev(k)=Philoc(k,L-1)/R(k)*RTpr !Xloc(k,L-1)
	enddo
!!!DO I=1,Nx
 !!! Xloc(I,L-1)=PARold(I)-(PAR(I)-PARx(I))*geom(L)/geom(L-1)
!!!  enddo
!!!	do k=1,Nx
!!!	PARprev(k)=Xloc(k,L-1)
!!!	enddo
	 call check_borderX(PARprev,L-1,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	 if (maink.ne.0) then
	 RTpr=(1.0d0-sumLprev)/R(NK)
DO I=1,Nx
  !Philoc(I,L-1)=PARold(I)-(PAR(I)*R(I)/RTcur-PARx(I)*R(I)/RTx)*geom(L)/geom(L-1)
  RTpr=RTpr+PARold(I)/R(I) !-(Philoc(I,L-1)-PARold(I))/R(NK)
  enddo
  RTpr=1.0d0/RTpr
	do k=1,Nx
	PARprev(k)=PARold(K)/R(k)*RTpr !Xloc(k,L-1)
	enddo
	call check_borderX(PARprev,L-1,Lmax,Nx,id,mainK,-1)	   !-1 = checking PAR
	 goto 222
	 endif

  call NEWTalf_chem(Nx,Nalf,PARprev,PARalfprev,DPalf,FU,DDalf,L-1,EPS,EPF,KEY,ghom)
  DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
     PARprev(K)=PARalfprev(K)
	 enddo
CALL inhomsol(Nx,Nalf,PARprev,FU,L-1,Wprev,ghom)
endif
!enddo
!!!!
	 !!!call check_borderX(PAR,L,Lmax,Nx,id,mainK,-1)           !-1 = checking   PAR
 	 !!! if (maink.ne.0) then
	!!! goto 222
	!!! endif
	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
	 !	if ((kkk.ge.150).or.(L.eq.Lmax))	 	 call NEWTalf_chem(Nx,Nalf,PAR,PARalf,DPalf,FU,DDalf,L,EPS,EPF,KEY,ghom)
!	  	  print *,I,PAR(5:6),AL_up(2,L),AL_al(2,L),chem_al(1,L),chem_al(2,L)
!	if ((kkk.lt.150).and.(L.ne.Lmax))  call ALFA_step (Nx,Nalf,X,PARalf,DPx,I,L)
    DO  K=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(K)=PARalf(K)
	  enddo

	  DO  I=Nx+1,Nx+Nalf+Nbonds                                                       
      PAR(I)=PARalf(I)
	  enddo         
!	   write(*,*)'x=',PAR(1)
	  CALL inhomsol(Nx,Nalf,PAR,FU,L,Wprev,ghom) 

	  write(*,*)'!!!Layer = ',L

!	 IF(KK.GT.100) EPS=1.0D-5 
 !IF(KK.GT.100) EPF=5.0D-1 !GOTO 7  	                                            
      A=A/Nx                                                             
      IF ((A.GT.EPS).and.(kk.lt.10)) GOTO 5
                                            
      A=0.0d0                                                              
      DO 6 I=1,Nx 
	  IF(KK.GT.180) write(*,*)'200',DABS(FU(I))
!	  if (L.eq.2) A=A+DABS(FU1(I))                                                       
   6  A=A+DABS(FU(I))
   print *,A
   !stop
  ! if (dabs(y_upch(1,L)*teta_chem(1,L)/m_up(L)-y_dwnch(1,L+1)*teta_chem(1,L+1)/m_dwn(L+1)))  A=A+1000.0
   write(*,*)'A=sum of FU =',A
   !if(A.lt.1000) stop  
   !   A=A+dabs(y_upch(1,L)*teta_chem(1,L)/m_up(L)-y_dwnch(1,L+1)*teta_chem(1,L+1)/m_dwn(L+1))                                              
      KEY=1                                                             
      IF(A.GT.EPF)  KEY=-1
	  !if (L.eq.13) print *,'epf=',EPF
	   
	 ! IF ((A.GT.EPS).and.(kk.lt.10)) GOTO 5
	  !  r=0
      !  DO  I=1,Nx
	 !   if (PAR(I).eq.1.0d0) r=r+1
     !   enddo
	  !  if (r.gt.1) GOTO 5
  endif          !HH  	
99	  DO  I=1,Nx
	!  if (I.le.Nx) then                                                        
      PARx(I)=PAR(I)
	!  else
	!  PARalf(I)=PAR(I)					!output values as input for the next iteration in micelle_energy
	!  endif 
	  enddo                                  
   7  CONTINUE                    
                           
      RETURN                                                            
      END                                                             
!********************************************************************************************************************


!*********************************************************************************************************************
!<NIKITA>
	  ! ==================================================================
!  REPLACEMENTS FOR IMSL ROUTINES
! ==================================================================

SUBROUTINE dlslrg(N, A, LDA, B, IPATH, X)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDA, IPATH
    REAL*8, INTENT(INOUT) :: A(LDA, N)
    REAL*8, INTENT(IN) :: B(N)
    REAL*8, INTENT(OUT) :: X(N)
    
    REAL*8 :: LU(LDA, N), R(N)
    INTEGER :: IPIV(N), INFO, I, J
    
    ! Copy A and B because Gaussian elimination modifies them
    LU = A
    R = B
    
    ! Solve using basic Gaussian Elimination with Partial Pivoting
    CALL local_dgesv(N, LU, LDA, IPIV, R, INFO)
    
    if (INFO /= 0) then
        print *, "Error: dlslrg matrix is singular or solver failed."
    endif
    X = R
END SUBROUTINE dlslrg

SUBROUTINE dzreal(F, EPS, ETA, NSIG, XGUESS, N, MAXIT, X, XFN, INFO)
    IMPLICIT NONE
    EXTERNAL :: F
    REAL*8 :: F
    REAL*8, INTENT(IN) :: EPS, ETA, NSIG, XGUESS
    INTEGER, INTENT(IN) :: N, MAXIT
    REAL*8, INTENT(INOUT) :: X(N)    ! Array matching X1
    REAL*8, INTENT(OUT) :: XFN(N)  ! Array matching Xmin
    INTEGER, INTENT(OUT) :: INFO(N) ! Array matching INFO (fixes rank mismatch)
    
    REAL*8 :: a, b, fa, fb, c, fc
    INTEGER :: iter, i
    
    ! Process each root tracking slot requested by N
    do i = 1, N
        a = X(i) - 1.0d-1
        b = X(i) + 1.0d-1
        
        fa = F(a)
        fb = F(b)
        
        if (fa * fb > 0.0d0) then
            a = X(i) - 1.0d1
            b = X(i) + 1.0d1
            fa = F(a)
            fb = F(b)
        endif
        
        INFO(i) = 0
        do iter = 1, MAXIT
            c = (a + b) / 2.0d0
            fc = F(c)
            
            if (dabs(fc) < EPS .or. (b - a) / 2.0d0 < ETA) then
                X(i) = c
                XFN(i) = fc
                exit
            endif
            
            if (fa * fc < 0.0d0) then
                b = c
                fb = fc
            else
                a = c
                fa = fc
            endif
            
            if (iter == MAXIT) INFO(i) = 1
        enddo
        X(i) = c
        XFN(i) = fc
    enddo
END SUBROUTINE dzreal


! Subroutine that emulates standard LAPACK DGESV behavior cleanly
SUBROUTINE local_dgesv(N, A, LDA, IPIV, B, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDA
    REAL*8, INTENT(INOUT) :: A(LDA, N)
    REAL*8, INTENT(INOUT) :: B(N)
    INTEGER, INTENT(OUT) :: IPIV(N)
    INTEGER, INTENT(OUT) :: INFO
    INTEGER :: I, J, K, P
    REAL*8 :: MAXVAL, TEMP, FACTOR
    
    INFO = 0
    DO I = 1, N
        IPIV(I) = I
    END DO
    
    DO K = 1, N - 1
        MAXVAL = DABS(A(K, K))
        P = K
        DO I = K + 1, N
            IF (DABS(A(I, K)) > MAXVAL) THEN
                MAXVAL = DABS(A(I, K))
                P = I
            END IF
        END DO
        
        IF (MAXVAL < 1.0D-20) THEN
            INFO = K
            RETURN
        END IF
        
        IF (P /= K) THEN
            IPIV(K) = P
            DO J = K, N
                TEMP = A(K, J)
                A(K, J) = A(P, J)
                A(P, J) = TEMP
            END DO
            TEMP = B(K)
            B(K) = B(P)
            B(P) = TEMP
        END IF
        
        DO I = K + 1, N
            FACTOR = A(I, K) / A(K, K)
            A(I, K) = FACTOR
            DO J = K + 1, N
                A(I, J) = A(I, J) - FACTOR * A(K, J)
            END DO
            B(I) = B(I) - FACTOR * B(K)
        END DO
    END DO
    
    IF (DABS(A(N, N)) < 1.0D-20) THEN
        INFO = N
        RETURN
    END IF
    
    ! Back substitution
    B(N) = B(N) / A(N, N)
    DO I = N - 1, 1, -1
        TEMP = B(I)
        DO J = I + 1, N
            TEMP = TEMP - A(I, J) * B(J)
        END DO
        B(I) = TEMP / A(I, I)
    END DO
END SUBROUTINE local_dgesv

!</NIKITA>
