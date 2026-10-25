

Subroutine mol_prop(NKK,type_mol,Xaq,Xorg)
      implicit none
      REAL*8 QT,RT,T,TMo,TT,TTMo,R(150),RMo(150),RS(150,10),RSMo(150,10),QSMo(150,10),Q(150),QMo(150),QS(150,10),vol_ref, q_ref
	  integer j, NKK,NK1
	  real*8 LIS(150,10),AL(10),ALI(10),ET(10,10),GP(150),RSTMo(150,10),X(150),X1(150),X2(150),XX(10),Y(150),GAM1(150),GAM2(150)
	  real*8 ZZ,ZR,ZZMo,ZRMo,JX,IJ,Xcomp(150),Xseg(150),Xaq(150),Xorg(150),RE
	  integer b,M,NK,NG,Nmol,I,type_mol(150),Nbonds
	  integer NKMo,NGMo,JXMo,IJMo
	  COMMON /C/ T,TT,R,RS,Q,QS,vol_ref,q_ref
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                      
	  REAL*8 BLIS(150,10),BAL(10),BALI(10),BET(10,10),BGP(150)
	  COMMON /H/ BLIS,BAL,BALI,BET,BGP  			   ! Dummy, apart from  BET=ET
      COMMON /D89/ ZZMo,ZRMo,NKMo,NGMo,JXMo,IJMo                                    
	  COMMON /C89/ QT,RT,TMo,TTMo,RMo,RSMo,RSTMo,QMo,QSMo
	  COMMON /H89/ LIS,AL,ALI,ET,GP
      COMMON /A89/ X,X1,X2,XX,Y,GAM1,GAM2

    open(22,file='slosh.dat')

! swap between common blocks      
	  ZZmo=ZZ		 
	  ZRMo=ZR
	  NKMo=NKK
	  NGMo=NG
	  TMo=T
	  TTMo=TT
! components geometry parameters
 b=0
	  DO I=1,NK
	 if (type_mol(I).eq.0) then
	  b=b+1
	   do  J=1,NG
	   QSMo(b+Nmol,J)=QS(I,J)
	   RSMo(b+Nmol,J)=RS(I,J)
	   enddo ! groups
	 else
	   do M=1,Nmol
	   if (type_mol(I).eq.M) then
	   do J=1,NG
	   QSMo(M,J)=QSMo(M,J)+QS(I,J)
	   RSMo(M,J)=RSMo(M,J)+RS(I,J)
	   enddo  ! groups
	   endif
	   enddo ! chains
	  endif
	   if (Nmol.eq.0) then
	   do J=1,NG
       QSMo(I,J)=QS(I,J)
	   RSMo(I,J)=RS(I,J)
	   enddo ! groups
	   endif 
	  enddo	! segments

	  do i=1,NKK
      RMo(i)=0.d0
	  QMo(i)=0.d0
	  do j=1,NG
	  RMo(i)=RMo(i)+RSMo(i,j)
  	  QMo(i)=QMo(i)+QSMo(i,j) 
      enddo 
	  enddo
! group ebergetica swapping and calculation of component's GP
	        CALL PARAM (T,TT)
			ET=BET                                                       

! change component numeration to normal:  holes are the last component and last group
      NK1=NKK-1                                                          
      do j=1,NG
	  RSMo(NK1,j)=RSMo(NKK,j)
 	  QSMo(NK1,j)=QSMo(NKK,j)
  	  RSMo(NKK,j)=0.0
  	  QSMo(NKK,j)=0.0
	  enddo
	  RSMo(NKK,NG)=1.0
	  QSMo(NKK,NG)=1.0
	  
	  RMo(NK1)=RMo(NKK)
  	  QMo(NK1)=QMo(NKK)
	  RMo(NKK)=1.0d0
	  QMo(NKK)=1.0d0
   	  
	  GP(NKK)=0.0

      do I=1,NK1                                                      
      GP(I)=0.0                                                         
      do J=1,NG                                                       
      AL(J)=QSMo(I,J)/QMo(I)                                                
	  enddo                                                          
      CALL IKSMo
      do J=1,NG                                                      
      GP(I)=GP(I)-ZZ*QSMo(I,J)*DLOG(XX(J))
	  enddo                                
      enddo ! components

	  X(1)=0.15d0
	  X(2)=0.25d0
	  X(3)=1.0D0-X(1)-X(2)
	  call LVE_flash

! change component numeration back: water is the last component
      RE=X1(NKK)
	  X1(NKK)=X1(NK1)
	  X1(NK1)=RE				 
      RE=X2(NKK)
	  X2(NKK)=X2(NK1)
	  X2(NK1)=RE
! determine the aqueous and organic phases
     if (X1(NKK).GT.X2(NKK))then
	 Xaq=X1
	 Xorg=X2
	 else				 
	 Xaq=X2
	 Xorg=X1
	 endif
 return 
 end
!**************************************************************************************************
      SUBROUTINE ACTCF(DX,GAM)                                          
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER I,J,L                                                     
      REAL*8 BB,DX(150),GAM(150),A                                        
      INTEGER NK,NG,JX,IJ                                               
      REAL*8 X1(150),X2(150),XX(10),Y(150),GAM1(150),GAM2(150),X(150),ET(10,10),GP(150),LIS(150,10),AL(10),ALI(10)
	  real*8 QT,RT,T,TT,ZZ,ZR,RMo(150),RSMo(150,10),RSTMo(150,10),QMo(150),QSMo(150,10)
      COMMON /A89/ X,X1,X2,XX,Y,GAM1,GAM2                               
	  COMMON /C89/ QT,RT,T,TT,RMo,RSMo,RSTMo,QMo,QSMo
      COMMON /D89/ ZZ,ZR,NK,NG,JX,IJ
      COMMON /H89/ LIS,AL,ALI,ET,GP                                     
      RT=0.0                                                            
      QT=0.0                                                            
      DO 1 I=1,NK                                                       
      RT=RT+RMo(I)*DX(I)                                                  
   1  QT=QT+QMo(I)*DX(I)                                                  
      DO 2 I=1,NK                                                       
      A=RMo(I)/RT                                                         
      BB=A*QT/QMo(I)                                                      
   2  GAM(I)=1-A+DLOG(A)-ZR*QMo(I)*(1-BB+DLOG(BB))                        
      DO 3 L=1,NG                                                       
      A=0.0                                                             
      DO 4 J=1,NK                                                       
   4  A=A+QSMo(J,L)*DX(J)                                                 
   3  AL(L)=A/QT                                                        
      CALL IKSMo                                                          
      DO 5 I=1,NK                                                       
      DO 6 L=1,NG                                                       
   6  GAM(I)=GAM(I)+ZZ*QSMo(I,L)*DLOG(XX(L))                              
   5  GAM(I)=DEXP(GAM(I)+GP(I))                                         
      RETURN                                                            
      END                                                                
!**************************************************************************************************
      SUBROUTINE IKSMo                                                    
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER I,K                                                       
      REAL*8 BB,RXX(10),RE                                        
      INTEGER NK,NG,JX,IJ                                               
	  REAL*8 ZR,ZZ,X(150),X1(150),X2(150),XX(10),Y(150),GAM1(150),GAM2(150),LIS(150,10),AL(10),ALI(10),ET(10,10),GP(150)
      COMMON /A89/ X,X1,X2,XX,Y,GAM1,GAM2                               
      COMMON /D89/ ZZ,ZR,NK,NG,JX,IJ                                    
      COMMON /H89/ LIS,AL,ALI,ET,GP                                     
      RE=1.0                                                            
      J1=0                                                              
       DO 1 K=1,NG                                                       
      XX(K)=RE                                                          
   1  RE=RE-1.0D-04                                                     
  19  J1=J1+1                                                           
      DO 2 K=1,NG                                                       
   2  XX(K)=DABS(XX(K))                                                 
      DO 3 K=1,NG                                                       
      BB=0.0                                                            
      DO 4 I=1,NG                                                       
    4 BB=BB+AL(I)*XX(I)*ET(K,I)                                         
   3  RXX(K)=1.0/BB                                                     
      DO 6 K=1,NG                                                       
   6  XX(K)=(XX(K)+RXX(K))/2.0                                          
      DO 7 K=1,NG                                                       
      IF (DABS(XX(K)-RXX(K)).GT.1.0D-09) GOTO 19                       
   7  CONTINUE                                                          
  20  RETURN                                                            
      END                                                               
!**************************************************************************************************
     SUBROUTINE LVE_flash                                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER I                                                         
      CHARACTER*9 EX1(16),EX2(16)
      REAL*8 N2,SL,SS,WG,WK,PRES,X1R(150),X2R(150),RE,N1,beta,bmin,bmax,g,dg,AA,BB,VL
      INTEGER NK,NG,JX,IJ,k,sw                                          
      REAL*8 ZR,QT,RT,T,TT,P,CP,CV,CO,RO,X(150),X1(150),X2(150),ZZ,PP,Pcur,XX(10),Y(150),GAM1(150),GAM2(150),R(150),RS(150,10),RST(150,10),Ztot(150)
      Real*8 Q(150),QS(150,10),KD(150),KDLn(150),nL1(150),nL2(150),XL1(150),XL2(150)         
      COMMON /A89/ X,X1,X2,XX,Y,GAM1,GAM2                               
      COMMON /C89/ QT,RT,T,TT,R,RS,RST,Q,QS                             
      COMMON /D89/ ZZ,ZR,NK,NG,JX,IJ                                    
 !   open(22,file='slosh.dat')
      NT=NK-1                                                           
	EX1(1)='BuOH'
	EX1(2)='H2O'
	CO=14.244d-6
	sw=0
 ! initialize phase compositions and pressure
 	PP=1.01325d5									  !	5	7	  8
	Pcur=1.0d0*1.0d6							  !	6	8	  7	 10.0d0*1.0d6

	X1(NK) = 0.001d0          !0.20				  ! 0.001	 0.2
      RE=1.0d0-X1(NK)
	X1(1) = 0.1d0*RE			  ! new
	X1(2)=0.8d0*RE
	X1(3)=RE-X1(1)-X1(2)

	X2(NK) = 0.001d0           !0.30				  ! 0.001	 0.3
      RE=1.0d0-X2(NK)
	X2(1) = 0.8d0*RE			  ! new
	X2(2)=0.1d0*RE
 	X2(3)=RE-X2(1)-X2(2)

      Ztot(NK)=X1(NK)
!	Ztot(1)=X1(1)+0.5d0				 !	the added constant value (0.5;0.2;0.8) controls the gross composition stored in X(1:2), see below
!	Ztot(2)=
!!	Ztot(3)=1.0d0-Ztot(1)-Ztot(2)-Ztot(NK)

 ! set the gross fraction of holes to pressure
 35	P=Pcur
      do i=1,nt
		Ztot(i)=X(i)*(1.0d0-Ztot(nk))
	enddo

      CALL ACTCF(Ztot,GAM2)

	  RE = dexp(-P*CO/(8.3143*T))/GAM2(NK)
	  if(Dabs(Ztot(nk)-RE).gt.1.0d-8) then
	  Ztot(NK)=RE
	  goto 35			 
	  endif

      PRES=1.0D-10                                                      
      IJ=.0                                                             
      CALL ACTCF(X1,GAM1)                                               
      CALL ACTCF(X2,GAM2)                                                
	 SL=0.0d0
  19  IJ=IJ+1
      SS=0.0d0                                                          
      DO 18 I=1,NK
	KD(I)=GAM1(I)/GAM2(I)                                                      
  18  SS=SS+DABS(X1(I)*GAM1(I)-X2(I)*GAM2(I))                           

!	 set the boundaries on beta - the phase fraction 
	bmin=1.0d-10
	do i=1,nk
	if (KD(i).gt.1.0d0) then
	RE= (KD(i)*Ztot(i)-1.0d0)/(KD(i)-1.0d0)
	if (RE.gt.bmin) bmin=RE
	endif
	enddo
      bmax=1.0d0-1.0d-10
	do i=1,nk
	if (KD(i).lt.1.0d0) then
	RE= (1.0d0-Ztot(i))/(1.0d0-KD(i))
	if (RE.lt.bmax) bmax=RE
	endif
	enddo
	beta=(bmin+bmax)/2.0d0
! iterate over beta
	if (g(beta,Ztot,KD, nk)*g(bmin,Ztot,KD, nk).lt.0) beta=bmin
51	RE = beta-g(beta,Ztot,KD, nk)/dg(beta,Ztot,KD, nk)
	if (dabs(RE-beta).gt.(pres*1.0d-3))then
	beta=RE
	goto 51
	endif
! calculate feed flow (nL1,nL2) and check material balance for feed and fractions in the phases - X1&X2
      RE=0.0d0
	N1=0.0d0
	N2=0.0d0 
      do i=1,nk
	nL1(i)=(1.0d0-beta)*Ztot(i)/(1.0d0-beta+beta*KD(i))
	nL2(i)=beta*KD(i)*Ztot(i)/(1.0d0-beta+beta*KD(i))
	N1=N1+nL1(i)
	N2=N2+nL2(i) 
	RE=RE+Ztot(i)-nL1(i)-nL2(i)
	enddo
	RE=0.0d0
	do i=1,nk
	X1(i)=nL1(i)/N1
	X2(i)=nL2(i)/N2
	RE=RE+Ztot(i)-(1.0d0-beta)*X1(i)-beta*X2(i)
	enddo

      CALL ACTCF(X1,GAM1)                                               
      CALL ACTCF(X2,GAM2)                                                

	if (ij.gt.300) then
      PRINT 52                                                          
	write (6,52) SS
  52  FORMAT (2X,'NUMBER OF ITER MORE THAN 300',E10.4)                  
	return
	endif

	print *, ss,sl
 	if (DABS(SS-SL).gt.pres)then
	SL=SS 
	goto 19
	endif

	P=- 8.3143*T*DLOG(GAM1(NK)*X1(NK))/CO

      if(dabs(1.0d0-P/PP).lt.1.0d-5)goto 78
	if(sw.eq.0)then
	sw=1
	yy1=PP-P
	xx1=Pcur
	Pcur=PP*0.5d0
	goto 35
	endif
	if(sw.eq.1)then
	sw=2
	yy2=PP-P
	xx2=Pcur
      AA=(yy1-yy2)/(xx1-xx2)
	BB= (yy1+yy2)/2.0d0-AA*(xx1+xx2)/2.0d0
	Pcur=-BB/AA
	goto 35
	endif
	if(sw.gt.1)then
	sw=sw+1
	yy1=yy2
	xx1=xx2
	yy2=PP-P
	xx2=Pcur
      AA=(yy1-yy2)/(xx1-xx2)
	BB= (yy1+yy2)/2.0d0-AA*(xx1+xx2)/2.0d0
	Pcur=-BB/AA
	goto 35
	endif
 78	RE=-8.3143d-6*T/CO
      P=RE*DLOG(GAM1(NK)*X1(NK))
 
	  write(22,'(/,1X,"ss= ",E10.4,1X,"sl= ",E10.4,1X)') ss, sl
      write(22,'(1X,"ij=",I4,2X,"sw=",I4,2X,2(E16.10,1X))') ij,sw,P,(RE*DLOG(GAM2(NK)*X2(NK)))	  
      write(22,'(/,1X,"X1,X2,Z",12(E14.8,1X))')X1(1:3),X2(1:3),Ztot(1:3)
!  output for equlibrium liquid phases
      write (22, 1) T                                                   
  1   FORMAT (/38X,'TEMPERATURE  =  ',F7.2)                             
      write (22, 2)                                                     
  2   FORMAT (' ',73X,'MOL.VOLMS./COMPR.FACTORS',/,16X,'Z(gross)',12X,'XL1',12X,'XL2',9X,'P(MPA)',8X,'PHASE1',5X,'PHASE2')

      AA2=.0d0                                                          
      AA3=.0d0                                                          
      DO 889 I=1,nk                                                    
      AA3=AA3+R(I)*X1(I)                                                
  889 AA2=AA2+R(I)*X2(I)                                                
      VL1=(1.0d6)*CO*AA3/(1.0-X1(NK))                                 
      VL2=(1.0d6)*CO*AA2/(1.0-X1(NK))                                  
      ZL1=(1.0d-6)*VL1*P/(T*8.3143d0)                                        
      ZL2=(1.0d-6)*VL2*P/(T*8.3143d0)                                          

!      P1=P*1.0d-6
!     P1=P*9.8692E-06                                                   
      write (22, 3) EX1(1),X(1),X1(1)/(1-X1(NK)),X2(1)/(1-X2(NK)),P,VL1,VL2,EX1(2),X(2),X1(2)/(1-X1(NK)),X2(2)/(1-X2(NK)),ZL1,ZL2
  3   FORMAT (A8,6X,3(E14.8,1X),3(E12.6,1X)/A8,6X,3(E14.8,1X),13X,2(E12.6,1X) /A8,6X,3(E14.8,1X))


	  close (22)

      RETURN                                                            
      END                                                               
!**************************************************************************************************
!---------------------------------------------------------------------------
      real*8 function g(beta, Z, KD, nk)
	implicit none
	integer i, nk
	real*8 beta, Z(150),KD(150)

	g=0.0d0
	do i=1,nk
	 g =g+Z(i)*(KD(i)-1.0d0)/(1.0d0-beta+beta*KD(i))
	enddo
	return
	end
!--------------------------------------------------------------------------      
!---------------------------------------------------------------------------
      real*8 function dg(beta, Z, KD, nk)
	implicit none
	integer i, nk
	real*8 beta, Z(150),KD(150)

	dg=0.0d0
	do i=1,nk
	dg =dg-Z(i)*((KD(i)-1.0d0)/(1.0d0-beta+beta*KD(i)))**2
	enddo
	return
	end
!--------------------------------------------------------------------------
      real*8 function VL(X,nk)
	implicit none
	integer i, nk 
	real*8 QT,RT,T,TT,R(150),Q(150),RS(150,10),RST(150,10),QS(150,10)
	real*8 X(10),AA
      COMMON /C89/ QT,RT,T,TT,R,RS,RST,Q,QS 
      AA=0.0d0
	do i=1,nk 
      AA=AA+R(I)*X(I)
	enddo
	VL=AA/(1.0-X(NK))
	return 
	end      
!--------------------------------------------------------------------------


!**************************************************************************************************
 SUBROUTINE XcomptoXseg(Xcomp,Xseg)
      IMPLICIT NONE
	  character*7 name_seg(150)
	  real*8 ZZ,ZR,JX,IJ,Xcomp(150),Xseg(150)
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  integer b,M,NK,NG,Nmol,I,chem(150),type_mol(150),Nbonds,rank_seg(150)
	  REAL*8 teta_chem(150,150),sum
sum=0.0d0
b=0
	  DO I=1,NK
	  !num_seg(I)=num_segg(I)	  !real from integer
	   if (Nmol.ne.0) then
	   do M=1,Nmol
	     if (type_mol(I).eq.M) then 
	    Xseg(I)=Xcomp(M) !*num_seg(I)
	    endif
	   enddo  !M
	   if (type_mol(I).eq.0) then
		b=b+1
		Xseg(I)=Xcomp(b+Nmol) !*num_seg(I)
	    endif
	   else
	  Xseg(I)=Xcomp(I) !*num_seg(I)
	   endif
	  enddo
	  DO I=1,NK 
	  sum=sum+Xseg(I)
	  enddo 
	  DO I=1,NK 
	  Xseg(I)=Xseg(I)/sum
	  enddo

Xseg(nk)=0.0d0
do I=1,NK-1
Xseg(nk)=Xseg(nk)+Xseg(I)
enddo
Xseg(nk)=1.0d0-Xseg(nk)

 return
 end
 !**************************************************************************************************
 real*8 FUNCTION gauss(x,mediana,height,sigma,left_border,right_border)
IMPLICIT none
REAL*8	x,mediana,height,sigma,left_border,right_border

gauss=dexp(-0.5d0*((x-mediana)/sigma)**2)
if (x.le.mediana) gauss=gauss*(height-left_border)+left_border 
if (x.gt.mediana) gauss=gauss*(height-right_border)+right_border
 END FUNCTION
!C********************************************************************************
integer FUNCTION seg_ch(bond)				!min number of segment with chem number=bond
IMPLICIT none
character*7 name_seg(150)
integer bond,I,k,Nmol,NK,NG,M,Nbonds, J 
real*8 	teta_chem(150,150)
REAL*8 ZZ,ZR,JX,IJ
 common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
integer chem(150),type_mol(150),rank_seg(150), rank_branch(150), n_branch, type_branch(150)
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
common /branch/rank_branch, n_branch, type_branch
seg_ch=0 
k=0
do M=1,Nmol	
DO  I=1,NK
if ( (type_mol(I).eq.M).and.(rank_seg(I+1).gt.1) ) then
k=k+1
if (k.eq.bond) seg_ch=I
endif
ENDDO

DO  I=1,NK
if	( (type_mol(I).eq.M) .and. (type_branch(I) .ne. 0) .and. (rank_branch(I) .eq. 1) ) then !(chem(J) .eq. 3)
	do J = 1, NK
	if ( (type_branch(J) .eq. I) .and. ( (rank_branch(J) .eq. 1) .or.  (type_branch(J+1) .eq. I) ) ) then 
	k=k+1
	if (k .eq. bond) seg_ch=J 
	endif
	enddo
endif
enddo
enddo !(M)

 END FUNCTION
!C********************************************************************************
 integer FUNCTION seg_max_ch(bond)				!max number of segment with chem number=bond
IMPLICIT none
character*7 name_seg(150)
integer bond,I,k,Nmol,NK,NG,M,Nbonds, J 
real*8 	teta_chem(150,150)
REAL*8 ZZ,ZR,JX,IJ
 common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
integer chem(150),type_mol(150),rank_seg(150), rank_branch(150), n_branch, type_branch(150)
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
common /branch/rank_branch, n_branch, type_branch
seg_max_ch=0 
k=0
do M=1,Nmol	
DO  I=1,NK
if ( (type_mol(I).eq.M).and.(rank_seg(I+1).gt.1) ) then
k=k+1
if (k.eq.bond) seg_max_ch=I+1
endif
ENDDO

DO  I=1,NK
if	( (type_mol(I).eq.M) .and. (type_branch(I) .ne. 0) .and. (rank_branch(I) .eq. 1) ) then !(chem(J) .eq. 3)
	do J = 1, NK
	if ( (type_branch(J) .eq. I) .and.  (rank_branch(J) .ne. 1) ) then 
	k=k+1
	if (k .eq. bond) seg_max_ch=J 
	endif
	enddo
endif
enddo
enddo !(M)

 END FUNCTION
!C********************************************************************************
subroutine total(X,QT,RT) !,NK,Q,R,QT,RT)   !QT, RT - output
IMPLICIT NONE
REAL*8 X(150) !,QT,RT,Q(150),R(150)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds  
integer k,NK,NG,Nmol,Nbonds
QT=0.0d0
RT=0.0d0
 DO k=1,NK
 QT=QT+Q(k)*X(k)
 RT=RT+R(k)*X(k)
 enddo
return
 end
 !C********************************************************************************
subroutine local_chem(X,L,QT) !,Nbonds,L,QT)			   !writes teta_chem in common
IMPLICIT NONE
character*7 name_seg(150)
integer L,k,NK,NG,Nmol,Nbonds,seg_ch,chem(150),type_mol(150),rank_seg(150)
real*8 	teta_chem(150,150),X(150,150)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds  

do k=1,Nbonds	
teta_chem(k,L)=X(seg_ch(k),L)/QT
ENDDO
 return
 end
!C********************************************************************************
real*8 function Rout(L)
IMPLICIT NONE
CHARACTER*3 shape
common /layer/shape,TL,Lmax,TLmin 
integer L,LL,Lmax                           !L - layer
real*8 TL(150),TLmin
 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
 Rout=r_solid !+TL(Lmax+1) !0.0d0
 if (L.gt.1) then
do LL=1,L-1
Rout=Rout-TL(LL)						!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
enddo
endif
 return
 end
!C********************************************************************************
real*8 function m_up(L)
IMPLICIT NONE
CHARACTER*3 shape
integer L,LL,I,NG,Nmol,NK, khom,Nbonds,Lmax                           !L - layer
real*8 Ntot,Vol,pi,TL(150),volume,TLmin    
REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150)
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ,XlocL(150),Tprev,Tcurrent
common /layer/shape,TL,Lmax,TLmin                              !TL - thickness of layer
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                
COMMON /C/T,TT,R,RS,Q,QS,vol_ref, q_ref 
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
common /kind/khom
REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
!LL=L
pi = 3.1415926d0
 Tprev=r_solid !0.0d0
 if (L.gt.1) then
do LL=1,L-1
Tprev=Tprev-TL(LL)				   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
enddo
endif
Tcurrent=Tprev-TL(L)			   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 if ((r_solid-Tcurrent).gt.1.0d-50) then
 DO I=1,NK                                                       
 if (khom /= 0) then
  XlocL(I)=X(I) 
  else
  XlocL(I)=Xloc(I,L)                                                 
  endif                                                 
 enddo  
 call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output

Vol=volume(L)
Ntot=Vol/RT
if (shape.eq.'sph')   m_up=4.0d0*pi*(Tcurrent**2)/(QT*Ntot*6.0d0)
if (shape.eq.'pln')  m_up=Vol/TL(L)/(QT*Ntot*6.0d0)
!if (shape.eq.'sph')   m_up=4.0d0*pi*(TL**2)*(LL**2)/(QT*Ntot*6.0d0) !(q_ref*QT*Ntot)                 !volume of layer
!if (shape.eq.'cyl') m_up=4.0d0*pi*(TL**2)*(LL**2)/(QT*Ntot*6.0d0) !(q_ref*QT*Ntot)
!if (shape.eq.'pln')  m_up=1.0d0/TL(L)/6.0d0 !6.0d0					 
!m_up=m_up*QT*Ntot
!DO I=1,NK                                                       
! if (khom) then
!  XlocL(I)=X(I) 
!  else
!  XlocL(I)=Xloc(I,L)                                                 
!  endif                                                 
! enddo  
! call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
!m_up=m_up/QT
else
m_up=0.0d0
endif
return
end
!******************************************************************************************************************
real*8 function m_dwn(L)
IMPLICIT NONE
CHARACTER*3 shape
integer L,LL,I,NG,Nmol,NK,khom,Nbonds,Lmax                            !L - layer
real*8 Ntot,Vol,pi,TL(150),volume,TLmin    
REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150)
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ,XlocL(150),Tprev
common /layer/shape,TL,Lmax,TLmin                                !TL - thickness of layer
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref 
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
common /kind/khom
REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
pi = 3.1415926d0
Tprev=r_solid !0.0d0
 if (L.gt.1) then
do LL=1,L-1
Tprev=Tprev-TL(LL)						!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
enddo
endif
!LL=L                                                        
 DO I=1,NK                                                       
 if (khom /= 0) then
  XlocL(I)=X(I) 
  else
  XlocL(I)=Xloc(I,L)
  endif                                                 
 enddo  
 call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
Vol=volume(L)

Ntot=Vol/RT
if (shape.eq.'sph')   m_dwn=4.0d0*pi*(Tprev**2)/(QT*Ntot*6.0d0) 
if (shape.eq.'pln')  m_dwn=Vol/TL(L)/(QT*Ntot*6.0d0)
!if (shape.eq.'sph')   m_dwn=4.0d0*pi*(TL**2)*((LL-1)**2)/(QT*Ntot*6.0d0) !(q_ref*QT*Ntot)                 !volume of layer
!if (shape.eq.'cyl') m_dwn=4.0d0*pi*(TL**2)*((LL-1)**2)/(QT*Ntot*6.0d0) !(q_ref*QT*Ntot)
!if (shape.eq.'pln')  m_dwn=1.0d0/TL/6.0d0 !6.0d0

if ((L.eq.1).and.(r_solid.lt.1.0d-50))  m_dwn=0.0d0
return
end
!******************************************************************************************************************
real*8 function lat(L)
integer L                            !L - layer
real*8 m_dwn,m_up
lat=1.0d0-m_dwn(L)- m_up(L)
!if (lat.lt.0.0d0) write(*,*)'Lat<0:',lat
return
end
!******************************************************************************************************************
real*8 function volume(L)				   !volume of layer L
CHARACTER*3 shape
integer L,Lmax,layer                            !L - layer
real*8 LL,TL(150),pi,Tcurrent,Tprev,TLmin  
common /layer/shape,TL,Lmax,TLmin
REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid  
pi = 3.1415926d0
Tprev=r_solid !0.0d0
do layer=1,L !-1
Tprev=Tprev-TL(layer)
enddo
Tcurrent=Tprev+TL(L)
!LL=L
if (shape.eq.'sph')  volume=4.0d0*pi/3.0d0*(Tcurrent**3-Tprev**3)
if (shape.eq.'pln')  volume=TL(L) 
!if (shape.eq.'sph')  volume=4.0d0*pi/3.0d0*(TL**3)*(LL**3-(LL-1)**3) 
!if (shape.eq.'cyl')  volume=4.0d0*pi*(TL**2)*(LL**2-(LL-1)**2)
!if (shape.eq.'pln')  volume=TL     !a_pl*TL
return
end
!******************************************************************************************************************
  SUBROUTINE calcALloc(PAR,L,ALloc,sumX)
      IMPLICIT NONE
	  character*7 name_seg(150)
	  real*8 X(150),ALloc(150),PAR(150),ZZ,ZR,JX,IJ,ALchem(150)
	  COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  integer NK,NG,Nmol,I,K,S,chem(150),type_mol(150),Nbonds,rank_seg(150),L
	  REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
	  REAL*8 teta_chem(150,150),sumX(150)
   X(NK)=1.0d0
	    DO  I=1,NK-1                                                  
        X(I)=PAR(I)
		X(NK)=X(NK)-X(I)
		enddo
   call total(X,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
  DO s=1,NG
  ALloc(s)=0.0d0
  DO k=1,NK
  ALloc(s)=ALloc(s)+QS(k,s)*X(k)/QT
  enddo 
  enddo
  call calcALchem(PAR,L,ALchem,sumX)

  DO I=1,Nbonds
  ALloc(I+NG)=ALchem(I)
  enddo
  return
  end
  !******************************************************************************************************************
    SUBROUTINE calcALchem(PAR,L,ALloc,sumX)
      IMPLICIT NONE
	  character*7 name_seg(150)
	  real*8 X(150),ALloc(150),PAR(150),ZZ,ZR,JX,IJ,m_up
	  CHARACTER*3 shape
	  COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
	  COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
	  common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	  integer NK,NG,Nmol,I,K,S,chem(150),type_mol(150),rank_seg(150),Nbonds,Lay,L,Lmax,seg_ch, seg_max_ch
	  REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
	   REAL*8 teta_chem(150,150),sumX(150),volume,TLmin  
	   REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
	   real*8 Vol,TL(150),VolLL,RTLL
	  COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
	  common /layer/shape,TL,Lmax,TLmin  
	
   X(NK)=0.0d0
	    DO  I=1,NK-1                                                  
        X(I)=PAR(I)
		X(NK)=X(NK)+X(I)
		Xloc(I,L)=X(I)
	!	print *,'l',Xloc(I,L)
		enddo
		X(NK)=1.0d0-X(NK)
		Xloc(NK,L)=X(NK)

    call total(X,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
	Vol=volume(L)
		do k=1,Nbonds
		sumX(k)=0.0d0
		enddo

do k=1,Nbonds
ALloc(k)=X(seg_ch(k))/QT/ZZ
teta_chem(k,L)=ALloc(k)*ZZ
enddo 

do k=1,Nbonds
	do Lay=1,L
	if (Lay.eq.L) then
	RTLL=RT
	VolLL=Vol
	else
	VolLL=volume(Lay)

	RTLL=0.0d0
	DO s=1,NK
	RTLL=RTLL+R(s)*Xloc(s,Lay)
	enddo
   endif
	sumX(k)=sumX(k)+(Xloc(seg_max_ch(k),Lay)-Xloc(seg_ch(k),Lay))/ZZ/RTLL*VolLL
		enddo
	sumX(k)=sumX(k)/QT/m_up(L)*RT/Vol
	enddo 
	!print *,'jjjjjj1',sumX
!		if (L.eq.21) print *,sumX(11)
  return
  end
  !******************************************************************************************************************
   SUBROUTINE calcAL_al_and_NG(L)
   IMPLICIT NONE
   character*7 name_seg(150)
   integer L,M,NK,NG,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150)
   integer s,KEY,KEYs(150),KEYall,KEYngup,KEYngdwn,seg_ch,seg_max_ch
   real*8 m_up,m_dwn,lat,ZZ,ZR,JX,IJ,sumup,sumdwn
   REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150)
   REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150)   
   common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
   common /chem_or/chem_up,ALchemold
   common /chemdwn/chem_dwn,chem_al
   COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
   COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
   COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
   common /for_check_depend/KEY
   common /key_for_minmax/KEYs,KEYall,KEYngup,KEYngdwn
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid

  KEY=0
  	KEYngup=0
   KEYngdwn=0
  !print *,'AL_al',ALloc(1,L),m_up(L)/lat(L),(ALloc(1,L)-AL_up(1,L)),m_dwn(L)/lat(L),(ALloc(1,L)-AL_dwn(1,L))

  sumup=0.0d0
  sumdwn=0.0d0
  DO  s=1,NG-1
  AL_al(s,L)=ALloc(s,L)+m_up(L)/lat(L)*(ALloc(s,L)-AL_up(s,L))+m_dwn(L)/lat(L)*(ALloc(s,L)-AL_dwn(s,L))
	if (AL_al(s,L).le.0.0d0) then
!	print *,'AL_al small',s
	!AL_al(s,L)=1.0d-100
	KEYs(s)=-1
	endif	   
  enddo
do s=1,NG-1
sumup=sumup+ALloc(s,L)-AL_up(s,L)
if ((L.ne.1).or.(r_solid.gt.1.0d-50))  sumdwn=sumdwn+ALloc(s,L)-AL_dwn(s,L)
enddo
! print *, 'd',sumdwn,ALloc(NG,L)+sumdwn,AL_dwn(1:(NG-1),L),ALloc(NG,L)
  if (Nmol /= 0) then
	DO  M=1,Nbonds
	 chem_al(M,L)=teta_chem(M,L)/ZZ+m_up(L)/lat(L)*(teta_chem(M,L)/ZZ-chem_up(M,L))+m_dwn(L)/lat(L)*(teta_chem(M,L)/ZZ-chem_dwn(M,L)) 
	if (chem_al(M,L).le.0.0d0) then
!	print *,'chem_al small',M,L,chem_al(M,L),teta_chem(M,L)/ZZ,chem_up(M,L),chem_dwn(M,L)
	!!stop
	!chem_al(M,L)=1.0d-150
	!KEY=1
	endif
	sumup=sumup+teta_chem(M,L)/ZZ-chem_up(M,L)
	sumup=sumup+teta_chem(M,L)/ZZ*Xloc(seg_max_ch(M),L)/Xloc(seg_ch(M),L)-chem_dwn(M,L+1)
if ((L.ne.1).or.(r_solid.gt.1.0d-50))	sumdwn=sumdwn+teta_chem(M,L)/ZZ-chem_dwn(M,L)
if ((L.ne.1).or.(r_solid.gt.1.0d-50))	sumdwn=sumdwn+teta_chem(M,L)/ZZ*Xloc(seg_max_ch(M),L)/Xloc(seg_ch(M),L)-chem_up(M,L-1)
	AL_up(NG,L)=AL_up(NG,L)+chem_up(M,L)+chem_dwn(M,L+1)

 if ((L.ne.1).or.(r_solid.gt.1.0d-50)) AL_dwn(NG,L)=AL_dwn(NG,L)+chem_dwn(M,L)+chem_up(M,L-1)
	if ((L.eq.1).and.(r_solid.lt.1.0d-50))  AL_dwn(NG,L)=AL_dwn(NG,L)+chem_dwn(M,L)*2.0d0
	AL_al(NG,L)=AL_al(NG,L)+chem_al(M,L)*2.0d0
	enddo
  endif
AL_up(NG,L)=ALloc(NG,L)+sumup
if ((L.ne.1).or.(r_solid.gt.1.0d-50))  AL_dwn(NG,L)=ALloc(NG,L)+sumdwn
if ((L.eq.1).and.(r_solid.lt.1.0d-50))  	AL_dwn(NG,L)=0.0d0
	if ((r_solid-L).lt.1.0d-100) AL_up(NG,L)=0.0d0

AL_al(NG,L)=ALloc(NG,L)+m_up(L)/lat(L)*(ALloc(NG,L)-AL_up(NG,L))+m_dwn(L)/lat(L)*(ALloc(NG,L)-AL_dwn(NG,L))
	 	if (AL_up(NG,L).le.0.0d0) then
		!	print *,'AL_up NG small',AL_up(NG,L),L
		KEYngup=-1
		! KEY=1
		! stop								    
		 endif														 
		
		if (AL_al(NG,L).le.0.0d0) then
	!	print *,'AL_al NG small'
		!AL_al(NG,L)=1.0d-150 				 
		KEYall=1
		endif
		if ((AL_dwn(NG,L).le.0.0d0).and.((L.ne.1).or.(r_solid.gt.1.0d-50))) then
	!	print *,'AL_dwn NG small',AL_dwn(NG,L),L
		KEYngdwn=-1
		endif 
	  ! KEY=0
	 !	print *,'AL_up NG',AL_up(NG,1)
  return							  
  end
  !******************************************************************************************************************
  SUBROUTINE check_chemdwn(PAR,L,KEY)
   IMPLICIT NONE
   CHARACTER*3 shape
   character*7 name_seg(150)
   integer L,LL,Lay,NK,NG,Nmol,Nbonds,I,k,s,Lmax,KEY,chem(150),type_mol(150),rank_seg(150),phase,seg_ch, seg_max_ch
   real*8 m_up,ZZ,ZR,JX,IJ,QT,TL(150),RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref,VolLL,RTLL,TLmin  
  REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150),ALchemold(150,150)
   REAL*8 teta_chem(150,150),PAR(150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150)  
   common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
   real*8 Vol,pi,volume,XlocL(150)
   common /chem_or/chem_up,ALchemold
   common /chemdwn/chem_dwn,chem_al
  COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
   COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
   common /layer/shape,TL,Lmax,TLmin  
   COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
    COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
	REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
	common /ph/phase
   pi = 3.1415926d0
  KEY=0
    Xloc(NK,L)=1.0d0
	    DO  I=1,NK-1                                                  
        Xloc(I,L)=PAR(I)
		Xloc(NK,L)=Xloc(NK,L)-Xloc(I,L)
		enddo

DO LL=L,Lmax
  DO k=1,NK
  XlocL(k)=Xloc(k,LL)
  ENDDO
  call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output

  Vol=volume(LL)
  do k=1,Nbonds
  chem_dwn(k,LL+1)=0.0d0
  enddo

	do k=1,Nbonds	
		do Lay=1,LL
		if (Lay.eq.LL) then
	RTLL=RT
	VolLL=Vol
	else
  	VolLL=volume(Lay)

	RTLL=0.0d0
	DO s=1,NK
	RTLL=RTLL+R(s)*Xloc(s,Lay)
	enddo
   endif
		chem_dwn(k,LL+1)=chem_dwn(k,LL+1)+(Xloc(seg_max_ch(k),Lay)-Xloc(seg_ch(k),Lay))/ZZ/RTLL*VolLL
		enddo
  	chem_dwn(k,LL+1)=chem_dwn(k,LL+1)/QT/m_up(LL)*RT/Vol+chem_up(k,LL)
!	if (LL.eq.5) print *,'lllllllllllllll',chem_dwn(5,6)

	if (chem_dwn(k,LL+1).le.0.0d0) KEY=1
	enddo !Nbonds
ENDDO	  !LL
!call calcAL_al_and_NG(L,phase)
!do s=1,NG
!if ((AL_al(s,L).eq.1.0d-150).or.(AL_up(s,L).eq.1.0d-150).or.(AL_dwn(s,L).eq.1.0d-150) )   KEY=1
!enddo
if (L.eq.Lmax) Key=0
!Key=0

  return
  end
  !******************************************************************************************************************
  SUBROUTINE check_dependALF(PARth,L,KEY)
  IMPLICIT NONE
  CHARACTER*3 shape
  character*7 name_seg(150)
  integer L,M,NK,NG,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150)
   integer k,I,Nx,Nalf,LL,KEY(150),Lmax
   real*8 ZZ,ZR,JX,IJ,PARth(150),DPx(150),DPalf(150),EPS,EPF,TL(150),XlocL(150),sumX(150)
   REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150)
   REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150),TLmin  
   REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)   
   common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
   common /chem_or/chem_up,ALchemold
   common /chemdwn/chem_dwn,chem_al
   COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
   COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
   COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
   common /forNEWT/Nx,Nalf,DPx,DPalf,EPS,EPF
   common /layer/shape,TL,Lmax,TLmin
   REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
  
   !common /for_check_depend/KEYY

   do M=1,Nbonds
  KEY(M)=0
  enddo
  
		DO  I=NK,Nx+Nalf/2		 
	    AL_up(I-Nx,L)=PARth(I)									
		enddo
 		if ((L.eq.1).and.(r_solid.lt.1.0d-50)) then
			DO  I=1,NG
			AL_dwn(I,L)=0.0d0
			enddo
		else
		DO  I=NK+Nalf/2+Nbonds,Nx+Nalf+Nbonds
		AL_dwn(I-Nx-Nalf/2-Nbonds,L)=PARth(I)
	    enddo
		endif

  if (Nmol /= 0) then
		do M=1,Nbonds
		chem_up(M,L)=PARth(Nx+Nalf/2+M)
		enddo

DO LL=L,Lmax-1
  DO k=1,NK
  XlocL(k)=Xloc(k,LL)
  ENDDO
 if ((LL.eq.1).and.(r_solid.lt.1.0d-50)) then	   
		DO k=1,Nbonds
		 chem_dwn(k,LL)=0.0d0 
        enddo	
 endif
call calcALchem(XlocL,LL,Xchem,sumX)
do k=1,Nbonds
chem_dwn(k,LL+1)=sumX(k)+chem_up(k,LL)
if (chem_dwn(k,LL+1).le.0.0d0) KEY(k)=1
!if (chem_dwn(k,LL+1).le.0.0d0) print *,sumX(k),LL,chem_up(k,LL)
enddo

ENDDO	  !LL
endif	  !Nmol

! DO LL=L,Lmax
 !call calcAL_al_and_NG(LL)
 !if (KEYY) KEY=1
!ENDDO  !LL

  return
  end
  !******************************************************************************************************************
  SUBROUTINE check_borderALF(PARth,L,Nsolve,ALlocnew,sumXX,id,mainK,switch)
  IMPLICIT NONE
  integer L,NK,NG,Nmol,Nbonds,CHANGE(150),Nsolve,Nfix,No_or(150)
   integer kk,I,Nx,Nalf,D,switch,allK,sup,chemup,sdwn,KEYid(100),id(100),KEYch(100),mainK
   real*8 ZZ,ZR,JX,IJ,PARth(150),PAR(150),ALlocnew(150),DPx(150),DPalf(150),EPS,EPF,sumXX(150)
   real*8 alf_max(150),alf_min(150)
   REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)    
   COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
   COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
   common /forNEWT/Nx,Nalf,DPx,DPalf,EPS,EPF
   common /without_or/ Nfix,No_or,CHANGE
   common /min_maxALF/alf_max,alf_min
   !switch=1 - went from initial profile, =0 - corrections   PARth, =-1 - checking   PARth after interation in NEWT
   kk=0
   if (switch.gt.0) then	!from initial profile
   do D=1,Nx
   PARth(D)=Xloc(D,L)
   enddo
   endif
   if (switch.ge.0) then
	call minmaxALF(PARth,Nx,Nalf,Nsolve,ALlocnew,sumXX,L,KEYid,KEYch,sup,chemup,sdwn,allK,1)
	call id_to_change(Nsolve,Nbonds,allK,sup,chemup,sdwn,KEYid,KEYch,id,mainK,'before iter in NEWTalf')

	 do D=1,Nsolve+Nbonds
	 PAR(D)=PARth(CHANGE(D))
	 enddo 

	 do while (mainK.ne.0)
	 kk=kk+1
	! if (L.eq.17) print *,sup,chemup,sdwn,allK
!	 if (L.eq.12) print *,PAR(1:15) !id(1:15)
!	if (L.eq.17) print *,KEYid(1:15)
!	if (L.eq.16) stop
	  	if (kk.gt.100000) then
!if (L.eq.15) then
			print *,'many steps for initial Alfa in NEWT',PAR(1)
			if (allK.ne.0) print *,'allK=',allK
			if (sup.ne.0) print *,'sup=',sup
			if (chemup.ne.0) print *,'chemup=',chemup
			if (sdwn.ne.0) print *,'sdwn=',sdwn
			write (*,'(6(5X,A7))')'1...Nsolve+Nbond','id','KEYid'
			 DO I=1,Nsolve+Nbonds
			if (id(I).ne.0) print *,i,id(I),KEYid(I)
			enddo
	!		stop
		mainK=1
		return
			endif

	 DO I=1,Nsolve+Nbonds
	 PAR(I)=PARth(CHANGE(I))
	 if (PAR(I).le.(alf_min(I)+1.0d-150)) then
	 PAR(I)=alf_min(I)+0.001d0*(alf_max(I)-alf_min(I))
	 elseif (PAR(I).ge.(alf_max(I)-1.0d-150)) then
	 PAR(I)=alf_max(I)-0.001d0*(alf_max(I)-alf_min(I))
	 else
		if (id(I).eq.1) then
		PAR(I)=PAR(I)*1.01d0
		elseif (id(I).eq.(-1)) then
		PAR(I)=PAR(I)*0.99d0 	
		endif
	endif
		PARth(CHANGE(I))=PAR(I)
	 enddo !I	

	 call minmaxALF(PARth,Nx,Nalf,Nsolve,ALlocnew,sumXX,L,KEYid,KEYch,sup,chemup,sdwn,allK,0)
	 call id_to_change(Nsolve,Nbonds,allK,sup,chemup,sdwn,KEYid,KEYch,id,mainK,'before iter in NEWTalf')
	 enddo !while

	elseif (switch.eq.(-1)) then
	call minmaxALF(PARth,Nx,Nalf,Nsolve,ALlocnew,sumXX,L,KEYid,KEYch,sup,chemup,sdwn,allK,1)	!0!!!!!!!!!!!!!!!!!!!
	call id_to_change(Nsolve,Nbonds,allK,sup,chemup,sdwn,KEYid,KEYch,id,mainK,'after iter in NEWTalf ')
	 endif !switch
  return
  end
  !******************************************************************************************************************
	SUBROUTINE	id_to_change(Nsolve,Nbonds,allK,sup,chemup,sdwn,KEYid,KEYch,id,mainK,switch)
	!output: id(D)=1 to increase D, =-1 to decrease D, =0 not to change D
	!allK=1 (increase all) or =0 (not matter)
	!sup=1 (increase all Al_up), =-1 (decrease) or 0
	!chemup 1,-1,0 (all chem_up)
	!sdwn 1,-1,0 (all Al_dwn)
	!mainK=0 when all id(D)=0, otherwise 1
	IMPLICIT NONE
	CHARACTER*22 switch
	integer	 D,Nsolve,Nbonds,allK,sup,chemup,sdwn,KEYid(100),id(100),KEYch(100),n,nup,ndwn,nchemup,mainK,KEYid_only,num_plus,num_min
!	integer
	KEYid_only=0
	num_plus=0
	num_min=0
	n=0
	nup=0
	ndwn=0
	nchemup=0
	mainK=0
!	if ( (allK.ne.0).and.((allK.eq.chemup).or.(allK.eq.sup).or.(allK.eq.sdwn)) )	 allK=0
	do D=Nsolve/2+1,Nsolve/2+Nbonds
	if ((KEYid(D).ne.0).or.(chemup.ne.0).or.(KEYch(D-Nsolve/2).ne.0)) KEYid_only=1
	enddo
	if (KEYid_only.eq.1) then
			do D=1,Nsolve+Nbonds
			id(D)=0
			if ((D.le.(Nsolve/2+Nbonds)).and.(D.gt.Nsolve/2)) then
				 if (KEYid(D).ne.0) then
				 id(D)=KEYid(D)
				 elseif (KEYch(D-Nsolve/2).ne.0) then
				 id(D)= KEYch(D-Nsolve/2)
				 elseif (chemup.ne.0) then
				 id(D)=chemup
				 endif
			if (id(D).ne.0)  mainK=1
			endif
			enddo
	else
!		enddo
!	else !all KEYid=0
!if (mainK.lt.0) then
	do D=1,Nsolve+Nbonds
	id(D)=0
	  if (KEYid(D).eq.1) then	 
	  id(D)=1
!	  num_plus=num_plus+1
	  elseif (KEYid(D).eq.(-1)) then
	  id(D)=-1
!	  num_min=num_min+1
!	  endif
!	 enddo
!	 	do D=1,Nsolve+Nbonds
!		if (KEYid(D).eq.0) then
	else
		if (D.le.Nsolve/2) then
		 if  (  (sup.eq.1).or.((allK.eq.1).and.(sup.eq.0).and.(num_min.ge.num_plus))  )	then
		 id(D)=1
	!	 num_plus=num_plus+1
	!	 sup=0
		 elseif (sup.eq.(-1)) then
		 id(D)=-1
	!	  num_min=num_min+1
	!	 sup=0
		 endif
	 !! elseif (D.le.(Nsolve/2+Nbonds)) then
	!!	if (  (chemup.eq.1).or.(KEYch(D-Nsolve/2).eq.1).or.((allK.eq.1).and.(chemup.eq.0).and.(KEYch(D-Nsolve/2).eq.0).and.(num_min.ge.num_plus))  )	then
	!!	id(D)=1
	!!	 num_plus=num_plus+1
	!	allK=0
	!chemup=0
	!!	elseif ((chemup.eq.(-1)).or.(KEYch(D-Nsolve/2).eq.(-1)) ) then
	!!	id(D)=-1
	!!	 num_min=num_min+1
	!	chemup=0
	!!	endif
	  elseif (D.gt.(Nsolve/2+Nbonds)) then   !D>Nsolve/2+Nbonds
	    if (  (sdwn.eq.1).or.((allK.eq.1).and.(sdwn.eq.0).and.(num_min.ge.num_plus))  )	then
		id(D)=1
	!	 num_plus=num_plus+1
		elseif (sdwn.eq.(-1)) then
		id(D)=-1
	!	 num_min=num_min+1
		endif
	  endif !D
	  endif !KEIid
!	  enddo
!	!  if ((allK.ne.0).and.(num_min.ge.num_plus))  then
!	  D=Nsolve/2+1
!	  do while ((D.le.(Nsolve/2+Nbonds)).and.(num_min.ge.num_plus)   )
!	  if (KEYid(D).ne.(-1)) then
!	  if (id(D).eq.(-1)) num_min=num_min-1
!	  if (id(D).ne.1) num_plus=num_plus+1
!	  id(D)=1
!	  endif
!	  D=D+1
	 if ((D.le.Nsolve/2).and.(sup.ne.0).and.(id(D).eq.sup)) nup=sup
	 if ((D.gt.Nsolve/2).and.(D.le.(Nsolve/2+Nbonds)).and.(chemup.ne.0).and.(id(D).eq.chemup)) nchemup=chemup
	 if ((D.gt.(Nsolve/2+Nbonds)).and.(sdwn.ne.0).and.(id(D).eq.sdwn)) ndwn=sdwn

	  if (id(D).ne.0)  mainK=1
	  if ((allK.ne.0).and.(id(D).eq.allK)) n=allK
		enddo
		endif
	
	if (mainK.lt.0) then
	do D=1,Nsolve+Nbonds
	
	id(D)=0
	  if (KEYid(D).eq.1) then	 
	  id(D)=1
	 ! num_plus=num_plus+1
	  elseif (KEYid(D).eq.(-1)) then
	  id(D)=-1
	!  num_min=num_min+1
	  else	! KEYid(D)=0
	  if (D.le.Nsolve/2) then
		 if  (  (sup.eq.1).or.((allK.eq.1).and.(sup.eq.0))  )	then
		 id(D)=1
	!	 sup=0
		 elseif (sup.eq.(-1)) then
		 id(D)=-1
	!	 sup=0
		 endif
	  elseif (D.le.(Nsolve/2+Nbonds)) then
		if (  (chemup.eq.1).or.(KEYch(D-Nsolve/2).eq.1).or.((allK.eq.1).and.(chemup.eq.0).and.(KEYch(D-Nsolve/2).eq.0))  )	then
		id(D)=1
	!	allK=0
	!chemup=0
		elseif ((chemup.eq.(-1)).or.(KEYch(D-Nsolve/2).eq.(-1)) ) then
		id(D)=-1
	!	chemup=0
		endif
	  else   !D>Nsolve/2+Nbonds
	    if (  (sdwn.eq.1).or.((allK.eq.1).and.(sdwn.eq.0))  )	then
		id(D)=1
		elseif (sdwn.eq.(-1)) then
		id(D)=-1
		endif
	  endif !D
	  endif !KEIid
	 if ((D.le.Nsolve/2).and.(sup.ne.0).and.(id(D).eq.sup)) nup=sup
	 if ((D.gt.Nsolve/2).and.(D.le.(Nsolve/2+Nbonds)).and.(chemup.ne.0).and.(id(D).eq.chemup)) nchemup=chemup
	 if ((D.gt.(Nsolve/2+Nbonds)).and.(sdwn.ne.0).and.(id(D).eq.sdwn)) ndwn=sdwn

	  if (id(D).ne.0)  mainK=1
	  if ((allK.ne.0).and.(id(D).eq.allK)) n=allK
	  enddo !D
	 endif
	if 	( (KEYid_only.eq.0).and. ((allK.ne.n).or.(sdwn.ne.ndwn).or.(chemup.ne.nchemup).or.(sup.ne.nup)) ) then
	print *,'different borders ALF ',switch
	if (allK.ne.n) print *,'allk',allk,n
	if (sdwn.ne.ndwn) print *,'sdwn',sdwn,ndwn
	if (chemup.ne.nchemup) print *,'chemup key',chemup,nchemup
	if (chemup.ne.nchemup) print *,'chem id=',(id(Nsolve/2+D),D=1,Nbonds)
	if (chemup.ne.nchemup) print *,'chem KEYid=',(KEYid(Nsolve/2+D),D=1,Nbonds)
	if (sup.ne.nup) print *,'sup',sup,nup
!	stop
	endif
   return
  end
  !******************************************************************************************************************
SUBROUTINE sum_of_chemup(PARth,Nx,Nalf,sumXX,L,chemup,KEYch)
 IMPLICIT NONE
 character*7 name_seg(150)
 REAL*8	 PARth(150),sum_chemup,sumXX(150),NG_dwn,alfNG_L1,lat,m_up,m_dwn,delta
 integer Nx,Nalf,M,chemup,L,I,KEYch(100)
 integer NK,NG,Nmol,Nbonds,chem(150),type_mol(150),rank_seg(150)
 REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),teta_chem(150,150)
 REAL*8 ZZ,ZR,JX,IJ
 COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
 COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
 common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
	 NG_dwn=0.0d0
	  do I=1,NG-1
	  NG_dwn=NG_dwn+Al_dwn(I,L+1)
	  enddo
	 
	sum_chemup=0.0d0
  	do M=1,Nbonds
	KEYch(M)=0
	sum_chemup=sum_chemup+2.0d0*PARth(Nx+Nalf/2+M)+sumXX(M)
	enddo
	alfNG_L1=m_dwn(L+1)/lat(L+1)*(sum_chemup+NG_dwn-1.0d0)+ALloc(NG,L+1)/lat(L+1)-m_up(L+1)/lat(L+1)*AL_up(NG,L+1)
	sum_chemup=sum_chemup+NG_dwn

	if ((chemup.eq.(-1)).and.(alfNG_L1.le.0.0d0)) then !sumup>1, sup=-1,but chemup can be changed
	chemup=0
	delta=0.0d0
	M=1
		do while ( ((dabs(alfNG_L1)/2.0d0).ge.delta).and.(M.le.Nbonds) )
		delta=delta+0.1d0*PARth(Nx+Nalf/2+M)
		KEYch(M)=1
		M=M+1
		enddo
		if (M.le.Nbonds) then
		do I=M,Nbonds
		KEYch(I)=-1
		enddo
		endif
	elseif (alfNG_L1.le.0.0d0) then
	chemup=1
	else
	chemup=0
	endif
	if (sum_chemup.ge.1.0d0) chemup=-1

		if (chemup.eq.(-1)) then
	do M=1,Nbonds
		if (PARth(Nx+Nalf/2+M).gt.1.0d-2) then
		chemup=0
		KEYch(M)=-1
		else
		KEYch(M)=0
		endif
	enddo
	endif

	if (chemup.eq.1) then
!	do M=1,Nbonds
!		if (teta_chem(M,L+1).ge.teta_chem(M,L)) then
!		chemup=0
!		KEYch(M)=1
!		else
!		KEYch(M)=0
!		endif
!	enddo
	endif
	!print *,alfNG_L1,sum_chemup
	if ((sum_chemup.ge.1.0d0).and.(alfNG_L1.le.0.0d0)) then
	print *,'different borders for chem in sum_of_chemup'
	 stop
	endif

  return
  end
  !******************************************************************************************************************
  SUBROUTINE minmaxALF(PARth,Nx,Nalf,Nsolve,ALlocnew,sumXX,L,KEYid,KEYch,sup,chemup,sdwn,allK,switch)
 IMPLICIT NONE
 character*7 name_seg(150)
 REAL*8	 PARth(150),sumup,sumdwn,sumal
 integer Nfix,No_or(150),CHANGE(150),s,Nx,Nalf,M,D,L,KEY(150),KEYid(100),KEYch(100),sup,chemup,sdwn,allK,switch,KEYngup,KEYngdwn
 integer NK,NG,Nmol,Nsolve,Nbonds,chem(150),type_mol(150),rank_seg(150),KEYs(150),KEYall
REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150) 
REAL*8 ZZ,ZR,JX,IJ,m_dwn,m_up,lat,alf_max(150),alf_min(150)
REAL*8 sumXX(150),sumX(150),ALlocnew(150),ALlocsol(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150)
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
common /min_maxALF/alf_max,alf_min
common /without_or/ Nfix,No_or,CHANGE
common /key_for_minmax/KEYs,KEYall,KEYngup,KEYngdwn
 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
 !switch=1 in the start of NEWTalf, =0 corrections all variables, =-D correction variable number D
  
	if (switch /= 0) then			 !the first time in border !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do D=Nsolve/2+1,Nsolve/2+Nbonds
	  sumX(D)=sumXX(D-Nsolve/2)
		enddo 

	  DO D=1,Nsolve/2+Nbonds
	  if (D.le.Nsolve/2) ALlocsol(D)=ALloc(CHANGE(D)-Nx,L)
	  if (D.gt.Nsolve/2) ALlocsol(D)=ALlocnew(D-Nsolve/2+NG)
	  
	  if (D.le.Nsolve/2) then
	  alf_max(D)=min((Allocsol(D)/m_up(L)),1.0d0)
	  alf_min(D)=max(1.0d-100,(Allocsol(D)/m_up(L)-(1.0d0-m_up(L))/m_up(L)))
	if ((L.ne.1).or.(r_solid.gt.1.0d-50))    alf_max(D+Nsolve/2+Nbonds)=min((Allocsol(D)/m_dwn(L)),1.0d0)
	if ((L.ne.1).or.(r_solid.gt.1.0d-50))    alf_min(D+Nsolve/2+Nbonds)=max(0.0d0,(Allocsol(D)/m_dwn(L)-(1.0d0-m_dwn(L))/m_dwn(L)))
	   if ((alf_min(D+Nsolve/2+Nbonds).ge.alf_max(D+Nsolve/2+Nbonds)).and.((L.ne.1).or.(r_solid.gt.1.0d-50))) then
		 	print *,'alf_min>alf_max in NEWTalf',alf_min(D+Nsolve/2+Nbonds),alf_max(D+Nsolve/2+Nbonds)
			stop
		endif
	  elseif (Nmol /= 0) then																			!*chem_up(D-Nsolve/2,L+1)
	  alf_min(D)=max(1.0d-100,-sumX(D),(teta_chem(D-Nsolve/2,L+1)/ZZ/m_dwn(L+1)-m_up(L+1)/m_dwn(L+1)*chem_up(D-Nsolve/2,L+1)-sumX(D)-lat(L+1)/m_dwn(L+1)*0.5d0),(Allocsol(D)/m_up(L)-m_dwn(L)*chem_dwn(D-Nsolve/2,L)/m_up(L)-lat(L)/m_up(L)*0.5d0))
	  alf_max(D)=min(1.0d0,(1.0d0-sumX(D)),(Allocsol(D)/m_up(L)-m_dwn(L)*chem_dwn(D-Nsolve/2,L)/m_up(L)),(teta_chem(D-Nsolve/2,L+1)/ZZ/m_dwn(L+1)-sumX(D)-m_up(L+1)/m_dwn(L+1)*chem_up(D-Nsolve/2,L+1) ))	  	 
	  endif																																							  !m_up(L+1)/m_dwn(L+1)*chem_up(D-Nsolve/2,L+1)
		if (alf_min(D).ge.alf_max(D)) then
		print *,'alf_min>alf_max in NEWTalf',D,alf_min(D),alf_max(D),(1.0d0-sumX(D)),(Allocsol(D)/m_up(L)-m_dwn(L)*chem_dwn(D-Nsolve/2,L)/m_up(L)),(teta_chem(D-Nsolve/2,L+1)/ZZ/m_dwn(L+1)-m_up(L+1)/m_dwn(L+1)*chem_up(D-Nsolve/2,L+1)-sumX(D) )	  	 
		print *,-sumX(D),(teta_chem(D-Nsolve/2,L+1)/ZZ/m_dwn(L+1)-m_up(L+1)/m_dwn(L+1)*chem_up(D-Nsolve/2,L+1)-sumX(D)-lat(L+1)/m_dwn(L+1)*0.5d0),(Allocsol(D)/m_up(L)-m_dwn(L)*chem_dwn(D-Nsolve/2,L)/m_up(L)-lat(L)/m_up(L)*0.5d0)
		stop
		endif		
	  enddo
	  endif !switch=1	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 
		do s=1,Nalf/2
	   	 Al_up(s,L)=PARth(Nx+s)				
		 Al_dwn(s,L)=PARth(Nx+Nalf/2+Nbonds+s)
		 enddo	
		 if (Nmol /= 0) then
		 do s=Nalf/2+1,Nalf/2+Nbonds
	   	 chem_up(s-Nalf/2,L)=PARth(Nx+s)
		 enddo 
		 endif
	!	print *,AL_up(1,L),L,switch,alf_min(1)
	  if (switch.ge.0) then	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  do D=1,Nsolve+Nbonds
	  	! print *,0.1d0*(alf_max(D)-alf_min(D)),D

	 KEYid(D)=0
	 	 if (((L.ne.1).or.(r_solid.gt.1.0d-50)).or.(D.le.(Nsolve/2+Nbonds))) then
	 if (PARth(CHANGE(D)).ge.alf_max(D)) KEYid(D)=-1
	 if (PARth(CHANGE(D)).le.alf_min(D)) KEYid(D)=1
	 endif
	 enddo 

	 allK=0
	 sup=0
	 chemup=0
	 sdwn=0

   call check_dependALF(PARth,L,KEY)	 !check chem_dwn	KEY(M)=1 to increase M

	 DO M=1,Nbonds
	 if ( (KEYid(Nsolve/2+M).eq.(-1)).and.(KEY(M).eq.1) ) then
	 print *,'different borders for chem',M,alf_max(Nsolve/2+M)
	 stop
	 elseif (KEY(M).eq.1) then
	 KEYid(Nsolve/2+M)=1
	 endif
	enddo

!!	 call  calc_sum_NG(L,AL_up,AL_dwn,ALloc,teta_chem,chem_up,chem_dwn,sumup,sumdwn,sumal,0,PARth)
!!	if (sumup.ge.1.0d0) sup=-1
!	if (sumup.ge.1.0d0) chemup=-1			!sumup>1, sup=-1,but chemup can be changed in sum_of_chemup
!!	if  ((L.ne.1).and.(sumdwn.ge.1.0d0)) sdwn=-1 	
if (Nmol /= 0) 	 call sum_of_chemup(PARth,Nx,Nalf,sumXX,L,chemup,KEYch)

  KEYall=0			!calc in calcAL_al_and_NG(L)
  DO  s=1,NG-1
  KEYs(s)=0
  enddo

	  call calcAL_al_and_NG(L)			  !KEYs(s)=-1 to decrease both s,KEYall=1 to increase all
!	  print *,'keys',keys(1),Al_al(1,L)
	  if  (((L.ne.1).or.(r_solid.gt.1.0d-50)).and.(KEYngdwn.eq.(-1)))    sdwn=-1 
	!  print *,sdwn
	  	  if ((KEYngup.eq.(-1)).and.(Nmol /= 0))	chemup=-1	
	if (KEYngup.eq.(-1))  sup=-1
   do D=1,Nsolve/2
   if ( (KEYid(D).ne.1).and.(KEYs(change(D)-Nx).eq.(-1)) ) KEYid(D)=-1
   if ( ((L.ne.1).or.(r_solid.gt.1.0d-50)).and.(KEYid(D+Nbonds+Nsolve/2).ne.1).and.(KEYs(change(D)-Nx).eq.(-1)) )  KEYid(D+Nbonds+Nsolve/2)=-1
		if (((KEYid(D).eq.1).and.(KEYs(change(D)-Nx).eq.(-1))).or.(((L.ne.1).or.(r_solid.gt.1.0d-50)).and.(KEYid(D+Nbonds+Nsolve/2).eq.1).and.(KEYs(change(D)-Nx).eq.(-1))))	then
		print *,'different borders for alfa nonchem',D
	!	stop
		endif
   enddo
   if (KEYall.eq.1) allK=1

   else		 !switch<0       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   D=-switch
   KEYid(D)=0
		if ( (D.gt.Nsolve/2).and.(D.lt.(Nsolve/2+Nbonds+1)) )	then
		 call check_dependALF(PARth,L,KEY)
		chemup=0
		call sum_of_chemup(PARth,Nx,Nalf,sumXX,L,chemup,KEYch)
		if ((KEY(D-Nsolve/2).eq.1).or.(chemup.eq.1).or.(KEYch(D-Nsolve/2).eq.1)) KEYid(D)=1
		if ((chemup.eq.(-1)).or.(KEYch(D-Nsolve/2).eq.(-1)))	 KEYid(D)=-1
		endif
  KEYngup=0
  KEYngdwn=0
  KEYall=0			!calc in calcAL_al_and_NG(L)
  KEYs(change(D)-Nx)=0
	 call calcAL_al_and_NG(L)
	 if ((D.le.Nsolve/2).and.(KEYs(change(D)-Nx).eq.(-1))) KEYid(D)=-1
	 if ((D.gt.(Nsolve/2+Nbonds)).and.(KEYs(change(D-Nsolve/2-Nbonds)-Nx).eq.(-1))) KEYid(D)=-1
	 if ((KEYngup.eq.(-1)).or.(KEYngdwn.eq.(-1)))  KEYid(D)=-1
	 if (KEYall.eq.1) KEYid(D)=1
	! print *,KEYs(change(D)-Nx)
	 if (PARth(CHANGE(D)).ge.alf_max(D)) KEYid(D)=-1
	 if (PARth(CHANGE(D)).le.alf_min(D)) KEYid(D)=1
   endif !switch
  return
  end
  !******************************************************************************************************************
   SUBROUTINE minmaxALF_fromTL(PARth,Nx,Nalf,Nsolve,L,switch)
 IMPLICIT NONE
 CHARACTER*3 shape
  character*7 name_seg(150)
 REAL*8	 PARth(150),sumup,sumdwn,sumal,TL(150),TLmin,TLold(150),T,Rout
 integer Lmax,Nfix,No_or(150),CHANGE(150),s,Nx,Nalf,M,D,L,KEY(150),KEYid(100),KEYch(100),sup,chemup,sdwn,allK,switch,KEYngup,KEYngdwn
 integer NK,NG,Nmol,Nsolve,Nbonds,chem(150),type_mol(150),rank_seg(150),KEYs(150),KEYall
REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150) 
REAL*8 ZZ,ZR,JX,IJ,m_dwn,m_up,lat,alf_max(150),alf_min(150)
REAL*8 sumXX(150),sumX(150),ALlocnew(150),ALlocsol(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150),minn,QT,RT
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
common /min_maxALF/alf_max,alf_min
common /without_or/ Nfix,No_or,CHANGE
common /key_for_minmax/KEYs,KEYall,KEYngup,KEYngdwn
	  common /layer/shape,TL,Lmax,TLmin
  common /tau_old/TLold
   switch=0
    if ((TL(L)-TLmin).lt.(1.0d-150)) 	then
	switch=1
	else
   if ((lat(L).le.1.0d-20).or.(lat(L).gt.1.0d0)) then
	 switch=1
	 else
	  KEYall=0			!calc in calcAL_al_and_NG(L)
	DO  s=1,NG-1
	KEYs(s)=0
	enddo
	 call calcAL_al_and_NG(L)
	 if (KEYall.ne.0) 		 switch=1
	  DO  s=1,NG-1
	if (KEYs(s).ne.0) 		 switch=1
	enddo
	!print *,KEYall,KEYs(1:NG-1)
	 endif
	  ! if ((shape.eq.'pln').or.(L.ne.1)) then
	  ! if (TL(L).gt.TLold(L)*1.5d0) 	switch=1
	   !else
	    if (TL(L).gt.TLold(L)*4.0d0) 	switch=1
	  ! endif

	endif
  
    return
  end
  !******************************************************************************************************************
  SUBROUTINE check_borderX(PAR,L,Lmax,Nx,id,mainK,switch)
  IMPLICIT NONE
  character*7 name_seg(150)
REAL*8 PAR(150)
integer Nx,L,I,D,kk,switch,KEYx(100),Lmax,allK,KEYxend(100),KEYxL(100),KEYid(100),id(100),mainK
integer NK,NG,Nmol,Nbonds,chem(150),type_mol(150),rank_seg(150)
REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150) 
REAL*8 ZZ,ZR,JX,IJ
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150)
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al	
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
 !switch=1 - went from initial profile, =0 - corrections   PAR, =-1 - checking   PAR after interation in NEWTx
 if (switch.ge.0) then
call minmaxX(PAR,L,Lmax,Nx,allK,KEYx,KEYxL,KEYxend,KEYid,1)
call id_to_changeX(Nx,allK,KEYx,KEYxL,KEYxend,KEYid,id,mainK,'before iter in NEWTx')  
	 do while (mainK.ne.0)
	 kk=kk+1
	 DO I=1,Nx
	 if (PAR(I).le.1.0d-150) then
	 PAR(I)=1.0d-5
	 elseif (PAR(I).ge.(1.0d0-1.0d-150)) then
	 PAR(I)=1.0d0-1.0d-5
	 else
		if (id(I).eq.1) then
		PAR(I)=PAR(I)*1.001d0
		elseif (id(I).eq.(-1)) then
		PAR(I)=PAR(I)*0.999d0 	
		endif
	endif
	 enddo !I
		    if (kk.gt.100000) then
			print *,'many steps for initial Xloc in NEWTx',PAR(1:Nx)
			print *,'allK=',allK
			write (*,'(6(6X,A7))')'1-Nx','id','KEYid','KEYx','KEYxL','KEYxend'
			do D=1,Nx
			print *,D,id(D),KEYid(D),KEYx(D),KEYxL(D),KEYxend(D)
			enddo	
			stop
			endif
call minmaxX(PAR,L,Lmax,Nx,allK,KEYx,KEYxL,KEYxend,KEYid,0)
	 call id_to_changeX(Nx,allK,KEYx,KEYxL,KEYxend,KEYid,id,mainK,'before iter in NEWTx')
	 enddo !while
		if (switch.gt.0) then
		Xloc(NK,L)=0.0d0
		DO I=1,Nx
		Xloc(I,L)=PAR(I)
		Xloc(NK,L)=Xloc(NK,L)+PAR(I)
		enddo
		Xloc(NK,L)=1.0d0-Xloc(NK,L)
		endif
 elseif (switch.lt.0) then
 call minmaxX(PAR,L,Lmax,Nx,allK,KEYx,KEYxL,KEYxend,KEYid,0)
 call id_to_changeX(Nx,allK,KEYx,KEYxL,KEYxend,KEYid,id,mainK,'after iter in NEWTx ')
 endif
return
  end
!******************************************************************************************************************
  SUBROUTINE minmaxX(PARx,L,Lmax,Nx,allK,KEYx,KEYxL,KEYxend,KEYid,switch)
  IMPLICIT NONE
  character*7 name_seg(150)
REAL*8 PAR(150),PARx(150),XlocL(150),min_sumX(150),max_sumX(150),min_sumR_1,max_sumR_1,NG_dwn,sumdwn,sumchemup,sumXloc,sumR_1,chem_min(150),chem_max(150)
integer Nx,L,LL,I,M,D,K,switch,seg_ch,KEYx(100),Lmax,allK,KEYxend(100),KEYxL(100),KEYid(100), rank_branch(150), n_branch, type_branch(150)
integer NK,NG,Nmol,Nbonds,chem(150),type_mol(150),rank_seg(150), seg_max_ch
REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150) 
REAL*8 ZZ,ZR,JX,IJ,m_dwn,m_up,lat,alfup_max(150),alfdwn_max(150),alfup_min(150),alfdwn_min(150)
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150),sumX(150),sumXX(150),ALlocnew(150),ALlocsol(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al	
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
common /min_maxX/max_sumX,min_sumX,max_sumR_1,min_sumR_1
common /branch/rank_branch, n_branch, type_branch
 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
!switch=1 in the start of NEWTx, =0 corrections all variables, =-I checking variable number I

if (switch /= 0) then
 NG_dwn=0.0d0
	  do I=1,NG-1
	  NG_dwn=NG_dwn+Al_dwn(I,L+1)
	  enddo
 sumdwn=0.0d0
 sumchemup=0.0d0
do M=1,Nbonds
if ((L.ne.1).or.(r_solid.gt.1.0d-50))	sumdwn=sumdwn+chem_dwn(M,L)+chem_up(M,L-1)
sumchemup=sumchemup+chem_up(M,L)
max_sumX(M)=min(1.0d0,(teta_chem(M,L+1)/ZZ/m_dwn(L+1)))	  !1.0d0-chem_up(M,L)	  !m_up(L+1)/m_dwn(L+1)*chem_up(M,L+1)
min_sumX(M)=max(-1.0d0,(teta_chem(M,L+1)/ZZ/m_dwn(L+1)-m_up(L+1)/m_dwn(L+1)-1.0d0-lat(L+1)/m_dwn(L+1)*0.5d0))	 !-chem_up(M,L)	!-m_up(L+1)/m_dwn(L+1)*chem_up(M,L+1)
	if (min_sumX(M).ge.max_sumX(M))	then
	print *,'min_sumX>max_sumX in NEWTx',M,min_sumX(M),max_sumX(M)
	stop
	endif
enddo
max_sumR_1=min((1.0d0-NG_dwn),(1.0d0-m_dwn(L)/m_up(L)*sumdwn),(lat(L+1)/m_dwn(L+1)-ALloc(NG,L+1)/m_dwn(L+1)+m_up(L+1)/m_dwn(L+1)*AL_up(NG,L+1)+1.0d0-NG_dwn))
min_sumR_1=max((-NG_dwn),(-1.0d0/m_up(L)-m_dwn(L)/m_up(L)*sumdwn),(-ALloc(NG,L+1)/m_dwn(L+1)+m_up(L+1)/m_dwn(L+1)*AL_up(NG,L+1)) )
	if (min_sumR_1.ge.max_sumR_1)	then
	print *,'min_sumR_1>max_sumR_1 in NEWTx',M,min_sumR_1,max_sumR_1
	stop
	endif
endif !switch

allK=0
sumXloc=0.0d0
Xloc(NK,L)=0.0d0
	  DO  I=1,Nx                                                       
      PAR(I)=PARx(I)
	  Xloc(I,L)=PARx(I)
	  Xloc(NK,L)=Xloc(NK,L)+PARx(I)

	  sumXloc=sumXloc+PARx(I)
	  KEYxend(I)=0
	  KEYxL(I)=0
	  KEYx(I)=0
	  KEYid(I)=0
	  if (PAR(I).ge.1.0d0) KEYid(I)=-1
	  if (PAR(I).le.1.0d-100) KEYid(I)=1
	  enddo
Xloc(NK,L)=1.0d0-Xloc(NK,L)
if (sumXloc.ge.1.0d0) allK=-1
 	
  call calcALloc(PAR,L,ALlocnew,sumXX)

	 DO D=1,NG-1
	 alfup_max(D)=min((Allocnew(D)/m_up(L)),1.0d0)
	  alfup_min(D)=max(0.0d0,(Allocnew(D)/m_up(L)-(1.0d0-m_up(L))/m_up(L)))
	if ((L.ne.1).or.(r_solid.gt.1.0d-50))  alfdwn_max(D)=min((Allocnew(D)/m_dwn(L)),1.0d0)
	if ((L.ne.1).or.(r_solid.gt.1.0d-50))  alfdwn_min(D)=max(0.0d0,(Allocnew(D)/m_dwn(L)-(1.0d0-m_dwn(L))/m_dwn(L)))
		if ((alfup_min(D).ge.1.0d0).or. (((L.ne.1).or.(r_solid.gt.1.0d-50)).and.(alfdwn_min(D).ge.1.0d0)) ) then
			do k=1,Nx
			if ((QS(k,D).gt.0.0d0).and.(KEYid(k).ne.1)) 	KEYid(k)=-1
			enddo
		 elseif	((alfup_max(D).le.0.0d0).or.(((L.ne.1).or.(r_solid.gt.1.0d-50)).and.(alfdwn_max(D).le.0.0d0)) ) then
			do k=1,Nx
			if ((QS(k,D).gt.0.0d0).and.(KEYid(k).ne.(-1))) 	KEYid(k)=1
			enddo
		 endif
	  enddo

   do D=1,Nbonds
   ALlocsol(D)=ALlocnew(D+NG)							   !-m_up(L+1)/m_dwn(L+1)*chem_up(D,L+1)
   sumX(D)=sumXX(D)
   chem_min(D)=max(-sumX(D),(teta_chem(D,L+1)/ZZ/m_dwn(L+1)-m_up(L+1)/m_dwn(L+1)*chem_up(D,L+1)-sumX(D)-lat(L+1)/m_dwn(L+1)*0.5d0),(Allocsol(D)/m_up(L)-m_dwn(L)*chem_dwn(D,L)/m_up(L)-lat(L)/m_up(L)*0.5d0))
   chem_max(D)=min((1.0d0-sumX(D)),(Allocsol(D)/m_up(L)-m_dwn(L)*chem_dwn(D,L)/m_up(L)),(teta_chem(D,L+1)/ZZ/m_dwn(L+1)-sumX(D)-m_up(L+1)/m_dwn(L+1)*chem_up(D,L+1) ))	  	 
	if ((chem_min(D).ge.1.0d0).and.(KEYid(seg_ch(D)).ne.1)) KEYid(seg_ch(D))=-1		!.or.(chem_min(D).ge.chem_max(D))							   !-m_up(L+1)/m_dwn(L+1)*chem_up(D,L+1)
	!	if ((chem_min(D).ge.1.0d0)) KEYid(seg_ch(D)+1)=1  !???
	if ((chem_max(D).le.0.0d0).and.(KEYid(seg_ch(D)).ne.(-1))) KEYid(seg_ch(D))=1
	if ((chem_max(D).le.0.0d0).and.(KEYid(seg_max_ch(D)).ne.1)) KEYid(seg_max_ch(D))=-1 
	!if ((switch.lt.0).and.(chem_max(D).le.0.0d0)) KEYid(seg_ch(D)+1)=-1
	if 	(((L.ne.1).or.(r_solid.gt.1.0d-50)).and.(chem_min(D).ge.chem_max(D))) then
		  if ( (Xloc(seg_max_ch(D),L)/Xloc(seg_ch(D),L)*teta_chem(D,L)/ZZ-teta_chem(D,L+1)/ZZ).ge.(m_up(L)*chem_up(D,L-1)+lat(L+1)*0.5d0) ) then
		if (KEYid(seg_max_ch(D)).ne.1)  KEYid(seg_max_ch(D))=-1
		  elseif ((Xloc(seg_max_ch(D),L)/Xloc(seg_ch(D),L)*teta_chem(D,L)/ZZ-teta_chem(D,L+1)/ZZ).le.(m_up(L)*chem_up(D,L-1)-m_up(L)*chem_up(D,L+1)-lat(L+1)*0.5d0) ) then
			if (KEYid(seg_max_ch(D)).ne.(-1))   KEYid(seg_max_ch(D))=1
		  elseif ( (Xloc(seg_max_ch(D),L)/Xloc(seg_ch(D),L)*teta_chem(D,L)/ZZ).le.(m_up(L)*chem_up(D,L-1)) ) then
			if (KEYid(seg_max_ch(D)).ne.(-1))   KEYid(seg_max_ch(D))=1
		  endif
	endif
	enddo	  
 !if (L.eq.4)  print *,(teta_chem(2,L+1)/ZZ/m_dwn(L+1)-sumX(2)-m_up(L+1)/m_dwn(L+1)*chem_up(2,L+1) ),KEYid(seg_ch(2)),switch

if ((switch.lt.0).and.(allK.ne.0)) KEYid(-switch)=allK
if ((switch.lt.0).and.(KEYid(-switch).ne.0)) return

	  do LL=L+1,Lmax-1
	  DO  I=1,Nx                                                       
      XlocL(I)=Xloc(I,LL)
	  enddo
	  call  calcALchem(XlocL,LL,Xchem,sumX)			   !-chem_up<sumX<1-chem_up
		do M=1,Nbonds
		if (sumX(M).ge.(1.0d0-chem_up(M,LL))) then			  !1.0d0-chem_up(M,LL)
			 if (KEYx(seg_ch(M)).ne.(-1)) KEYx(seg_ch(M))=1
			 if (KEYx(seg_max_ch(M)).ne.1) KEYx(seg_max_ch(M))=-1
	!		 print *,'LL=',LL,M
		elseif (sumX(M).le.(-chem_up(M,LL))) then			   !-chem_up(M,LL)
			if (KEYx(seg_ch(M)).ne.1) KEYx(seg_ch(M))=-1
			if (KEYx(seg_max_ch(M)).ne.(-1)) KEYx(seg_max_ch(M))=1
	!		print *,'LL=',LL,M
		endif
		enddo

	  enddo
	
 if ((switch.lt.0).and.(KEYx(-switch).ne.0))  then
 KEYid(-switch)=KEYx(-switch)
 return
 endif
	 
sumR_1=0.0d0
do M=1,Nbonds
sumR_1=sumR_1+sumXX(M)
if (sumXX(M).le.min_sumX(M)) then
KEYxL(seg_ch(M))=-1
KEYxL(seg_max_ch(M))=1
elseif (sumXX(M).ge.max_sumX(M)) then
KEYxL(seg_ch(M))=1
KEYxL(seg_max_ch(M))=-1
endif 
enddo !M
if (sumR_1.le.min_sumR_1) then
	do I=1,Nx
	if (rank_seg(I).eq.1) KEYxend(I)=-1
	if ((rank_seg(I).ne.0).and.(rank_seg(I+1).le.1)) KEYxend(I)=1
	if (rank_branch(I) .eq. 1) KEYxend(I)=-1
	if ((rank_branch(I) .gt. 1) .and. (type_branch(I+1) .ne. type_branch(I))) KEYxend(I)=1 
	enddo  !I
elseif (sumR_1.ge.max_sumR_1) then
	do I=1,Nx
	if (rank_seg(I).eq.1) KEYxend(I)=1
	if ((rank_seg(I).ne.0).and.(rank_seg(I+1).le.1)) KEYxend(I)=-1
	if (rank_branch(I) .eq. 1) KEYxend(I)=1
	if ((rank_branch(I) .gt. 1) .and. (type_branch(I+1) .ne. type_branch(I))) KEYxend(I)=-1
	enddo  !I
endif  !sumR_1

  if ((switch.lt.0).and.(KEYxL(-switch).ne.0))  then
 KEYid(-switch)=KEYxL(-switch)
 elseif ((switch.lt.0).and.(KEYxend(-switch).ne.0))  then
 KEYid(-switch)=KEYxend(-switch)
 endif

  return
  end
  !******************************************************************************************************************
  	SUBROUTINE	id_to_changeX(Nx,allK,KEYx,KEYxL,KEYxend,KEYid,id,mainK,switch)
	!output: id(D)=1 to increase D, =-1 to decrease D, =0 not to change D
	!allK=-1 (decrease all) or =0 (not matter)
	!KEYid(I)=1 (Xloc<min,increase Xloc(I)), =-1 (Xloc>max,decrease) or 0
	!KEYx(I) 1,-1,0 (Xloc(I))		!for two segments bonded chemically
	!KEYxL(I) 1,-1,0 (Xloc(I))		!for two segments bonded chemically
	!KEYxend(I) 1,-1,0 (Xloc(I))	!for the first and last segments of chains
	!mainK=0 when all id(D)=0, otherwise 1
	IMPLICIT NONE
	CHARACTER*20 switch
	integer	 D,Nx,allK,KEYx(100),KEYxL(100),KEYxend(100),KEYid(100),id(100),n,mainK
	n=0
	mainK=0

	do D=1,Nx
	id(D)=0
	  if (KEYid(D).ne.0) then	 
	  id(D)=KEYid(D)
	  else	! KEYid(D)=0
	!	if (KEYxend(D).ne.0) then
	!	id(D)=KEYxend(D)
	!	else   !KEYxend(D)=0
			if ((KEYx(D)*KEYxL(D)).lt.0) then
			print *,'different borders Xloc ',switch,D
		!	stop
			id(D)=KEYxL(D)
		!	id(D)=0
		!	mainK=1
			elseif (KEYx(D).ne.0) then
			id(D)=KEYx(D)
			elseif (KEYxL(D).ne.0) then	 !KEYx(D)=0
			id(D)=KEYxL(D)
			elseif (KEYxend(D).ne.0) then !KEYxL(D)=0
			id(D)=KEYxend(D)
			elseif (allK.ne.0) then 	 
			id(D)=allK
			endif
	!	endif
	  endif
	  if (id(D).ne.0)  mainK=1
	  if ((allK.ne.0).and.(allK.eq.id(D)) ) n=allK
	  enddo !D

	if 	 (allK.ne.n) then
	print *,'different borders Xloc ',switch
	print *,allk,n
	stop
	endif
   return
  end
  !******************************************************************************************************************
SUBROUTINE calc_sum_NG(L,AL_up,AL_dwn,ALloc,teta_chem,chem_up,chem_dwn,sumup,sumdwn,sumal,switch,PARth)   
   IMPLICIT NONE
   real*8  sumdwn,sumup,sumal,m_up,m_dwn,lat,ZZ,ZR,JX,IJ,DPx(150),DPalf(150),EPS,EPF
   REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),PARth(150)
   REAL*8 AL_up(10,150),AL_dwn(10,150),ALloc(10,150) 
   integer k,s,L,NG,Nbonds,Nmol,NK,switch,Nx,Nalf
   COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
    common /forNEWT/Nx,Nalf,DPx,DPalf,EPS,EPF
    REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid

  IF (switch /= 0) then
		do s=1,Nalf/2
	   	 Al_up(s,L)=PARth(Nx+s)				
		 Al_dwn(s,L)=PARth(Nx+Nalf/2+Nbonds+s)
		 enddo	
		 do s=Nalf/2+1,Nalf/2+Nbonds
	   	 chem_up(s-Nalf/2,L)=PARth(Nx+s)
		 enddo
  endif
sumdwn=0.0d0
sumup=0.0d0
sumal=0.0d0

	do s=1,NG-1
	sumdwn=sumdwn+AL_dwn(s,L)
	sumup=sumup+AL_up(s,L)
	sumal=sumal+ALloc(s,L)-m_up(L)*AL_up(s,L)-m_dwn(L)*AL_dwn(s,L)
	enddo
	DO k=1,Nbonds
	if ((L.ne.1).or.(r_solid.gt.1.0d-50))	sumdwn=sumdwn+chem_dwn(k,L)+chem_up(k,L-1)
	sumup=sumup+chem_up(k,L)+chem_dwn(k,L+1)
	sumal=sumal+2.0d0*(teta_chem(k,L)/ZZ-m_up(L)*chem_up(k,L)-m_dwn(L)*chem_dwn(k,L))
	enddo
	sumal=sumal/lat(L)
  return
  end
  !******************************************************************************************************************
  SUBROUTINE correct_ALF_forNG(L,AL_up,AL_dwn,ALloc,teta_chem,chem_up,chem_dwn,sumup,sumdwn,sumal,switch,PARth)   
   IMPLICIT NONE
   real*8  sumdwn,sumup,sumal,m_up,m_dwn,ZZ,ZR,JX,IJ,DPx(150),DPalf(150),EPS,EPF,alfup_max(150), alfup_min(150),alfdwn_max(150),alfdwn_min(150)
   REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),PARth(150)
   REAL*8 AL_up(10,150),AL_dwn(10,150),ALloc(10,150) 
   integer k,s,L,NG,Nbonds,Nmol,NK,switch,Nx,Nalf
   COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
    common /forNEWT/Nx,Nalf,DPx,DPalf,EPS,EPF
	 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid

	do s=1,NG-1
	   alfup_max(s)=min((ALloc(s,L)/m_up(L)),1.0d0)
	  alfup_min(s)=max(1.0d-150,(ALloc(s,L)/m_up(L)-(1.0d0-m_up(L))/m_up(L)))
	  alfdwn_max(s)=1.0d0
	  alfdwn_min(s)=0.0d0
	if ((L.ne.1).or.(r_solid.gt.1.0d-50))  alfdwn_max(s)=min((ALloc(s,L)/m_dwn(L)),1.0d0)
	if ((L.ne.1).or.(r_solid.gt.1.0d-50))  alfdwn_min(s)=max(1.0d-150,(ALloc(s,L)/m_dwn(L)-(1.0d0-m_dwn(L))/m_dwn(L)))
	 enddo
  	call calc_sum_NG(L,AL_up,AL_dwn,ALloc,teta_chem,chem_up,chem_dwn,sumup,sumdwn,sumal,switch,PARth) 
		k=0
		do while ((sumal.ge.1.0d0).or.(sumdwn.ge.1.0d0).or.(sumup.ge.1.0d0))
		k=k+1
		if  ((sumal.ge.1.0d0).and.(sumdwn.ge.1.0d0).and.(sumup.ge.1.0d0)) then
		print *,'wrong initial sums'
		stop
		endif
	if (sumal.ge.1.0d0) then 
		do s=1,NG-1
		if ((sumup.lt.1.0d0).and.(AL_up(s,L).lt.alfup_max(s))) AL_up(s,L)=AL_up(s,L)*1.0001d0 ! +(sumal-1.0d0)/((NG-1.0d0)) !*1.001d0
		if ((sumdwn.lt.1.0d0).and.(AL_dwn(s,L).lt.alfdwn_max(s))) AL_dwn(s,L)=AL_dwn(s,L)*1.0001d0 !+(sumal-1.0d0)/((NG-1.0d0)) !*1.001d0
		
		if ((AL_up(s,L).ge.1.0d0).or.(AL_dwn(s,L).ge.1.0d0)) then
		print *,'initial AL large',AL_up(s,L),AL_dwn(s,L),s,l
		stop
		endif
		enddo
	elseif 	(sumdwn.ge.1.0d0) then
		do s=1,NG-1
		if (AL_dwn(s,L).gt.alfdwn_min(s)) AL_dwn(s,L)=AL_dwn(s,L)*0.9999d0 !-(sumdwn-1.0d0)/(NG-1.0d0) !*0.999d0
		if (AL_dwn(s,L).le.1.0d-150) then
		print *,'initial AL_dwn small',AL_dwn(s,L),s,l
		!stop
		AL_dwn(s,L)=1.0d-10
		elseif (AL_dwn(s,L).ge.(1.0d0-1.0d-150)) then
		AL_dwn(s,L)=1.0d0-1.0d-10
		endif
		enddo
	else
		do s=1,NG-1
	if (AL_up(s,L).gt.alfup_min(s))	AL_up(s,L)=AL_up(s,L)*0.9999d0 !-(sumup-1.0d0)/(NG-1.0d0) !*0.999d0
		if (AL_up(s,L).le.1.0d-150) then
		print *,'initial AL_up small',AL_up(s,L),s,l
		!stop
		AL_up(s,L)=1.0d-10
		elseif (AL_up(s,L).ge.(1.0d0-1.0d-150)) then
		AL_up(s,L)=1.0d0-1.0d-10
		endif
		enddo
	endif
	call calc_sum_NG(L,AL_up,AL_dwn,ALloc,teta_chem,chem_up,chem_dwn,sumup,sumdwn,sumal,switch,PARth) 
	enddo !while

  return
  end
  !******************************************************************************************************************
