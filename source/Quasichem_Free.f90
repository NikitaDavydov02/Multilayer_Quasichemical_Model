! Quasichem  Free Energy
!*************************************************************************************************************
SUBROUTINE PARAM(T,TT)                                                  
IMPLICIT NONE
REAL*8 ZZ,ZR,JX,IJ,W(10,10),WT(10,10),WTT(10,10),REDT(10,10)
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150),C1,C2,T,TT 
integer Nmol,NG,I,J,NK,Nbonds                                         
COMMON /B/ W,WT,WTT,REDT
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                      
COMMON /H/ BLIS,AL,ALI,ET,GP
 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid                                       
	C1=TT/T                                                           
    C2=C1-1.0d0                                                         
    C1=DLOG(C1)-C2                                                    
DO 5 I=1,NG  
ET_solid(I)=DEXP(-W_solid(I)-WT_solid(I)*C2-WTT_solid(I)*C1)                                                       
DO 5 J=1,NG                                                       
   5  ET(I,J)=DEXP(-W(I,J)-WT(I,J)*C2-WTT(I,J)*C1)                  
RETURN                                                            
END
!*************************************************************************************************************
real*8 FUNCTION pair(s,L,direct,t)
IMPLICIT none
REAL*8 m_up,m_dwn,lat
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
CHARACTER*2  direct
integer s,L,t
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150)                                     
COMMON /H/ BLIS,AL,ALI,ET,GP
COMMON /S_or/Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/AL_up,AL_dwn,AL_al,ALloc

if (direct.eq.'up') then
pair=m_up(L)*AL_up(s,L)*AL_dwn(t,L+1)*XX_up(s,L)*XX_dwn(t,L+1)*ET(s,t)
elseif (direct.eq.'dw')  then
	if (L.eq.1)	 pair=0.0d0
	if (L.ne.1)	 pair=m_dwn(L)*AL_dwn(s,L)*AL_up(t,L-1)*XX_dwn(s,L)*XX_up(t,L-1)*ET(s,t)		 !pairs between groups s in L and t in L-1
else
pair=lat(L)*AL_al(s,L)*AL_al(t,L)*XX_al(s,L)*XX_al(t,L)*ET(s,t)		             !pairs between groups s in L and t in L
if (s.eq.t) pair=pair/2.0d0
endif
!A_upb(s,t,L)=m_up(L)*AL_up(s,L)*AL_dwn(t,L+1)*ET(s,t)		                             !pairs based on Bragg-Williams
!IF (L.ne.1) A_dwnb(s,t,L)=m_dwn(L)*AL_dwn(s,L)*AL_up(t,L-1)*ET(s,t)		 
!A_alb(s,t,L)=lat(L)*AL_al(s,L)*AL_al(t,L)*ET(s,t)		             					 
!if (s.eq.t) A_alb(s,t,L)=A_alb(s,t,L)/2

 END FUNCTION
!*************************************************************************************************************
 SUBROUTINE	 sum_pairs(A_up,A_dwn,A_al,L,NG,A_updwnsum,A_alsum)		 !A_updwnsum - fraction of transverse pairs s-t in L, A_alsum - lateral pairs 
 IMPLICIT NONE
 REAL*8	 A_up(10,10,150),A_dwn(10,10,150),A_al(10,10,150),A_updwnsum(150,150,150),A_alsum(150,150,150),A_sum
 integer s,t,L,NG

DO s=1,NG
 A_sum=0.0d0
 DO t=s,NG
 if (t.ne.s) A_updwnsum(s,t,L)=A_up(t,s,L)+A_up(s,t,L)+A_dwn(t,s,L)+A_dwn(s,t,L)
 if (t.eq.s) A_updwnsum(s,t,L)=A_up(s,t,L)+A_dwn(t,s,L)
 A_alsum(s,t,L)=A_al(t,s,L) !+A_al(s,t,L)
 A_sum=A_updwnsum(s,t,L)+A_alsum(s,t,L)
 A_updwnsum(s,t,L)= A_updwnsum(s,t,L)/A_sum
 A_alsum(s,t,L)= A_alsum(s,t,L)/A_sum
 enddo
 enddo
 return
 END
!*************************************************************************************************************
 SUBROUTINE	 sum_chempairs(chem_up,chem_dwn,chem_al,L,Lmax,Nbonds,Ach_updwnsum,Ach_alsum)	
  !Ach_updwnsum - fraction of transverse chem bonds M in L, Ach_alsum - lateral bonds 
 IMPLICIT NONE
 REAL*8	 A_sum,Ach_updwnsum(150,150),Ach_alsum(150,150),m_up,m_dwn,lat
 real*8   chem_up(150,150),chem_dwn(150,150),chem_al(150,150)
 integer M,L,Lmax,Nbonds

do M=1,Nbonds
 if (L.eq.1) then
 Ach_updwnsum(M,L)=m_up(L)*chem_up(M,L)+m_dwn(L)*chem_dwn(M,L)*2.0d0+m_dwn(L+1)*chem_dwn(M,L+1)
 elseif (L.eq.Lmax) then
 Ach_updwnsum(M,L)=m_up(L)*chem_up(M,L)*2.0d0+m_dwn(L)*chem_dwn(M,L)+m_up(L-1)*chem_up(M,L-1)
 else
 Ach_updwnsum(M,L)=m_up(L)*chem_up(M,L)+m_dwn(L)*chem_dwn(M,L)+m_dwn(L+1)*chem_dwn(M,L+1)+m_up(L-1)*chem_up(M,L-1)
 endif
 Ach_alsum(M,L)=lat(L)*chem_al(M,L) !/2.0d0
 A_sum=Ach_updwnsum(M,L)+Ach_alsum(M,L)
 Ach_updwnsum(M,L)=Ach_updwnsum(M,L)/A_sum 
 Ach_alsum(M,L)=Ach_alsum(M,L)/A_sum
 enddo
 return
 END
!*************************************************************************************************************
SUBROUTINE comb_part(R,Q,NK,RT,QT,ZR,GGA)			 	!combinatorial part of chem pot			 !,X
IMPLICIT NONE
REAL*8 A,BB,RT,QT,ZR,R(150),Q(150),GGA(150)	   !,X(150)
integer k,NK
do k=1,NK
 A=R(k)/RT                                                         
 BB=A*QT/Q(k)                                        
 GGA(k)=-ZR*Q(k)*(1.0d0-BB+DLOG(BB))+DLOG(A)+1.0d0-A		!+DLOG(X(k))
enddo    
RETURN                                                            
END
!*************************************************************************************************************
SUBROUTINE chemparts(L,chempart_segm,chempart_short)                                         ! contribution of bonds in chem pot        
IMPLICIT NONE
CHARACTER*3 shape
character*7 name_seg(150) 
real*8   chem_up(150,150),chem_dwn(150,150),chem_al(150,150),ALchemold(150,150),teta_chem(150,150),chempart_segm(150),chempart(150), sum, chem(150),chempart_short(150)
REAL*8 ZZ,ZR,JX,IJ,m_up,m_dwn,XlocL(150),TL(150),TLmin
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
integer   L,NK,NG,I,M,Nmol,Nbonds,k,Lmax,LL,chemm(150),type_mol(150),rank_seg(150), seg_max_ch, seg_ch 
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds  
common /layer/shape,TL,Lmax,TLmin
common /chem/chemm,type_mol,name_seg,teta_chem,rank_seg
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al
do k=1,NK
XlocL(k)=Xloc(k,L)
chem(k) = chemm(k)
enddo
call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
DO  I=1,NK
chempart_segm(I)= chem(I)*(m_up(L)*(  dlog(XX_up(NG,L)*AL_up(NG,L)/(XX_al(NG,L)*AL_al(NG,L)))) -dlog(Xloc(I,L)/QT/ZZ))
IF  (L.ne.1) then
	chempart_segm(I)=chempart_segm(I)+chem(I)*m_dwn(L)*(  dlog(XX_dwn(NG,L)*AL_dwn(NG,L)/(XX_al(NG,L)*AL_al(NG,L)))  ) 
endif
chempart_short(I) = chempart_segm(I)
enddo
do	k=1,Nbonds
sum = 0.0d0
	DO LL=L,Lmax
	sum = sum+0.5d0*dlog(chem_up(k,LL)/chem_dwn(k,LL+1))
	ENDDO
	   !if (chem_al(k,L).le.0.0d0) then 
		!print *,'chemal',k,L,chem_al(k,L)
		!endif
	chempart_short(seg_max_ch(k)) = chempart_short(seg_max_ch(k)) + 0.5d0*dlog(chem_al(k,L))
	chempart_short(seg_ch(k)) = chempart_short(seg_ch(k)) + 0.5d0*dlog(chem_al(k,L))

 chempart_segm(seg_max_ch(k)) = chempart_segm(seg_max_ch(k)) + 0.5d0*dlog(chem_al(k,L)) - sum
 chempart_segm(seg_ch(k)) = chempart_segm(seg_ch(k)) + 0.5d0*dlog(chem_al(k,L)) + sum
 enddo

RETURN                                                            
END
!*************************************************************************************************************
SUBROUTINE IKS_or(phase,L)                                         !L - number of layer, shape - type of lattice           
IMPLICIT NONE  
integer NK,NG,k,I,J1,Nmol,Nbonds
integer L ,phase,Lmax 
CHARACTER*3 shape                                       
REAL*8 RXX_al(10),RXX_up(10),TL(150),RE,BB_al,BB_up,BB_dwn
REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150),XX_uphom(150),TLmin
REAL*8 ZZ,ZR,JX,IJ
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150) 
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150),XXX_up(10),XXX_dwn(10),XXX_al(10)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                      
COMMON /H/ BLIS,AL,ALI,ET,GP         
common /layer/shape,TL,Lmax,TLmin
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
common /hom/XX_uphom
!phase 1 - L+1 reflects L; 0 - inhom; -1 - L=1 (or L-1 is homog)

! write(*,*)'start of IKS_or',L
! write(*,*)shape,phase                 
	RE=1.0d0                                                            
    J1=0                                                              
DO K=1,NG
    RE=XX(K)                    !homogeneous value as initial point                                                       
	XX_al(K,L)=RE !-1.0D-04 
	if (phase.eq.1) then
	XX_up(K,L)=RE !-1.0D-04   
	endif   
enddo                                                   
 !  1  RE=RE-1.0D-04                                                     
  19  J1=J1+1                                                                                                   
DO 3 K=1,NG                                                       
      BB_al=0.0d0
	  BB_dwn=0.0d0
	  BB_up=0.0d0 
	                                  
DO 4 I=1,NG                                                       
	if (phase.eq.1) then
	 BB_up=BB_up+AL_up(I,L)*XX_up(I,L)*ET(K,I)
	endif
   4 BB_al=BB_al+AL_al(I,L)*XX_al(I,L)*ET(K,I)
if (phase.eq.1) RXX_up(K)=1.0d0/BB_up                           
   3 RXX_al(K)=1.0d0/BB_al   
                           
DO 6 K=1,NG  
     IF (phase.eq.1) then
  XX_up(K,L)=(XX_up(K,L)+RXX_up(K))/2.0d0
	  endif                                      
   6  XX_al(K,L)=(XX_al(K,L)+RXX_al(K))/2.0d0                                         
 DO 7 K=1,NG  
	if (phase.eq.1) then
    IF (DABS(XX_up(K,L)-RXX_up(K)).GT.1.0D-13) GO TO 19
	 endif
	IF (DABS(XX_al(K,L)-RXX_al(K)).GT.1.0D-13) GO TO 19 
   7  CONTINUE 

if (J1.gt.1000) write(*,*)'Too many steps in IKS_or'
if (J1.gt.1000) stop  

 DO I=1,NG
 XXX_al(I)=XX_al(I,L)
 ENDDO

 call ITER_iks(phase,L,XXX_up,XXX_dwn,XXX_al)
 DO I=1,NG
 XX_up(I,L)=XXX_up(I)
 XX_dwn(I,L+1)=XXX_dwn(I)
 ENDDO

 if (phase.eq.1) then
DO I=1,NG
 XX_al(I,L+1)=XX_al(I,L)
 ENDDO
 endif
 if (phase.ne.(-1)) then
 do I=1,NG
 XXX_al(i)=XX_al(I,L-1)
 enddo
 if (phase.eq.1) call ITER_iks(phase-1,L-1,XXX_up,XXX_dwn,XXX_al)
 if (phase.ne.1) call ITER_iks(phase,L-1,XXX_up,XXX_dwn,XXX_al)
endif
 DO I=1,NG
 if (L-1 /= 0) XX_up(I,L-1)=XXX_up(I)
 !!if (L.eq.1) XX_uphom(I)=XX_dwn(I,L)
if (L-1 /= 0) XX_dwn(I,L)=XXX_dwn(I)
 ENDDO 
                                    
!write(*,*)'end of IKS_or'   
20  RETURN                                                            
END
!************************************************************************************************************
SUBROUTINE ITER_iks(phase,L,XX_up,XX_dwn,XX_al)			 !solving system of quasichem equations
IMPLICIT NONE 
integer NK,NG,I,Nmol,J,Nbonds
integer L ,phase,Lmax 
CHARACTER*3 shape                                       
REAL*8 TL(150),TLmin
REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150),XX_up(10),XX_dwn(10),XX_al(10)
REAL*8 ZZ,ZR,JX,IJ
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150) 
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                      
COMMON /H/ BLIS,AL,ALI,ET,GP         
common /layer/shape,TL,Lmax,TLmin
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
!phase 1 => L+1 reflects L; 0 => inhom; -1 => L=1 (or L-1 is homog)
integer KK
real*8 AAA(20,20),A
real*8 FFAAA(20),XX_or(20), XXR(20),AL_or(20)
!print *,'iterrr',L,phase
A=1.0d0
do i =1,2*NG
do j =1,2*NG
AAA(i,j)=0.0d0
enddo
FFAAA(i)=1.0d0
if (i.le.NG) XX_or(i) = XX_AL(i)
if (i.gt.NG) XX_or(i) = XX_AL(i-NG)
!A=A-1.0d-05
enddo
DO i =1,2*NG
if (i.le.NG) then
	if (phase.eq.1) then
	AL_or(i)=AL_up(i,L)							!L+1 reflects L
	else
	AL_or(i)=AL_dwn(i,L+1)						!i=1..NG       XX_up(L)*sum(XX_dwn(L+1)*AL_dwn(L+1)*ET)	 XX_up(L-1)*sum(XX_dwn(L)*AL_dwn(L)*ET)
	endif
endif
if (i.gt.NG) then
	if (L.eq.0) then
	AL_or(i)=AL_dwn(i-NG,L)							!0 is homog
	else 
	AL_or(i)=AL_up(i-NG,L)				    	!i=NG+1..2*NG  XX_dwn(L+1)*sum(XX_up(L)*AL_up(L)*ET)	 XX_dwn(L)*sum(XX_up(L-1)*AL_up(L-1)*ET)
	endif
endif
enddo

! print *,'kkk',AL_or
 do i =1,NG
 do j =1,NG
 aaa(i,j+NG) =AL_or(j)*ET(i,j)
 aaa(i+NG,j) =AL_or(j+NG)*ET(i,j) 
 enddo
 enddo

 kk=0  

77 kk=kk+1 
 do i=1,2*NG
 XXR(i)=XX_or(i)
 enddo
 	  
 do i =1, 2*NG
	  A=0.0d0
 do j = 1, 2*ng
  A=A+AAA(i,j)*XX_or(j)
 enddo

 XX_or(i)=FFAAA(i)/A
 enddo
!  print *,XX_or
 if (kk.gt.1000) then
  !print *,AL(1),AL(2),AL(3),L,L+1,AL(1)+AL(2)+AL(3)
 print *,AL_up(1,L),AL_up(2,L),AL_up(3,L),AL_up(4,L),L,L+1,AL_up(1,L)+AL_up(2,L)+AL_up(3,L) !+AL_up(4,L)
 print *,AL_dwn(1,L+1),AL_dwn(2,L+1),AL_dwn(3,L+1),AL_dwn(4,L+1),L,L+1,AL_dwn(1,L+1)+AL_dwn(2,L+1)+AL_dwn(3,L+1) !+AL_dwn(4,L+1)
 print *,'too many steps in XX_or'
 print *,XX_or
 stop
 endif
 do i=1,2*ng
 if(Dabs(XX_or(i)-XXR(i)).gt.1.0d-13) goto 77					!!!!!!!!!  -10
 enddo

 do i =1, 2*ng
	  A=0.0d0
 do j = 1, 2*ng
  A=A+AAA(i,j)*XX_or(j)
 enddo

 !print *, (A*XX_or(i)-FFAAA(i))
 enddo
! print *,XX_or	,'kk=   ',kk
 !stop
 DO i =1,2*NG
if (i.le.NG) XX_up(i)=XX_or(i)						!i=1..NG       XX_up(i,L)	     XX_up(i,L-1)
if (i.gt.NG) XX_dwn(i-NG)=XX_or(i)				    !i=NG+1..2*NG  XX_dwn(i,L+1)	 XX_dwn(i,L)
enddo
!print *,'kkk',XX_or
RETURN                                                            
END
!************************************************************************************************************
SUBROUTINE IKS                                            !quasichemical equations for homogeneous case           
IMPLICIT NONE 
integer NK,NG,k,I,J1,Nmol,Nbonds                                   
REAL*8 RXX(10),RE,BB
REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150)
REAL*8 ZZ,ZR,JX,IJ
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150) 
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                      
COMMON /H/ BLIS,AL,ALI,ET,GP                                       
	RE=1.0d0                                                            
    J1=0  
!	if (AL(4).eq.0.0d0) NG=NG-1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                            
DO 1 K=1,NG                                                       
    XX(K)=RE                                                          
   1  RE=RE-1.0D-04                                                     
  19  J1=J1+1                                                           
DO 2 K=1,NG                                                       
   2  XX(K)=DABS(XX(K))                                                 
DO 3 K=1,NG                                                       
      BB=0.0d0                                                            
DO 4 I=1,NG                                                  
    4 BB=BB+AL(I)*XX(I)*ET(K,I)                                         
   3  RXX(K)=1.0d0/BB                                                     
DO 6 K=1,NG                                                       
   6  XX(K)=(XX(K)+RXX(K))/2.0d0                                          
DO 7 K=1,NG                                                       
IF (DABS(XX(K)-RXX(K)).GT.1.0D-13) GO TO 19                       
   7  CONTINUE   
 ! 	if (AL(4).eq.0.0d0) NG=NG+1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                           
  20  RETURN                                                              
END
!******************************************************************************************************************
SUBROUTINE homogsol(ghom,ghom_forrmix)									!homogeneous phase     
! output ghom - chem pot of segment in hom.phase (bulk)
IMPLICIT NONE
character*7 name_seg(150)
integer NK,NG,Nmol,type_mol(150),M,k,I,s,khom,b,chemm(150),Nbonds,rank_seg(150)
real*8 ghom(150),chempart(150),A,BB,res,teta_ch(150)
REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150)
REAL*8 W(10,10),WT(10,10),WTT(10,10),REDT(10,10)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ,UM,VEC,HEC,GMC,GGA(150)
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150),GPmol(150)
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 teta_chem(150,150) 
real*8 chem(150),sum_Q(150),sum_R(150),sum_seg(150),ghom_formix(150),ghom_forrmix(150)
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /B/ W,WT,WTT,REDT                                 
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                      
COMMON /H/ BLIS,AL,ALI,ET,GP                                                                        
COMMON /Z/ UM,VEC,HEC,GMC,GGA
common /chem/chemm,type_mol,name_seg,teta_chem,rank_seg
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
common /kind/khom
! pure segment k  residual reference act coef.- GP(k)
khom=1

DO  I=1,NK+Nmol
!num_seg(I)=num_segg(I)	  !real from integer
chem(I)=chemm(I)		  !real from integer
sum_Q(I)=0.0d0
sum_R(I)=0.0d0
sum_seg(I)=0.0d0
GPmol(I)=0.0d0 
enddo
DO  k=1,NK
	 if (Nmol.ne.0) then
	   do M=1,Nmol
	    if (type_mol(k).eq.M) then
		sum_Q(M)=sum_Q(M)+Q(k)
		sum_R(M)=sum_R(M)+R(k)
		sum_seg(M)=sum_seg(M)+1.0d0
		endif
	   enddo
	 endif
	enddo

if (Nmol.ne.0) then
	   do M=1,Nmol
		DO  s=1,NG                                                       
		AL(s)=0.0d0
		enddo
DO  k=1,NK
 if (type_mol(k).eq.M) then
 DO  s=1,NG                                                       
 AL(s)=AL(s)+QS(k,s)/sum_Q(type_mol(k))      
 enddo
 endif
enddo                                        
CALL IKS

 DO  k=1,NK
 if (type_mol(k).eq.M) then
	 DO  s=1,NG                                                      
	 GPmol(M)=GPmol(M)-ZZ*QS(k,s)*DLOG(XX(s))
	 enddo
endif
enddo
	  enddo !Nmol

 DO  k=1,NK
if (type_mol(k) /= 0) then  
  A=R(k)/sum_R(type_mol(k))                                                         
 BB=A*sum_Q(type_mol(k))/Q(k)
GPmol(type_mol(k))=GPmol(type_mol(k))+ZR*Q(k)*(1.0d0-BB+DLOG(BB))-(1.0d0-A+DLOG(A))		!num_seg(k)/sum_seg(type_mol(k))
GPmol(type_mol(k))=GPmol(type_mol(k))+chem(k)/2.0d0*dlog(1.0d0/ZZ/sum_Q(type_mol(k)))
!print *,ZR*Q(k)*DLOG(BB)-DLOG(A*num_seg(k)),chem(k)/2.0d0*dlog(sum_chem(type_mol(k))/sum_Q(type_mol(k))),k
endif !type_mol(k)
enddo
 endif !Nmol
b=0
DO  k=1,NK
DO 16000 s=1,NG                                                        
16000 AL(s)=QS(k,s)/Q(k)      !NK - water                                         
CALL IKS

GP(k)=0.0d0     
 DO  s=1,NG                                                      
GP(k)=GP(k)-ZZ*QS(k,s)*DLOG(XX(s))
 enddo
if (type_mol(k).ne.0) GP(k)=GP(k)+1.0d0/2.0d0*chem(k)*dlog(1.0d0/ZZ/Q(k))		  
if (type_mol(k).eq.0) then
b=b+1
GPmol(Nmol+b)=GP(k)
endif  
enddo


call total(X,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
 DO  I=1,NK                                                    
	if (type_mol(I).ne.0) teta_ch(I)=X(I)/QT			   !without ZZ
 enddo

 do	I=1,NK
 if (teta_ch(I).le.0.0d0) then
 chempart(I)=0.0d0
 else
 chempart(I)=-chem(I)/2.0d0*dlog(teta_ch(I)/ZZ) 
 endif
 enddo
 call comb_part(R,Q,NK,RT,QT,ZR,GGA)		   !GGA=combinatorial part of chem pot 
 
DO 3 s=1,NG                                                       
      A=0.0d0                                                             
DO 4 k=1,NK                                                       
   4  A=A+QS(k,s)*X(k)                                                 
   3  AL(s)=A/QT                                                        
CALL IKS  
                                                        
DO  k=1,NK 
res=0.0d0                                                      
DO 6 s=1,NG  
res=res+ZZ*QS(k,s)*DLOG(XX(s))                                                  
   6  GGA(k)=GGA(k)+ZZ*QS(k,s)*DLOG(XX(s))              !residual part  
   enddo
 			 
	b=0
	DO  k=1,NK
    ghom(k)=GGA(k) +chempart(k)+GP(k)+dlog(X(k)) 
	if (type_mol(k) /= 0) then
	ghom_formix(k)=GGA(k)+chempart(k)+dlog(X(k)) 
	else
	b=b+1
	ghom_formix(k)=GGA(k)+chempart(k) +GPmol(Nmol+b)+dlog(X(k))
	endif
	enddo
  khom=0

	DO  k=1,NK+Nmol
  ghom_forrmix(k)=0.0d0
  enddo
 	DO  k=1,NK
	 if (Nmol.ne.0) then
	   do M=1,Nmol
	    if (type_mol(k).eq.M) ghom_forrmix(M)=ghom_forrmix(M)+ghom_formix(k) !*num_seg(k)
	   enddo
	 endif
	enddo
   do M=1,Nmol
	ghom_forrmix(M)=ghom_forrmix(M)+GPmol(M)
	enddo

b=0
	DO  k=1,NK
if (type_mol(k).eq.0) then
b=b+1
ghom_forrmix(Nmol+b)=ghom_formix(k) !*num_seg(k)
endif
	enddo
!	 print *,ghom_forrmix
!	 stop
return
end
!************************************************************************************************************
Subroutine inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)									!calc equations for layer L in inhomogeneous system
! input:Nx - NK-1
!  		Nalf- 2(NG-1)
!       PAR - Xloc,AL_up,chem_up, AL_dwn
!       L-layer
!       ghom - chem pot in bulk phase (from micelle_energy)
! output FU - chem pot of segment in layer L plus work of transfer minus chem pot(bulk)  1,NK-1				
!             AL_up in layer L  NK,NK-1+NG-1
!			  chem_up in layer L  NK-1+NG,NK-1+(NG-1)+Nbonds               
!			  AL_dwn in layer L	 NK-1+NG+Nbonds,NK-1+2*(NG-1)+Nbonds 
IMPLICIT NONE
CHARACTER*3 shape
character*7 name_seg(150)
integer L,M,NK,NG,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150)
integer I,s,Nx,Nalf,k,khom,phase,Lmax,LL,iter,final,final1,seg_ch,seg_max_ch
real*8 m_up,m_dwn,lat,TL(150),resup,resdwn,resadd,FUtau,loc_tens(150),tension,grand_potential,delta,Rout 
REAL*8 PAR(150),FU(150),ghom(150),XX_uphom(150),sumX(150),pair
REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150)
REAL*8 W(10,10),WT(10,10),WTT(10,10),REDT(10,10)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref,TLmin
REAL*8 ZZ,ZR,Phi,IJ,UM,VEC,HEC,GMC,GGA(150)
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150)
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150),XlocL(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ghom2(150),ALchemold(150,150)
REAL*8 Xold(150,150),AL_upold(10,150),AL_dwnold(10,150),AL_alold(10,150),ALlocold(10,150),A_updwnsum(150,150,150),A_alsum(150,150,150),Ach_updwnsum(150,150),Ach_alsum(150,150)
REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150),chempart_segm(150),chempart_short(150),GAMshort(150,150),pN(150) 
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /B/ W,WT,WTT,REDT                                 
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /D/ ZZ,ZR,NK,NG,Phi,IJ,Nmol,Nbonds                                      
COMMON /H/ BLIS,AL,ALI,ET,GP                                                                        
COMMON /Z/ UM,VEC,HEC,GMC,GGA
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al
common /layer/shape,TL,Lmax,TLmin
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
common /old/Xold,AL_upold,AL_dwnold,AL_alold,ALlocold
common /ph/phase			!1 - L+1 reflects L, 0 - inhom, -1 - L=1 (or L-1 is homog)
common /kind/khom
common /hom/XX_uphom
common /pairs/A_updwnsum,A_alsum,Ach_updwnsum,Ach_alsum
common /cycles/iter,final,final1
common /forcalc/GAMshort,chempart_short,ghom2
common /without_or/ Nfix,No_or,CHANGE
common /pressure/pN,dpT
integer Nfix,No_or(150),CHANGE(150)
real*8 chempartt(150,150),eq1(150),eq2(150),res(150),A,BB  
real*8 A_up(10,10,150),A_dwn(10,10,150),A_al(10,10,150),dpT(150,150),comb_mol(150),res_mol(150),chempart_mol(150),transfer_mol(150),pure_mol(150),total_mol(150)    
 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
common /chempot/comb_mol,res_mol,chempart_mol,transfer_mol,pure_mol,total_mol      
!write(*,*)'start of inhomsol'                                                       
khom=0
!writing Xloc and AL_or from PAR
	Xloc(NK,L)=0.0d0   
		DO  I=1,Nx                                          !setting values from PAR          
        Xloc(I,L)=PAR(I)   
	!	Xloc(NK,L)=Xloc(NK,L)+Xloc(I,L)	   !water
		Xloc(NK,L)=Xloc(NK,L)+Xold(I,L)-Xloc(I,L)
		enddo
!	Xloc(NK,L)=1.0d0-Xloc(NK,L)	
	Xloc(NK,L)=Xold(NK,L)+Xloc(NK,L)
		if (Xloc(NK,L).le.0.0d0) then
		print *,'Xloc(NK)<0 in inhomsol',Xloc(NK,L)
		stop
		elseif (Xloc(NK,L).gt.1.0d0) then
		print *,'Xloc(NK)>1 in inhomsol',Xloc(NK,L)
		stop
		endif
   

  !print *,QT,XlocL(1:NK)
!  call comb_part(R,Q,NK,RT,QT,ZR,GGA)		   !GGA=combinatorial part of chem pot 	   !,XlocL
if (Nmol /= 0) then
		do M=1,Nbonds
		chem_up(M,L)=PAR(Nx+Nalf/2+M)
		if (chem_up(M,L).le.0.0d0) then 
		print *,'chemup',M,L,chem_up(M,L)
		endif
		enddo
!calculating chem_dwn (dependant variables)
DO LL=L,Lmax
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
	if (chem_dwn(k,LL+1).le.0.0d0) then 
		print *,'chemdwn',k,LL+1,chem_dwn(k,LL+1)
		endif
enddo
ENDDO	  !LL
endif	  !Nmol
   DO k=1,NK
  XlocL(k)=Xloc(k,L)
  ENDDO
  call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
  call local_chem(Xloc,L,QT) !,Nbonds,L,QT)		!writes teta_chem in common

  DO s=1,NG
  ALloc(s,L)=0.0d0
  DO k=1,NK
  ALloc(s,L)=ALloc(s,L)+QS(k,s)*Xloc(k,L)/QT
  enddo
  enddo
  !print *,ALloc(1,L),ALloc(2,L),ALloc(3,L),QT 		
		DO  I=NK,Nx+Nalf/2		 
	    AL_up(I-Nx,L)=PAR(I)									 !setting values from PAR
		   if (AL_up(I-Nx,L).le.0.0d0) then 
		print *,'AL_up',I-Nx,L,AL_up(I-Nx,L)
		endif
		enddo
 	if ((L.eq.1).and.(r_solid.lt.1.0d-50)) then   
			DO  I=1,NG
			AL_dwn(I,L)=0.0d0
			enddo
		else
		DO  I=NK+Nalf/2+Nbonds,Nx+Nalf+Nbonds
		AL_dwn(I-Nx-Nalf/2-Nbonds,L)=PAR(I)
		 if (	AL_dwn(I-Nx-Nalf/2-Nbonds,L).le.0.0d0) then 
		print *,'AL_dwn',I-Nx-Nalf/2-Nbonds,L,AL_dwn(I-Nx-Nalf/2-Nbonds,L)
		endif
	    enddo
		endif
!print *,AL_up(1,L),AL_up(2,L),AL_up(3,L) 
!print *,'!',AL_al(1,L),Alloc(1,L) !chem_dwn(1:Nbonds,L+1)  !AL_up(1,L),AL_dwn(1,L) !AL_al(2,L),AL_al(3,L) 

call calcAL_al_and_NG(L)  	  !calc dependant variables
!print *,'!',AL_al(1,L),Alloc(1,L) !,chem_dwn(1:Nbonds,L+1)  !,AL_up(1,L),AL_dwn(1,L) !,AL_al(2,L),AL_al(3,L) 
if (final1.eq.0) then
 	 do k=1,NK
	 Xloc(k,Lmax+1)=Xloc(k,Lmax)
	 enddo !NK
	do k=1,Nbonds
	 teta_chem(k,Lmax+1)=teta_chem(k,Lmax)
	chem_up(k,Lmax+1)=chem_dwn(k,Lmax)
	 chem_dwn(k,Lmax+1)=chem_up(k,Lmax)
	 chem_al(k,Lmax+1)=chem_al(k,Lmax) 
	enddo  !Nbonds
	do s=1,NG
	ALloc(s,Lmax+1)=ALloc(s,Lmax)
	AL_al(s,Lmax+1)=AL_al(s,Lmax)
	AL_dwn(s,Lmax+1)=AL_up(s,Lmax)
	AL_up(s,Lmax+1)=AL_dwn(s,Lmax)
	enddo
  if (L.eq.Lmax) then
	call IKS_or(1,L)
else
call calcAL_al_and_NG(L+1)
if (L.eq.1)  call IKS_or(-1,L) 
if (L.ne.1) call IKS_or(0,L)
endif  !Lmax
else
 DO LL=1,Lmax
 !DO LL=L,L+1
 
 !call calcAL_al_and_NG(L+1)
 !if (L.eq.2) print *,AL_up(1,L),AL_up(2,L),AL_up(3,L) 
 !if (L.eq.Lmax) then
 call calcAL_al_and_NG(LL)
	 do k=1,NK
	 !!Xloc(k,Lmax+1)=Xloc(k,Lmax)
	 XlocL(k)=Xloc(k,Lmax+1)
	 enddo !NK
	  call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
  call local_chem(Xloc,Lmax+1,QT) !,Nbonds,L,QT)		!writes teta_chem in common

  DO s=1,NG
  ALloc(s,Lmax+1)=0.0d0
  DO k=1,NK
  ALloc(s,Lmax+1)=ALloc(s,Lmax+1)+QS(k,s)*Xloc(k,Lmax+1)/QT
  enddo
  enddo

	do k=1,Nbonds
	 teta_chem(k,Lmax+1)=teta_chem(k,Lmax)
	chem_up(k,Lmax+1)=0.0d0 !chem_dwn(k,Lmax)
	 chem_dwn(k,Lmax+1)=chem_up(k,Lmax)
	 chem_al(k,Lmax+1)=chem_al(k,Lmax) 
	enddo  !Nbonds
	do s=1,NG
!!	ALloc(s,Lmax+1)=ALloc(s,Lmax)
!	AL_al(s,Lmax+1)=AL_al(s,Lmax)
	AL_dwn(s,Lmax+1)=AL_up(s,Lmax)
	AL_up(s,Lmax+1)= 0.0d0 !AL_dwn(s,Lmax)
	call calcAL_al_and_NG(Lmax+1)
	enddo
  if (LL.eq.Lmax) then
	call IKS_or(1,LL)
else
call calcAL_al_and_NG(LL+1)
if (LL.eq.1)  call IKS_or(-1,LL) 
if (LL.ne.1) call IKS_or(0,LL)
endif  !Lmax
ENDDO  !LL
endif
!LL=L
!do L=LL,LL+1
!print *,'before'  		              
!equations for chem.pot. 
50 call comb_part(R,Q,NK,RT,QT,ZR,GGA)		   !GGA=combinatorial part of chem pot 
Phi=R(NK-1)/RT 
                                              
call chemparts(L,chempart_segm,chempart_short)	!chemical part of chem pot
!print *,'after'
resadd=0.0d0
 DO s=1,NG
 !print *,XX_up(s,L)*AL_up(s,L),XX_al(s,L),AL_al(s,L),XX_dwn(s,L)*AL_dwn(s,L),XX_al(s,L)*AL_al(s,L),s
   IF  (L.ne.1) then
   resadd=resadd-Alloc(s,L)*ZZ*(m_up(L)*dlog(XX_up(s,L)*AL_up(s,L)/(XX_al(s,L)*AL_al(s,L)))+m_dwn(L)*dlog(XX_dwn(s,L)*AL_dwn(s,L)/(XX_al(s,L)*AL_al(s,L))))
   else
   resadd=resadd-Alloc(s,L)*ZZ*m_up(L)*dlog(XX_up(s,L)*AL_up(s,L)/(XX_al(s,L)*AL_al(s,L)))
   endif
	enddo
	 do M=1,Nbonds
	  IF  (L.ne.1) then
   resadd=resadd-(teta_chem(M,L)+teta_chem(M,L)*Xloc(seg_max_ch(M),L)/Xloc(seg_ch(M),L))*(m_up(L)*dlog(XX_up(NG,L)*AL_up(NG,L)/(XX_al(NG,L)*AL_al(NG,L)))+m_dwn(L)*dlog(XX_dwn(NG,L)*AL_dwn(NG,L)/(XX_al(NG,L)*AL_al(NG,L))))
   else
   resadd=resadd-(teta_chem(M,L)+teta_chem(M,L)*Xloc(seg_max_ch(M),L)/Xloc(seg_ch(M),L))*m_up(L)*dlog(XX_up(NG,L)*AL_up(NG,L)/(XX_al(NG,L)*AL_al(NG,L)))
   endif
   enddo
 !if (final)  print *,resadd,L
 FUtau=0.0d0
DO k=1,NK
res(k)=0.0d0	 !residual part	of chem pot
!resup=0.0d0
!resdwn=0.0d0
!resal=0.0d0
   DO s=1,NG
   IF  (L.ne.1) then
   res(k)=res(k)+ZZ*QS(k,s)*(m_up(L)*dlog(XX_up(s,L)*AL_up(s,L)/Alloc(s,L))+m_dwn(L)*dlog(XX_dwn(s,L)*AL_dwn(s,L)/Alloc(s,L))+lat(L)*dlog(XX_al(s,L)*AL_al(s,L)/Alloc(s,L))  )!-dlog(Alloc(s,L)))
  ! res=res-Alloc(s,L)*ZZ*Q(k)*(m_up(L)*dlog(XX_up(s,L)*AL_up(s,L)/(XX_al(s,L)*AL_al(s,L)))+m_dwn(L)*dlog(XX_dwn(s,L)*AL_dwn(s,L)/(XX_al(s,L)*AL_al(s,L))))
   else	   
    res(k)=res(k)+ZZ*QS(k,s)*(m_up(L)*dlog(XX_up(s,L)*AL_up(s,L)/Alloc(s,L))+lat(L)*dlog(XX_al(s,L)*AL_al(s,L)/Alloc(s,L)) ) !-dlog(Alloc(s,L)))
	! res=res-Alloc(s,L)*ZZ*Q(k)*m_up(L)*dlog(XX_up(s,L)*AL_up(s,L)/(XX_al(s,L)*AL_al(s,L)))
  endif
  IF ((L.eq.1).and.(r_solid.gt.1.0d-50)) then
    res(k)=res(k)+ZZ*QS(k,s)*m_dwn(L)*dlog(AL_dwn(s,L)/(Alloc(s,L)*ET_solid(s)))  !YYY	
  endif

 ! FUtau=FUtau-ZZ*QS(k,s)*Xloc(k,L)*dlog(XX_al(s,L)*AL_al(s,L)/Alloc(s,L))
   enddo
  res(k)=res(k)+resadd*Q(k)
GAMloc(k,L)=GGA(k)+res(k)+GP(k) +chempart_segm(k)+dlog(Xloc(k,L)) !+resadd*Q(k)		!chem pot
GAMshort(k,L)=GGA(k)+res(k)+GP(k) +chempart_short(k)+dlog(Xloc(k,L)) !+resadd*Q(k)
!FUtau=FUtau+(ghom(k)-ghom(NK)*R(k)/R(NK)-GP(k)-(GGA(k)+dlog(Xloc(k,L)-GAMloc(NK,L)*R(k)/R(NK))))*Xloc(k,L)	
enddo

comb_mol=0.0d0
res_mol=0.0d0
chempart_mol=0.0d0
transfer_mol=0.0d0
pure_mol=0.0d0
total_mol=0.0d0
do M=1,Nmol
  do k=1,nk
  if (type_mol(k).eq.M) then
  comb_mol(M)=comb_mol(M)+GGA(k)+dlog(Xloc(k,L))
  res_mol(M)=res_mol(M)+res(k)
  chempart_mol(M)=chempart_mol(M)+chempart_short(k)
  transfer_mol(M)=transfer_mol(M)+chempart_segm(k)-chempart_short(k)
  pure_mol(M)= pure_mol(M)+GP(k)
  total_mol(M)=total_mol(M)+GAMloc(k,L)-GAMloc(nk-1,L)*R(k)/R(NK-1)
  endif
  enddo
enddo
M=Nmol
 do k=1,nk
  if (type_mol(k).eq.0) then
  M=M+1
  comb_mol(M)=comb_mol(M)+GGA(k)+dlog(Xloc(k,L))
  res_mol(M)=res_mol(M)+res(k)
  chempart_mol(M)=chempart_mol(M)+chempart_short(k)
  transfer_mol(M)=transfer_mol(M)+chempart_segm(k)-chempart_short(k)
  pure_mol(M)= pure_mol(M)+GP(k)
  total_mol(M)=total_mol(M)+GAMloc(k,L)
  if (k.ne.nk-1) total_mol(M)=total_mol(M)-GAMloc(nk-1,L)*R(k)/R(NK-1)
  endif
  enddo

!DO k=1,NK
!FUtau=FUtau+(ghom(k)-GP(K)-(GGA(k)+dlog(Xloc(k,L))))*Xloc(k,L)	
!FUtau=FUtau+(GAMloc(k,L)-GAMloc(NK,L)*R(k)/R(NK)+ghom(NK)*R(k)/R(NK)-GP(K)-(GGA(k)+dlog(Xloc(k,L))))*Xloc(k,L)	
!enddo
!if (shape.eq.'sph')  FUtau= FUtau/RT
! DO s=1,NG
 !if (shape.eq.'sph')  FUtau= FUtau-2.0d0/Rout(L)*AL_up(s,L)*dlog(AL_up(s,L)*XX_up(s,L)/AL_al(s,L)/XX_al(s,L))
 !enddo

FUtau=pN(L)-pN(Lmax) !+ghom(NK)/R(NK) !-pN(Lmax)
 if ((shape.eq.'sph').and. (L.ne.Lmax)) then
 !CALL grandpot_pressure(L,ghom,loc_tens,tension,grand_potential,delta)
 !FUtau=pN(L)-pN(L+1)+2.0d0*(-loc_tens(L)-pN(L)+pN(Lmax))*dlog(Rout(L+1)/Rout(L))
 !	-normall(L+1)+normall(L)+2.0d0*(tang(L)-normall(L))*dlog(Tcurrent+1/(Tcurrent)))
!! FUtau= FUtau !+ghom(NK)/R(NK)*Rout(Lmax)**2/Rout(L)**2
 do LL=L+1,Lmax
!!!CALL grandpot_pressure(LL,ghom,loc_tens,tension,grand_potential,delta)
!!!FUtau= FUtau-loc_tens(LL)*(Rout(LL)**2-Rout(LL-1)**2)/Rout(L)**2
!FUtau= FUtau+dpT(LL,L)
  enddo
!else
!FUtau= FUtau+ghom(NK)/R(NK)
  endif
   ! print *,'pn,fu=',L,pN(L),pN(Lmax),FUtau,TL(L)
!print *,L,FUtau
!print *,'before'  
DO k=1,Nx
FU(k)=(GAMloc(k,L)-GAMloc(NK,L)*R(k)/R(NK))-(ghom(k)-ghom(NK)*R(k)/R(NK)) !**2 +(GAMloc(NK,L)-ghom(NK))**2	!(GAMloc(k,L)-GAMloc(NK,L))-(ghom(k)-ghom(NK)) !	  
enddo
!FU(nx)=(GAMloc(nk,L)-GAMloc(NK-1,L)*R(nk)/R(NK-1))-(ghom(nk)-ghom(NK-1)*R(nk)/R(NK-1))
! print *,'after'
 !print *,ghom(NK-1)-GAMloc(NK-1,L)
!equations for orientations
 DO s=2,NG
eq1(s-1)=XX_up(s,L)*AL_up(s,L)*XX_al(s-1,L)*Al_al(s-1,L)-XX_al(s,L)*XX_up(s-1,L)*AL_al(s,L)*AL_up(s-1,L)  
 IF  (L.ne.1) then
 eq2(s-1)=XX_dwn(s,L)*XX_al(s-1,L)*AL_dwn(s,L)*AL_al(s-1,L)-XX_al(s,L)*XX_dwn(s-1,L)*AL_al(s,L)*AL_dwn(s-1,L)
 else if (r_solid.gt.1.0d-50) then
 eq2(s-1)=XX_al(s-1,L)*AL_dwn(s,L)*Al_al(s-1,L)*ET_solid(s-1)-XX_al(s,L)*AL_dwn(s-1,L)*AL_al(s,L)*ET_solid(s)  
 endif
 enddo

 if (Nmol /= 0) then 
 DO M=1,Nbonds
 !eq1(NG-1+M)=dlog(chem_up(M,L)*chem_dwn(M,L+1))+2.0d0*dlog(XX_al(NG,L)*AL_al(NG,L)*XX_al(NG,L+1)*AL_al(NG,L+1))-dlog(chem_al(M,L)*chem_al(M,L+1))-2.0d0*dlog(XX_up(NG,L)*AL_up(NG,L)*XX_dwn(NG,L+1)*AL_dwn(NG,L+1))
!if (L.ne.1) eq1(NG-1+M)=2.0d0*Wprev(M)+dlog(chem_dwn(M,L+1)/chem_up(M,L))
 eq1(NG-1+M)=chem_up(M,L)*chem_dwn(M,L+1)*(XX_al(NG,L)*AL_al(NG,L)*XX_al(NG,L+1)*AL_al(NG,L+1))**2-chem_al(M,L)*chem_al(M,L+1)*(XX_up(NG,L)*AL_up(NG,L)*XX_dwn(NG,L+1)*AL_dwn(NG,L+1))**2
 !if (Nfix.eq.0) eq1(NG-1+M)=chem_up(M,L)*chem_dwn(M,L+1)*(XX_al(NG,L)*AL_al(NG,L)*XX_al(NG,L+1)*AL_al(NG,L+1))**2-chem_al(M,L)*chem_al(M,L+1)*(XX_up(NG,L)*AL_up(NG,L)*XX_dwn(NG,L+1)*AL_dwn(NG,L+1))**2
 !if (Nfix) eq1(NG-1+M)=chem_up(M,L)*chem_dwn(M,L+1)*(XX_al(No_or(1),L)*AL_al(No_or(1),L)*XX_al(NG,L+1)*AL_al(NG,L+1))**2-chem_al(M,L)*chem_al(M,L+1)*(XX_up(No_or(1),L)*AL_up(No_or(1),L)*XX_dwn(NG,L+1)*AL_dwn(NG,L+1))**2
 ENDDO
 endif

 DO s=NK,Nx+NG-1+Nbonds
 FU(s)=eq1(s-Nx) 
 !write (108,'(I1,(4X,E12.6),4X,I1)')s-Nx,FU(s),L
! write(*,*)'FUalfup=',s-Nx,FU(s),L
 enddo

 if ((L.ne.1).or.(r_solid.gt.1.0d-50)) then
   DO s=Nx+NG+Nbonds,Nx+Nalf+Nbonds
   FU(s)=eq2(s-Nx-NG-Nbonds+1)   
  !  write (108,'(I1,(4X,E12.6),4X,I1)')s-Nx-NG+1,FU(s),L
   enddo
 endif

  if (final /= 0) then	  !information for final profiles
 DO s=1,NG
 DO k=1,NG
 A_up(s,k,L)=pair(s,L,'up',k)
 A_dwn(s,k,L)=pair(s,L,'dw',k)	 !down
 A_al(s,k,L)=pair(s,L,'al',k)
 enddo
 enddo
 call sum_pairs(A_up,A_dwn,A_al,L,NG,A_updwnsum,A_alsum)
 call sum_chempairs(chem_up,chem_dwn,chem_al,L,Lmax,Nbonds,Ach_updwnsum,Ach_alsum)
 endif
 !write(108,*)'!!!!!!!!!!!'
 !write(*,*)'end of inhomsol'
return
end
!************************************************************************************************************
Subroutine micelle_energy(sh,thick,r_max,NL,Xin,waterLmax,tot_energy)									
! input: shape - cyl or sph or pln
!!!		TL- layer's thickness
!       NL - number of layers
!       Xin - input for LVE
! output tot_energy
IMPLICIT NONE
CHARACTER*3 shape,sh 
integer L,Lmax,KEYit,iter,M,KEY,NK,NG,LAY,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150),bond_seg(150),seg_bond(150),NKK,KEYit1
integer J,I,s,Nx,Nalf,k,switch,kstop,phase,ii,jj,final,final1,NL,waterLmax,seg_max_ch,seg_ch,info(1)
real*8 loc_energy(150),loc_tens(150),grand_potential,tenss,tot_energy,energy2,tension,ALlocnew(150),sumXX(150),ALchemold(150,150)
real*8 pi,FUtau,delta,TLmin,r_max, FU_LVEP
real*8 Xsegtot,Xseg(100),g(100),gmx,thick(150),m_up,m_dwn,lat,TL(150),TLold(150)   
REAL*8 PAR(150),FU(150),DD(150,150),PARx(150),PARalf(150),DPx(150),DPalf(150),ghom(150),EPS,EPF,ghom2(150)
REAL*8 Xold(150,150),X(150),Xin(150),Xout1(150),Xout2(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150),PAR1(150),Xlocc(150,150)
REAL*8 W(10,10),WT(10,10),WTT(10,10),REDT(10,10)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ,UM,VEC,HEC,GMC,GGA(150)
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150),GAMshort(150,150),chempart_short(150)
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150),ghom_formix(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150)
REAL*8 AL_upold(10,150),AL_dwnold(10,150),AL_alold(10,150),ALlocold(10,150),A_updwnsum(150,150,150),A_alsum(150,150,150),Ach_updwnsum(150,150),Ach_alsum(150,150)
REAL*8 teta_chem(150,150),Cmpsi(300),chem_up(150,150),chem_dwn(150,150),chem_al(150,150),chem_upold(150,150)
real*8 Falfeq,PARalfeq(150),Fxeq,PARxeq(150),Xbulk(150),Rout,Y1(1),Ymin(1),pN(150),Pnormal,grand_potold
character*7 name_seg(150)
character*30 name_mol(150)
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /B/ W,WT,WTT,REDT                                 
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                      
COMMON /H/ BLIS,AL,ALI,ET,GP                                                                         
COMMON /Z/ UM,VEC,HEC,GMC,GGA
common /layer_one/PAR1
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al
common /layer/shape,TL,Lmax,TLmin
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
common /old/Xold,AL_upold,AL_dwnold,AL_alold,ALlocold
common /oldch/chem_upold
common /ph/phase				  !1 - L+1 reflects L, 0 - inhom, -1 - L=1 (or L-1 is homog)
common /conc/Cmpsi,kstop
common /pairs/A_updwnsum,A_alsum,Ach_updwnsum,Ach_alsum
common /cycles/iter,final,final1
common /forcalc/GAMshort,chempart_short,ghom2
common /bond_to_segment/bond_seg,seg_bond
common /eq/Falfeq,PARalfeq,Fxeq,PARxeq
common /energy/energy2
common /forNEWT/Nx,Nalf,DPx,DPalf,EPS,EPF
common /composition/Xbulk
 common /tau_old/TLold
 common /equib/Xinn,NNK,NKK,P,PP,Temp
external inhomsol,LVE_flash,mol_prop,FU_LVEp 
real*8 P,PP,Xinn(150),Temp,dpT(150,150) 
common /pressure/pN,dpT
real*8 comb_mol(150),res_mol(150),chempart_mol(150),transfer_mol(150),pure_mol(150),total_mol(150)    
common /chempot/comb_mol,res_mol,chempart_mol,transfer_mol,pure_mol,total_mol 
integer NNK
   pi = 3.1415926d0
!OPEN(UNIT=33,FILE='quasi.txt')
OPEN(UNIT=34,FILE='Xloc.txt')
OPEN(UNIT=35,FILE='Alfa_up.txt')
OPEN(UNIT=36,FILE='Alfa_dwn.txt')
OPEN(UNIT=37,FILE='Alfa_al.txt')
OPEN(UNIT=38,FILE='Alfa_loc.txt')
OPEN(UNIT=39,FILE='XX_or.txt')
!OPEN(UNIT=55,FILE='A_up.txt')
!OPEN(UNIT=56,FILE='A_al.txt')
OPEN(UNIT=110,FILE='graph.txt')
OPEN(UNIT=111,FILE='graph_pofile.txt')
OPEN(UNIT=120,FILE='Final_FU.txt')
OPEN(UNIT=121,FILE='Final_X.txt')
OPEN(UNIT=122,FILE='final chempot.txt')
OPEN(UNIT=123,FILE='final chempot_parts.txt')
OPEN(UNIT=135,FILE='Alfa_upfinal.txt')
OPEN(UNIT=136,FILE='Alfa_dwnfinal.txt')
OPEN(UNIT=137,FILE='Alfa_alfinal.txt')
OPEN(UNIT=138,FILE='Alfa_locfinal.txt')
!OPEN(UNIT=139,FILE='Alfa_totfinal.txt')
OPEN(UNIT=155,FILE='pairs_up+down_final.txt')
OPEN(UNIT=156,FILE='pairs_lat_final.txt')
OPEN(UNIT=200,FILE='local_tension.txt')
OPEN(UNIT=201,FILE='tension(r).txt')
OPEN(UNIT=202,FILE='chempot_solvent.txt')
OPEN(UNIT=210,FILE='local_free_energy.txt')
OPEN(UNIT=220,FILE='potential_water.txt')
!OPEN(UNIT=221,FILE='tension.txt')
OPEN(UNIT=112,FILE='init profiles.txt')
OPEN(UNIT=113,FILE='current profiles.txt')
OPEN(UNIT=114,FILE='current chempot.txt')
OPEN(UNIT=115,FILE='phase equilibrium.txt')
!OPEN(UNIT=114,FILE='parabola.txt')
Temp=T
NNK=NK
Xinn=Xin	 
Lmax=NL ! number of layers in system
phase = 0 !inhom layers

shape=sh
TL=thick
TLmin=0.2d0 !r_max

NKK=Nmol						       !NKK=N(chains)+N(monomeric molecules)
	  DO I=1,NK
	  if (type_mol(I).eq.0) NKK=NKK+1 
	  enddo
	 ! name_mol=name_seg
	  do I=1,NK
	  if (type_mol(I).ne.0)	then
	  jj=type_mol(I)
		if (rank_seg(I).eq.1)  then
		name_mol(jj)=trim(name_seg(I))
		else
		name_mol(jj)=trim(name_mol(jj))//trim(name_seg(I))
		endif
	  else
	  jj=jj+1
	   name_mol(jj)=trim(name_seg(I))
	  endif
	 ! if (name_mol(jj).ne.name_seg(I)) then
	 ! name_mol(jj)=name_mol(jj)//name_seg(I)
	 ! else
	 ! name_mol(jj)=name_seg(I)
	 ! endif
	  enddo

NKK=Nmol						       !NKK=N(chains)+N(monomeric molecules)
	  DO I=1,NK
	  if (type_mol(I).eq.0) NKK=NKK+1 
	  enddo
Call mol_prop(NKK,type_mol,Xout1,Xout2)		! Xout1 is the aqueous phase, Xout2 is the organic phase

!if (NKK.ge.3) then
!Y1(1)=Xin(NKK-1)		   ! holes
!P=1.01325d5

!call DZREAL (FU_LVEp, 1.0d-10, 1.0d-10, 1.0d-2, 1.0d-2, 1, 1000,Y1, Ymin, INFO)	
!Xin(NKK-1)=Ymin(1)
!if ((NKK-1).ne.1) Xin(1)=1.0d0-Xin(NKK-1)
! if ((NKK-2).ge.2) then
!	  do I=2,NKK-2
!	  Xin(1)=Xin(1)-Xin(I)
!	  enddo
!	  endif
!endif
if (waterLmax.eq.0) then
!XXXXXXX    call LVE(Xin,NK,Xout2,Xout1) 
P=PP 
!call LVE(Xin,NK,Xout2,Xout1)   !calculation of phase equilibrium	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elseif (waterLmax.eq.1) then
!    call LVE(Xin,NK,Xout1,Xout2)
Xbulk=Xout1 
! P=PP 
else
Xout2=Xbulk
Xout1=Xbulk
!call LVE(Xin,NK,Xout1,Xout2)
!P=PP 
!XXXXX   do I=1,NK
!XXX    Xout2(I)=Xbulk(I)
!XXXXXXX   enddo
!Xout1=Xbulk
!XXXXXXXXXXXXXXX    Xout1(1)=0.23822d0
!XXXXXXXXXXXXXXX    Xout1(2)=0.31229d0
!Xout1(3)=0.1d0
!XXXXXXXXXXXXXX    Xout1(3)=1.0d0-Xout1(1)-Xout1(2)
endif
NKK=Nmol						       !NKK=N(chains)+N(monomeric molecules)
	  DO I=1,NK
	  if (type_mol(I).eq.0) NKK=NKK+1 
	  enddo
print *,'compositions of eq phases'
print *,'phase 1',Xout1(1:NKK)
print *,'phase 2',Xout2(1:NKK)
print *,'phase 1',Xout1(1:NKK-2)/(1.0d0-Xout1(NKK-1)),Xout1(NKK)/(1.0d0-Xout1(NKK-1))
print *,'phase 2',Xout2(1:NKK-2)/(1.0d0-Xout2(NKK-1)),Xout2(NKK)/(1.0d0-Xout2(NKK-1))
! stop
write(115,'(10X,A40,/,16X,30(A18,4X))')'compositions of equilibrium phases',name_mol(1:NKK)
write(115,'(A10,30(4X,E18.12))')'phase 1',Xout1(1:NKK)
write(115,'(A10,30(4X,E18.12))')'phase 2',Xout2(1:NKK)
write(115,'(30(6X,A18))')name_mol(1:NKK-2),name_mol(NKK)
write(115,'(30(4X,E18.12))')Xout1(1:NKK-2)/(1.0d0-Xout1(NKK-1)),Xout1(NKK)/(1.0d0-Xout1(NKK-1))
write(115,'(30(4X,E18.12))')Xout2(1:NKK-2)/(1.0d0-Xout2(NKK-1)),Xout2(NKK)/(1.0d0-Xout2(NKK-1))
write(115,'(A2,30(4X,E18.12))')'P=',P
close (115)

!Xout1, Xout2 - molecular compositions of eq phases (for chains and monomeric molecules)
call XcomptoXseg(Xout1,X2)
	  	DO I=1,NK
		X(I)=X2(I)
		enddo  
call homogsol(ghom2,ghom_formix) 					 ! save chem pot of segments in the bulk phases
!Xout2(NKK)=1.0d0
!do I=1,NKK-1
!Xout2(I)=Xout2(I)*0.99d0
!!Xout2(NKK)=Xout2(NKK)-Xout2(I)
!enddo
call XcomptoXseg(Xout2,X)		  !X - segments'composition in bulk phase

call homogsol(ghom,ghom_formix)	  !ghom - chem.pot. of segments in homogeneous (bulk) phase 
!print *,-8.3143d0*Temp*ghom(NK-1)/R(NK-1)/14.244D-06,-8.3143d0*Temp*ghom2(NK-1)/R(NK-1)/14.244D-06

!stop


!guess initial profile of segments' concentrations  Xloc andorientations	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call init_profile(X2,X,Lmax)
! print *,ghom(NK-1)-ghom2(NK-1)
do L=1,Lmax
 phase=0
if (L.eq.1) phase=-1
if (L.eq.Lmax) phase=1
call IKS_or(phase,L)		   !solve quasichemical equations,calc XX_or(L)
enddo
!stop
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculation of FU for initial profiles
 !do k=1,Nx+Nalf+Nbonds
 ! Y(k)=0.0d0
 ! enddo

  do L=1,Lmax				      ! cycle over layers
	 phase=0
	if (L.eq.1) phase=-1
	if (L.eq.Lmax) phase=1
  do k=1,Nx
  PARx(k)=Xloc(k,L)						
  enddo
  do s=1,NG-1
  PARalf(NK-1+s)=AL_upold(s,L) !Alloc(s,L)					
  enddo	

  if (Nmol /= 0) then
 DO M=1,Nbonds
 PARalf(NK-1+NG-1+M)=chem_upold(M,L) !teta_chem(M,L)/ZZ	
 enddo	
  endif
 do s=1,NG-1				
  PARalf(NK+NG-2+Nbonds+s)=AL_dwnold(s,L) !Alloc(s,L)
  enddo	
 
    DO  I=1,Nx+Nalf+Nbonds
	  if (I.le.Nx) then                                                        
      PAR(I)=PARx(I)
	  else
	  PAR(I)=PARalf(I)					!creating unified array PAR
	  endif 
	 enddo
!	  print *,'Alfa=',m_up(L),lat(L),m_up(L+1),lat(L+1)
  CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom) 
  call calcALloc(PARx,L,ALlocnew,sumXX)

 do s=Nalf/2+1,Nalf/2+Nbonds
	 ALchemold(s-Nalf/2,L)=ALlocnew(s+NG-Nalf/2) !chem_up(s-Nalf/2,L) !Alloc(s,L)					
	enddo
	do s=1,NG
	ALlocold(s,L)=ALloc(s,L)
	enddo
 
 ! if (L.eq.Lmax) pN(L)=Pnormal(L,PAR,Nx,Nalf,DPalf,EPS,EPF,ghom)
 write (110,'(A3,I2)')'L=',L
   write (110,'(16(4X,E18.12))')Xloc(nk,L),(Xloc(k,L),k=1,NK-1)
   write (110,'(16(4X,E18.12))')(FU(s),s=1,Nx) 
!	 stop
  ! write (110,'(/)')
   write (110,'(16(4X,E18.12))')(FU(s),s=Nx+1,Nx+Nalf+Nbonds) 
 !  write (110,'(16(4X,E18.12))')(FU(s),s=Nx+Nalf/2+Nbonds+1,Nx+Nalf+Nbonds) 
 !stop
   enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	   ebd of layers cycle
 ! do k=1,Nx+Nalf+Nbonds
 ! Y(k)=FU(k)
 ! enddo
 ! start of iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 TLold=TL
 grand_potold=1.0d0
iter=0
 final1=0
70 final=0
KEYit=1
switch=0
KEYit1=0																						  
7 do while (KEYit /= 0)					       !KEYit=0 when Xloc, Alfa don't change between iter and iter+1	
8 phase=0	
!if ((iter).and.(L.eq.1))	switch=switch-1		!go up
!if ((iter.gt.2).and.(L.eq.Lmax)) switch=switch+1	!go down	   
iter=iter+1
write (110,'(A10,I5)')'iteration=',iter
write (111,'(A10,I5)')'iteration=',iter
KEYit=0
!KEYit1=0
 if (iter.gt.1) then
   do Lay=1,Lmax
   TLold(Lay)=TL(Lay)						       
   do k=1,NK
   Xold(k,Lay)=Xloc(k,Lay)
   enddo 
   do s=1,NG
   AL_upold(s,Lay)=AL_up(s,Lay)
   AL_dwnold(s,Lay)=AL_dwn(s,Lay)
   AL_alold(s,Lay)=AL_al(s,Lay)
   ALlocold(s,Lay)=ALloc(s,Lay)
   enddo 
  if (Nmol /= 0) then
  DO M=1,Nbonds
  chem_upold(M,Lay)=chem_up(M,Lay)	
  enddo	
  endif
  enddo	  !Lay
	do k=1,NK
   Xold(k,Lmax+1)=Xloc(k,Lmax)
   enddo 
   do s=1,NG
   AL_upold(s,Lmax+1)=AL_up(s,Lmax)
   AL_dwnold(s,Lmax+1)=AL_dwn(s,Lmax)
   AL_alold(s,Lmax+1)=AL_al(s,Lmax)
   ALlocold(s,Lmax+1)=ALloc(s,Lmax)
   enddo 
 endif

if (switch.eq.0) L=0					 !  switch=0:   start from 1st layer, switch>0 start from last
!if (switch) L=Lmax+1
gmx=0.0d0
Xsegtot=0.0d0
do M=1,Nmol
g(M)=0.0d0
Xseg(M)=0.0d0
enddo
do 	  !the last layer is Lmax where Xloc=X(homog) 
if (switch.eq.0) L=L+1
!if (switch) L=L-1
phase=0
if (L.eq.1) phase=-1
if (L.eq.Lmax) phase=1
 do k=1,Nx
  PARx(k)=Xloc(k,L)						
  enddo
write(*,*)'L=',L,'iter=',iter   
!call NEWTx_P(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)

!if ((iter.ne.1).and.((L.eq.1).or.(final1.eq.1).or.(shape.eq.'pln')) ) then

!else
!if (final1.eq.0) call NEWTx_canonic(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,DD,L,EPS,EPF,KEY,ghom)	 !if (mod(iter,15).eq.0) 	  !ASASASASASASASASA call NEWTx(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)
!if (L.ne.1) then
call NEWTx(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)
Xlocc=Xloc
if (final1.eq.1) pN(L)=Pnormal(L,Xlocc,Nx,Nalf,DPalf,EPS,EPF,ghom)			! final=0, final=1 -finished;   final1- switch for the normal pressure
 !print *,-ghom(NK-1)
!else
! if (L.eq.Lmax) pN(Lmax)=Pnormal(Lmax,Xloc,Nx,Nalf,DPalf,EPS,EPF,ghom)	  ! setting pressure in last layer for mexh equil calculation, see NEWT_mechanical below

 DO  I=1,Nx+Nalf+Nbonds
	  if (I.le.Nx) then                                                        
      PAR(I)=PARx(I)
	  else
	  PAR(I)=PARalf(I)					!creating unified array PAR
	  endif 
	 enddo
!if ((final1.eq.1).and.(L.lt.(Lmax-1))) call NEWT_mechanical(Nx,Nalf,PAR,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)	  !.and.(iter.ne.1)
!	print *,'Xloc=',Xloc(1,L),L 
!!		DO I=1,NK
!		X(I)=Xloc(I,Lmax)
!		enddo  
!call homogsol(ghom,ghom_formix) 
!if ((mod(iter,15).ne.0).and.(L.ne.Lmax)) call NEWT_mechanical(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)		!	if ((iter.gt.1).and.(L.lt.2))
!endif
!if (shape.eq.'pln') then
!if (final1.eq.0) call NEWTx(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)
!if (final1.eq.1) call NEWTx_P(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)	  ! mech equil from newt-raphson
!if (final1.eq.1) then
!call NEWTP(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)
!else
!if (iter.eq.1) call NEWTx(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)
!if (iter.ne.1) call NEWTx_P(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)
!endif
!call NEWTx(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)
!if (L.eq.Lmax) then
!do J=1,Lmax
 !do k=1,Nx
 ! PARx(k)=Xloc(k,L)						
 ! enddo
 !call NEWTP(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,FUtau,DD,L,EPS,EPF,KEY,ghom)
 !if (dabs(TL(J)-TLold(J))/TLold(J).gt.1.0d-5) KEYit=1	
 !enddo
 !endif

!print *,Al_up(1,L)+Al_up(2,L)+chem_up(1,L)+chem_up(2,L)+chem_up(3,L)+chem_dwn(1,L+1)+chem_dwn(2,L+1)+chem_dwn(3,L+1)
!else
!call NEWTx_canonic(Nx,Nalf,PARx,PARalf,DPx,DPalf,FU,DD,L,EPS,EPF,KEY,ghom)
!endif  !shape
 !call calcALloc(PARx,L,ALlocnew,sumXX)

! do s=Nalf/2+1,Nalf/2+Nbonds
!ALchemold(s-Nalf/2,L)=ALlocnew(s+NG-Nalf/2) !chem_up(s-Nalf/2,L) !Alloc(s,L)					
!enddo
		  
IF(KEY.eq.(-1)) then
write(*,*)'Bad precision after NEWT'
stop
endif
write(*,*)'!!!!!!!!!!!','phase=',phase
CALL grandpot_pressure(L,ghom,loc_tens,tension,grand_potential,delta)

! convergence tests below  (keyit)
if (mod(iter,20).eq.0) then
 if (dabs(TL(L)-TLold(L))/TLold(L).gt.1.0d-3) KEYit=1	
 do k=1,NK
 ! if (dabs(Xloc(k,L)-Xold(k,L))/Xloc(k,L).gt.1.0d-5) print *,'mole frac',k,dabs(Xloc(k,L)-Xold(k,L))/Xloc(k,L)
  if (dabs(Xloc(k,L)-Xold(k,L))/Xloc(k,L).gt.1.0d-3) KEYit=1					!>0.001%
  enddo
  do s=1,NG																		!5.0d-6
  if (dabs(AL_up(s,L)-AL_upold(s,L))/AL_up(s,L).gt.5.0d-5) KEYit=1
   !if (dabs(AL_up(s,L)-AL_upold(s,L))/AL_up(s,L).gt.5.0d-7) print *,'up',s,dabs(AL_up(s,L)-AL_upold(s,L))/AL_up(s,L)
   if (L.ne.1) then
  if (dabs(AL_dwn(s,L)-AL_dwnold(s,L))/AL_dwn(s,L).gt.5.0d-5) KEYit=1			 !5.0d-6
  ! if (dabs(AL_dwn(s,L)-AL_dwnold(s,L))/AL_dwn(s,L).gt.5.0d-7) print *,'dwn',s,dabs(AL_dwn(s,L)-AL_dwnold(s,L))/AL_dwn(s,L)
   endif
  enddo
  if (Nmol /= 0) then
 DO M=1,Nbonds
 if (dabs(chem_up(M,L)-chem_upold(M,L))/chem_up(M,L).gt.5.0d-5) KEYit=1			 !5.0d-6
 !if (dabs(chem_up(M,L)-chem_upold(M,L))/chem_up(M,L).gt.5.0d-7) print *,'chem',M,dabs(chem_up(M,L)-chem_upold(M,L))/chem_up(M,L)
 enddo
 endif
 if ((L.eq.Lmax).and.(dabs(grand_potential-grand_potold)/grand_potold.gt.1.0d-3)) KEYit=1
  if (L.eq.Lmax)	grand_potold=grand_potential
else
 if (dabs(TL(L)-TLold(L))/TLold(L).gt.1.0d-5) KEYit=1	
 !if (iter.lt.3) then 
 do k=1,NK
 ! if (dabs(Xloc(k,L)-Xold(k,L))/Xloc(k,L).gt.1.0d-5) print *,'mole frac',k,dabs(Xloc(k,L)-Xold(k,L))/Xloc(k,L)
  if (dabs(Xloc(k,L)-Xold(k,L))/Xloc(k,L).gt.1.0d-5) KEYit=1					!>0.001%
  enddo
  do s=1,NG																		!5.0d-6
  if (dabs(AL_up(s,L)-AL_upold(s,L))/AL_up(s,L).gt.5.0d-7) KEYit=1
   !if (dabs(AL_up(s,L)-AL_upold(s,L))/AL_up(s,L).gt.5.0d-7) print *,'up',s,dabs(AL_up(s,L)-AL_upold(s,L))/AL_up(s,L)
   if (L.ne.1) then
  if (dabs(AL_dwn(s,L)-AL_dwnold(s,L))/AL_dwn(s,L).gt.5.0d-7) KEYit=1			 !5.0d-6
  ! if (dabs(AL_dwn(s,L)-AL_dwnold(s,L))/AL_dwn(s,L).gt.5.0d-7) print *,'dwn',s,dabs(AL_dwn(s,L)-AL_dwnold(s,L))/AL_dwn(s,L)
   endif
  enddo
  if (Nmol /= 0) then
 DO M=1,Nbonds
 if (dabs(chem_up(M,L)-chem_upold(M,L))/chem_up(M,L).gt.5.0d-7) KEYit=1			 !5.0d-6
 !if (dabs(chem_up(M,L)-chem_upold(M,L))/chem_up(M,L).gt.5.0d-7) print *,'chem',M,dabs(chem_up(M,L)-chem_upold(M,L))/chem_up(M,L)
 enddo	
  endif
  endif
      	if ((L.eq.Lmax).and.(KEYit.ne.0).and.(KEYit1.eq.1))	then
!	 KEYit1=0
	endif

	 !if (KEYit1.ne.0) KEYit=1
 		write (110,'(3X,I3,50(4X,E18.12))')L,TL(L),pN(L),Xloc(nk,L),(Xloc(jj,L),jj=1,NK-1)
	if (L.eq.1) 	write (111,'((4X,A14,4X),3X,A3,30(4X,A14,4X))')'coord','L','thickness','pN(L)',name_seg(NK),(name_seg(jj),jj=1,NK-1)!'water_Xloc',(jj,'_Xloc',jj=1,NK-1)
		write (111,'((4X,E18.12),3X,I3,50(4X,E18.12))')Rout(L),L,TL(L),pN(L),Xloc(nk,L),(Xloc(jj,L),jj=1,NK-1)
		if (mod(iter,10).eq.0) then
		if ((iter.gt.1).and.(L.eq.1)) OPEN(UNIT=113,FILE='current profiles.txt')
		if ((iter.gt.1).and.(L.eq.1)) 	write (113,'((4X,A14,4X),3X,A2,30(4X,A14,4X))')'coord','L','thickness',name_seg(NK),(name_seg(jj),jj=1,NK-1)!'water_Xloc',(jj,'_Xloc',jj=1,NK-1)
		write (113,'((4X,E18.12),3X,I3,50(4X,E18.12))')Rout(L),L,TL(L),Xloc(nk,L),(Xloc(jj,L),jj=1,NK-1)
		if (L.eq.Lmax)			write (113,'(A10,I5)')'iter=',iter
		if (L.eq.Lmax)			close (113)
		if ((iter.gt.1).and.(L.eq.1)) OPEN(UNIT=114,FILE='current chempot.txt')
		if ((iter.gt.1).and.(L.eq.1)) 	write (114,'(A3,30(4X,A18,4X))')'L',(name_mol(jj),jj=1,NKK)!'water_Xloc',(jj,'_Xloc',jj=1,NK-1)
		write (114,'(I3,50(4X,E18.12))')L,(total_mol(jj),jj=1,NKK)
		if (L.eq.Lmax)			write (114,'(A10,I5)')'iter=',iter
		if (L.eq.Lmax)			close (114)
		endif

if (L.eq.Lmax)		write (111,'(A10,E18.12)')'grand_pot=',grand_potential
!if (L.eq.Lmax)		write (111,'(A10,E18.12,E18.12)')'pressure=',-GAMloc(NK,1),-GAMloc(NK,Lmax)
!if (L.eq.Lmax)		write (111,'(A10,E18.12)')'r=',2.0d0*0.206404351d0/(-GAMloc(NK,1)+GAMloc(NK,Lmax))
if (L.eq.1)   write (35,'(A2,50(5X,I3,A9,6X))')'L',(jj,'.non_chem', jj=1,NG),(M,'.chembond', M=1,Nbonds)
if (L.eq.1) write (36,'(A2,50(5X,I3,A9,6X))')'L',(jj,'.non_chem', jj=1,NG),(M,'.chembond', M=1,Nbonds)
if (L.eq.1) write (37,'(A2,50(5X,I3,A9,6X))')'L',(jj,'.non_chem', jj=1,NG),(M,'.chembond', M=1,Nbonds)
if (L.eq.1) write (38,'(A2,50(5X,I3,A9,6X))')'L',(jj,'.non_chem', jj=1,NG),(M,'.chembond', M=1,Nbonds)

	write (110,'(I2,56(4X,E17.11))')L,(FU(ii),ii=1,Nx),FUtau
	write (110,'(I2,56(4X,E17.11))')L,(FU(ii),ii=1+Nx,Nx+Nalf+Nbonds)
    write (35,'(I2,56(4X,E18.12))')L,(AL_up(s,L),s=1,NG),(chem_up(k,L),k=1,Nbonds) !AL_up(4,L),AL_up(5,L)
	  write (36,'(I2,56(4X,E18.12))')L,(AL_dwn(s,L),s=1,NG),(chem_dwn(k,L),k=1,Nbonds) !,AL_dwn(4,L),AL_dwn(5,L)
	  write (37,'(I2,56(4X,E18.12))')L,(AL_al(s,L),s=1,NG),(chem_al(k,L),k=1,Nbonds) !,AL_al(4,L),AL_al(5,L)
	   write (38,'(I2,56(4X,E18.12))')L,(ALloc(s,L),s=1,NG),(teta_chem(k,L)/ZZ,k=1,Nbonds) !,ALloc(4,L),ALloc(5,L)	
!	   write (55,'(I2,16(4X,E12.6))')L,(A_upsum(jj,jj,L), jj=1,NG),A_upsum(1,2,L),A_upsum(1,3,L),A_upsum(2,3,L),(Ach_upsum(M,L), M=1,Nbonds)
!	   write (56,'(I2,16(4X,E12.6))')L,(A_alsum(jj,jj,L), jj=1,NG),A_alsum(1,2,L),A_alsum(1,3,L),A_alsum(2,3,L),(Ach_alsum(M,L), M=1,Nbonds)
	if ((L.eq.Lmax).and.(KEYit.eq.0).and.(final1.eq.0).and.(KEYit1.eq.0))	then
	 final1=1
	 KEYit1=1
	 KEYit=1
	 !goto 7
	endif

 	if ((L.eq.Lmax).and.(KEYit.eq.0).and.(final1.eq.1).and.(KEYit1.eq.0))	then
	 final1=0
	 KEYit1=1
	 KEYit=1
	 !goto 7
	endif

!switch=0
!!!for canonical
 if (((shape.eq.'sph').and.(L.eq.Lmax)).or.(((shape.eq.'pln').and.(L.eq.Lmax)))) then
 do k=1,NK-1
 X(k)=Xloc(k,Lmax+1) 
 enddo
! 	 call NEWTalf_chem(Nx,Nalf,X,PARalf,DPalf,FU,DD,Lmax+1,EPS,EPF,KEY,ghom)
 X(nk)=1.0d0
do k=1,NK-1
X(nk)=X(nk)-X(k)
enddo
call homogsol(ghom,ghom_formix)
 endif

if ((switch.eq.0).and.(L.eq.Lmax)) goto 88
!if ((switch).and.(L.eq.1)) goto 88

enddo			        	!cycle L
write (110,'(/)')
88 write(*,*)'iter=',iter
switch=0
			        
!!  if (iter.gt.97000) then		 !.and.(KEYit.eq.0)
!!  write (*,*)'Too many iterations in micelle_energy'
  !stop
!!  endif 
enddo						!cycle iter
if (final1.eq.0) then
final1=1
!goto 70
!if (shape.eq.'sph') goto 70
endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if (shape.eq.'pln') then
 !!shape='sph'
 !X(1)=1.0d0-0.906d0 !X2(1)*1.001d0 !0.91d0

 !X(nk)=1.0d0
!do k=1,NK-1
!X(nk)=X(nk)-X(k)
!enddo
!call homogsol(ghom,ghom_formix)
!!do L=1,Lmax
!!Xloc(2,L)=1.0d0-Xloc(2,L)
!!Xloc(1,L)=1.0d0-Xloc(1,L)
!!enddo
 !!Goto 70
 !!endif
 write (135,'(A2,30(5X,I2,A9,6X))')'L',(jj,'.non_chem', jj=1,NG),(M,'.chembond', M=1,Nbonds)
 write (136,'(A2,30(5X,I2,A9,6X))')'L',(jj,'.non_chem', jj=1,NG),(M,'.chembond', M=1,Nbonds)
 write (137,'(A2,30(5X,I2,A9,6X))')'L',(jj,'.non_chem', jj=1,NG),(M,'.chembond', M=1,Nbonds)
 write (138,'(A2,30(5X,I2,A9,6X))')'L',(jj,'.non_chem', jj=1,NG),(M,'.chembond', M=1,Nbonds)
  write (39,'(A2,30(8X,I2,A7,5X))')'L',(s,'.XX_up ',s=1,NG),(s,'.XX_dwn',s=1,NG),(s,'.XX_al ',s=1,NG)

 write (155,'(A2,30(2X,A8,f3.2,f3.2))')'L',('non_chem',0.01*jj,0.01*jj, jj=1,NG),(('non_chem',0.01*s,0.01*k,k=s+1,NG),s=1,NG-1),('seg_chem',0.01*seg_ch(M),0.01*seg_max_ch(M), M=1,Nbonds)
write (156,'(A2,30(2X,A8,f3.2,f3.2))')'L',('non_chem',0.01*jj,0.01*jj, jj=1,NG),(('non_chem',0.01*s,0.01*k,k=s+1,NG),s=1,NG-1),('seg_chem',0.01*seg_ch(M),0.01*seg_max_ch(M), M=1,Nbonds)
write (121,'((4X,A14,4X),3X,A2,30(4X,A14,4X))')'coord','L','thickness',name_seg(NK),(name_seg(jj),jj=1,NK-1) !,'water_Xloc',(jj,'_Xloc',jj=1,NK-1)
write (122,'(A3,30(4X,A14,4X))')'L',(name_mol(jj),jj=1,NKK)
write (123,'(A3,4X,A4,30(4X,A14,4X))')'L','comp','combinat','residual','chempart','transfer','pure','total'
55 tenss=0.0d0
 final=1
 call check_TL(TL,TLmin,Lmax)
  do L=1,Lmax
	 phase=0						  !  0 intermediate, -1 down, 1 up
	if (L.eq.1) phase=-1
	if (L.eq.Lmax) phase=1
	 do k=1,Nx
  PARx(k)=Xloc(k,L)						
  enddo
  do s=1,NG-1
  PARalf(NK-1+s)=Al_up(s,L)					
  PARalf(NK+NG-2+s+Nbonds)=Al_dwn(s,L)				
  enddo	 
   if (Nmol /= 0) then
      DO M=1,Nbonds
      PARalf(NK-1+NG-1+M)=chem_up(M,L)	
      enddo	
     endif

	DO  I=1,Nx+Nalf+Nbonds
	  if (I.le.Nx) then                                                        
      PAR(I)=PARx(I)
	  else
	  PAR(I)=PARalf(I)					
	  endif 
	  enddo
CALL inhomsol(Nx,Nalf,PAR,FU,L,FUtau,ghom)

!CALL free_energy(L,loc_energy,tot_energy)
!CALL grandpot_pressure(L,ghom,loc_tens,tension,grand_potential,delta)
  ! if (dabs(delta/TL(1)).ge.1.0d-1) goto 55
 	write (121,'((4X,E18.12),3X,I2,50(4X,E18.12))')Rout(L),L,TL(L),Xloc(nk,L),(Xloc(jj,L),jj=1,NK-1)
	write (122,'(I3,50(4X,E18.12))')L,(total_mol(jj),jj=1,NKK)
	do jj=1,NKK
	write (123,'(I3,4X,I3,50(4X,E18.12))')L,jj,comb_mol(jj),res_mol(jj),chempart_mol(jj),transfer_mol(jj),pure_mol(jj),total_mol(jj)
	enddo
	write (120,'(I2,50(4X,E17.11))')L,(FU(ii),ii=1,Nx),FUtau
	write (120,'(I2,50(4X,E17.11))')L,(FU(ii),ii=1+Nx,Nx+Nalf+Nbonds)
     write (135,'(I2,57(4X,E18.12))')L,(AL_up(s,L),s=1,NG),(chem_up(k,L),k=1,Nbonds)
	   write (136,'(I2,57(4X,E18.12))')L,(AL_dwn(s,L),s=1,NG),(chem_dwn(k,L),k=1,Nbonds) 
	   write (137,'(I2,57(4X,E18.12))')L,(AL_al(s,L),s=1,NG),(chem_al(k,L),k=1,Nbonds) 
	   write (39,'(I2,50(4X,E18.12))')L,(XX_up(s,L),s=1,NG),(XX_dwn(s,L),s=1,NG),(XX_al(s,L),s=1,NG)
	   write (138,'(I2,57(4X,E18.12))')L,(ALloc(s,L),s=1,NG),(teta_chem(k,L)/ZZ,k=1,Nbonds) 
	   write (155,'(I2,50(4X,E12.6))')L,(A_updwnsum(jj,jj,L), jj=1,NG),((A_updwnsum(s,k,L),k=s+1,NG),s=1,NG-1),(Ach_updwnsum(M,L), M=1,Nbonds)
	   write (156,'(I2,50(4X,E12.6))')L,(A_alsum(jj,jj,L), jj=1,NG),((A_alsum(s,k,L),k=s+1,NG),s=1,NG-1),(Ach_alsum(M,L), M=1,Nbonds)
	   write(220,'(I2,(4X,E12.6))')L,7.6d0*(AL_dwn(1,L)-AL_up(1,L))*m_up(L)/ALloc(1,L)*Xloc(nk,L)      
   enddo                          !  end of L-cycle
   	write (121,'((4X,E18.12),3X,I2,50(4X,E18.12))')Rout(Lmax+1),Lmax+1,TL(Lmax+1),Xloc(nk,Lmax+1),(Xloc(jj,Lmax+1),jj=1,NK-1)

  final=2
  do L=1,Lmax
!CALL free_energy(L,loc_energy,tot_energy)
Xlocc=Xloc
 pN(L)=Pnormal(L,Xlocc,Nx,Nalf,DPalf,EPS,EPF,ghom)
CALL grandpot_pressure(L,ghom,loc_tens,tension,grand_potential,delta)
enddo
!!   tenss=tenss+tot_energy !*Lmax
 !  print *,'tot_en=',tot_energy
print *,'tension=',tension
!!print *,'tension2=',tenss
	print *,delta,TL(1)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  	write (120,'(15(4X,E17.11))')((ghom(k)-ghom2(k)),k=1,NK)
!  write (120,'(A21,E18.12)')'grand_pot/V+p=',grand_potential
	write (120,'(A25,E18.12)')'grand_pot/V (from tang)=',energy2
	write (120,'(A10,E18.12)')'tension=',tension
return
end
!************************************************************************************************************
!SUBROUTINE free_energy(L,local,tot_energy)
real*8 function free_energy(L)
 IMPLICIT NONE
 real*8 local(150),tot_energy
 CHARACTER*3 shape
 character*7 name_seg(150)
integer L,Lmax,M,NK,NG,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150)
integer I,s,k
real*8 m_up,m_dwn,lat,TL(150),TLmin,XlocL(150),chemm(150),volume,GGA(150)   
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150)
REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150)
 common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al
common /layer/shape,TL,Lmax,TLmin
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
COMMON /C/T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
   chemm=chem
    free_energy=0.0d0
  DO k=1,NK
  XlocL(k)=Xloc(k,L)
  ENDDO
 call total(XlocL,QT,RT)
! print *,'before',L,XX_up(1:NG,L),XX_dwn(1:NG,L),XX_al(1:NG,L)
 do s=1,NG
if (L.eq.1) then
free_energy=free_energy+m_up(L)*AL_up(s,L)*dlog(AL_up(s,L)*XX_up(s,L)/ALloc(s,L))+lat(L)*AL_al(s,L)*dlog(AL_al(s,L)*XX_al(s,L)/ALloc(s,L))
!free_energy=free_energy-m_up(L)*ALloc(s,L)*dlog(ALloc(s,L))-lat(L)*ALloc(s,L)*dlog(ALloc(s,L))
else
free_energy=free_energy+m_up(L)*AL_up(s,L)*dlog(AL_up(s,L)*XX_up(s,L)/ALloc(s,L))+lat(L)*AL_al(s,L)*dlog(AL_al(s,L)*XX_al(s,L)/ALloc(s,L))+m_dwn(L)*AL_dwn(s,L)*dlog(AL_dwn(s,L)*XX_dwn(s,L)/ALloc(s,L))
!local(L)=local(L)-m_up(L)*ALloc(s,L)*dlog(ALloc(s,L))-lat(L)*ALloc(s,L)*dlog(ALloc(s,L))-m_dwn(L)*ALloc(s,L)*dlog(ALloc(s,L))
!free_energy=free_energy-ALloc(s,L)*dlog(ALloc(s,L))
endif
 enddo
!print *,'free1=',free_energy
 if (Nmol /= 0) then
 DO M=1,Nbonds
 if (L.eq.1) then
 free_energy=free_energy+m_up(L)*chem_up(M,L)*dlog(chem_up(M,L))+lat(L)*chem_al(M,L)*dlog(chem_al(M,L))
 else
 free_energy=free_energy+m_up(L)*chem_up(M,L)*dlog(chem_up(M,L))+lat(L)*chem_al(M,L)*dlog(chem_al(M,L))+m_dwn(L)*chem_dwn(M,L)*dlog(chem_dwn(M,L))
 endif
 enddo
 endif
 free_energy=free_energy*ZZ*QT
 call comb_part(R,Q,NK,RT,QT,ZR,GGA)		   !GGA=combinatorial part of chem pot 
 !print *,'free2=',free_energy

 DO k=1,NK
 free_energy=free_energy+Xloc(k,L)*(dlog(Xloc(k,L))+GGA(k))
 
 if (type_mol(k).gt.0) then
	free_energy=free_energy-chemm(k)*Xloc(k,L)*dlog(Xloc(k,L)/QT/ZZ) 	
endif	
 enddo
 !print *,local(L)
 !print *,'free3=',free_energy

! if (dabs(local(L)).lt.1.0d-14) local(L)=0.0d0 
free_energy=free_energy*volume(L)/RT 
!print *,'free4=',free_energy,RT

!tot_energy=0.0d0
!do i=1,Lmax
!tot_energy=tot_energy+local(i) !/Lmax
!enddo
 !write (210,'(I2,(4X,E18.12))')L ,local(L)
return
end
!************************************************************************************************************
SUBROUTINE grandpot_pressure(L,ghom,local_tens,tension,energy,delta)
 IMPLICIT NONE
 real*8 local(150),energy,local2(150),energy2,geom(150),local_tens(150),normal(150),tang(150),ghom(150),tension,tens,r_s,ii,pressure_zero,zero_temp,sumR,eps,delta,grandpot,grand_corr
 CHARACTER*3 shape
 character*7 name_seg(150)
integer L,Lmax,M,Mol,NK,NG,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150),seg_ch,seg_max_ch
integer I,s,k,iter,final,final1,surf,layer
real*8 m_up,m_dwn,lat,TL(150),kT,vol_cell,LL,LLmax,req,adsorb,T_all,volume,pi,Tcurrent,a,normall(150),TLmin   
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ,QThom,RThom,GGA(150),GGAhom(150),XlocL(150)
REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
REAL*8 AL_up(10,150),AL_dwn(10,150),AL_al(10,150),ALloc(10,150),ALchemold(150,150)
REAL*8 teta_chem(150,150),chem_up(150,150),chem_dwn(150,150),chem_al(150,150)
REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150)
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150),chemm(150),Rout,pN(150)
COMMON /H/ BLIS,AL,ALI,ET,GP 
 common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
common /chem_or/chem_up,ALchemold
common /chemdwn/chem_dwn,chem_al
common /layer/shape,TL,Lmax,TLmin
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
COMMON /H_or/ AL_up,AL_dwn,AL_al,ALloc
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds 
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
common /cycles/iter,final,final1
common /pressure2/normal,tang
common /grand/local,grandpot
common /energy/energy2
common /pressure/pN,dpT
 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10),dpT(150,150) 
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid   
chemm=chem
Pi = 3.1415926d0
!kT=1.380649d0*T !mN/m*Ang^2             !/101325.0d0	 !atm*A^3
kT=8.314462d0*T !mN/m*Ang^2
!normal(Lmax)=0.0d0
vol_cell=14.244d1 !/6.02214076d0 !27.0d0  !Ang^3
a=(vol_cell**0.6666666d0)*pi/6.0d0   !Ang^2
T_all=0.0d0
do layer=1,Lmax+1
T_all=T_all+TL(layer)
enddo
 !normal(Lmax)=0.0d0
  normal(L)=0.0d0
  local(L)=0.0d0
  DO k=1,NK
  XlocL(k)=Xloc(k,L)
  ENDDO
 call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
 call total(X,QThom,RThom) !,NK,Q,R,QT,RT)		!QT, RT - output

 do s=1,NG
 normal(L)= normal(L)+QT*ZZ/RT*(ALloc(s,L)*dlog(AL_al(s,L)*XX_al(s,L)/ALloc(s,L))) !+ (ALloc(s,L)-AL_al(s,L))) !/(1.0d0+AL_al(s,L)*XX_al(s,L)*XX_al(s,L)/AL_up(s,L)/XX_up(s,L)/XX_dwn(s,L+1)))
 if (shape.eq.'sph')  normal(L)= normal(L)+2.0d0/Rout(L)*AL_up(s,L)*dlog(AL_up(s,L)*XX_up(s,L)/AL_al(s,L)/XX_al(s,L))
! normal(Lmax)= normal(Lmax)+AL(s)*dlog(XX(s))
 if (L.eq.1) then
local(L)=local(L)+m_up(L)*AL_up(s,L)*dlog(AL_up(s,L)*XX_up(s,L)/XX(s))+lat(L)*AL_al(s,L)*dlog(AL_al(s,L)*XX_al(s,L)/XX(s))
local(L)=local(L)-m_up(L)*ALloc(s,L)*dlog(ALloc(s,L))-lat(L)*ALloc(s,L)*dlog(ALloc(s,L))

else
local(L)=local(L)+m_up(L)*AL_up(s,L)*dlog(AL_up(s,L)*XX_up(s,L)/XX(s))+lat(L)*AL_al(s,L)*dlog(AL_al(s,L)*XX_al(s,L)/XX(s))+m_dwn(L)*AL_dwn(s,L)*dlog(AL_dwn(s,L)*XX_dwn(s,L)/XX(s))
!local(L)=local(L)-m_up(L)*ALloc(s,L)*dlog(ALloc(s,L))-lat(L)*ALloc(s,L)*dlog(ALloc(s,L))-m_dwn(L)*ALloc(s,L)*dlog(ALloc(s,L))
local(L)=local(L)-ALloc(s,L)*dlog(ALloc(s,L))
endif

 enddo
!print *,local(L)

 if (Nmol /= 0) then
 DO M=1,Nbonds
  normal(L)= normal(L)+Xloc(seg_ch(M),L)*dlog(chem_al(M,L))+0.5d0*(Xloc(seg_ch(M),L)-Xloc(seg_max_ch(M),L))*dlog(AL_up(NG,L)*XX_up(NG,L)*AL_dwn(NG,L+1)*XX_dwn(NG,L+1)*chem_al(M,L)/AL_al(NG,L)/XX_al(NG,L)/AL_al(NG,L+1)/XX_al(NG,L+1)/chem_dwn(M,L+1))
 if (L.eq.1) then
 local(L)=local(L)+m_up(L)*chem_up(M,L)*dlog(chem_up(M,L)*(X(seg_ch(M))/QThom/ZZ))+lat(L)*chem_al(M,L)*dlog(chem_al(M,L)*(X(seg_ch(M))/QThom/ZZ))
else
 local(L)=local(L)+m_up(L)*chem_up(M,L)*dlog(chem_up(M,L)*(X(seg_ch(M))/QThom/ZZ))+lat(L)*chem_al(M,L)*dlog(chem_al(M,L)*(X(seg_ch(M))/QThom/ZZ))+m_dwn(L)*chem_dwn(M,L)*dlog(chem_dwn(M,L)*(X(seg_ch(M))/QThom/ZZ))
 endif
 enddo
 endif
!print *,local(L)
 local(L)=local(L)*ZZ*QT
 call comb_part(R,Q,NK,RT,QT,ZR,GGA)		   !GGA=combinatorial part of chem pot 
 call comb_part(R,Q,NK,RThom,QThom,ZR,GGAhom)		   !GGA=combinatorial part of chem pot 
 DO k=1,NK
 local(L)=local(L)+Xloc(k,L)*(dlog(Xloc(k,L)/X(k))+GGA(k)-GGAhom(k))
 normal(L)= normal(L)+Xloc(k,L)/RT*(dlog(Xloc(k,L))+GGA(k))
! normal(Lmax)= normal(Lmax)+X(k)/QThom/ZZ*(dlog(X(k))+GGAhom(k))
	if (type_mol(k).ne.0) then
	 normal(L)= normal(L)-chemm(k)/RT*Xloc(k,L)*dlog(Xloc(k,L)/QT/ZZ)
	 !normal(Lmax)= normal(Lmax)-chemm(k)/2.0d0*X(k)/QThom/ZZ*dlog(X(k)/QThom/ZZ)
	 	local(L)=local(L)-chemm(k)*Xloc(k,L)*dlog(Xloc(k,L)/QT/ZZ) 
!	if (chem(k).eq.1)	local(L)=local(L)-Xloc(k,L)*dlog(Xloc(k,L)/QT/ZZ) !/QT !m_up(L)*Xloc(k,L)/QT/ZZ*dlog(Xloc(k,L)/QT/ZZ)-lat(L)*Xloc(k,L)/QT/ZZ*dlog(Xloc(k,L)/QT/ZZ)-m_dwn(L)*Xloc(k,L)/QT/ZZ*dlog(Xloc(k,L)/QT/ZZ)
!	if (chem(k).eq.2)	local(L)=local(L)-2.0d0*Xloc(k,L)*dlog(Xloc(k,L)/QT/ZZ)	  !/QT
!	if (chem(k).eq.3)	local(L)=local(L)-3.0d0*Xloc(k,L)*dlog(Xloc(k,L)/QT/ZZ)	   !/QT
!	if (chem(k).eq.4)	local(L)=local(L)-4.0d0*Xloc(k,L)*dlog(Xloc(k,L)/QT/ZZ) 	   !/QT
	endif
normal(L)= normal(L)-(ghom(k)-GP(k))*Xloc(k,L)/RT
 enddo
!
!print *,-normal(L),ghom(NK-1)/R(NK-1)-ghom2(NK-1)/R(NK-1)
 normal(L)= -ghom(NK-1)/R(NK-1)-normal(L) !(-normal(L))/RT !-ghom(NK-1)
! DO Mol=1,Nmol
! k=1
!	do while (type_mol(k).ne.Mol)
!	k=k+1
!	enddo
!	DO while ((type_mol(k).eq.Mol).and.(rank_seg(k+1).gt.1)) .or. ((type_mol(k).eq.Mol) .and. (type_branch(k) .ne. 0) .and. (
!		local(L)=local(L)+(Xloc(k,L)/QT)*dlog(X(k)/QThom/ZZ) !m_up(L)*chem_up(M,L)*dlog(X(k)/QThom/ZZ)+lat(L)*chem_al(M,L)*dlog(X(k)/QThom/ZZ)+m_dwn(L)*chem_dwn(M,L)*dlog(X(k)/QThom/ZZ)
!	k=k+1
!	enddo		
 !enddo
 LL=L
 LLmax=Lmax
 if (shape.eq.'sph') then
 geom(L)=volume(L)/(4.0d0*pi/3.0d0*T_all**3) !(3.0d0*LL**2-3.0d0*LL+1.0d0)/(1.0d0*LLmax**3)
 elseif (shape.eq.'cyl') then
 geom(L)=volume(L)/(pi*T_all**2) !(2.0d0*LL-1.0d0)/LLmax**2
 else
 geom(L)=volume(L)/T_all !1.0d0/LLmax
 endif

 local_tens(L)=local(L)/RT !+ghom(NK) !*kT/vol_cell		 !MPa

 local(L)=local(L)*geom(L)/RT

energy=0.0d0
do i=1,Lmax
energy=energy+local(i)
enddo
normal=pN
!normal(Lmax)=-ghom(NK)/R(NK)
!print *,'!!!',-ghom(NK),-ghom(NK-1)
tang(L)=-local_tens(L) !*QT/RT
!if (L.eq.1) print *,'pT(1)=',tang(L)
if (L.eq.Lmax) then

	do i=1,Lmax
	!normal(i)=-ghom(NK-1)
if (final.eq.2)	tang(i)=-local_tens(i)+normal(Lmax) !-ghom(NK)/R(NK) !+normal(Lmax) !-GAMloc(NK,i)
	enddo

!	endif	!shape
 energy2=0.0d0
do i=1,Lmax
energy2=energy2-tang(i)*geom(i) !/geom(Lmax) !*vol_cell/kT
enddo

!	normal_zero=normal(1)
!tang(L)=normal(L)-local_tens(L)
	if (final.eq.2) then

		if (shape.eq.'sph') then
		if (r_solid.lt.1.0d-50)	call surf_of_tension(delta,r_s,tens,grand_corr,normall)
			write (200,'(5(4X,A17),8X,A2,5(4X,A17))')'coord','local_energy','Pn(ln)','Pn','Pt','L','Pn[MPa]','Pt[MPa]'
		tension=tens !*3.0d-4				 !N/m
		grandpot=grand_corr*4.0d0*pi
		if (r_solid.gt.1.0d-50)  delta=r_solid 
		else	  !plane
			write (200,'(4(4X,A17),8X,A2,5(4X,A17))')'coord','local_energy','Pn','Pt','L','Pn[MPa]','Pt[MPa]'
		tension=energy*T_all !LLmax !*kT/vol_cell*3.0d-4 !*3.0d-7*1.0d9 ! N/m   !101325.0d0  !for plane
		grandpot=tension
		endif
	 call adsorption(delta,r_s,req,adsorb)
	 Tcurrent=TL(1) !0.0d0
	do i=1,Lmax
	Tcurrent=Tcurrent-TL(i)
	if (shape.eq.'sph') then
		if ((i.ne.Lmax).or.(r_solid.gt.1.0d-50)) write (200,'(5(4X,E18.12),5X,I2,5(4X,E18.12))')Tcurrent+delta, local_tens(i), normall(i),normal(i),tang(i),i,normal(i)*kT*1.0d1/vol_cell,tang(i)*kT*1.0d1/vol_cell !delta/TL(1)+i, normall(i),normal(i),tang(i),delta+Tcurrent,tang(i)*kT*1.0d1/vol_cell
	!	if ((i.eq.Lmax).and.(delta.lt.0.5d0)) write (200,'(f7.3,5(4X,E18.12))')-delta+i, local_tens(i),normal(i),tang(i)
		if (i.eq.Lmax) write (200,'((8X,A7,7X),(3X,A17,2X),(3X,A21,2X),(7X,A6,12X),(2X,A10,9X),(2X,A10,2X))')'tension', 'radius_of_tension','radius_equimol(water)','tolmen','adsorbtion','grandpot' !, 'pressure_in_center'
		if (i.eq.Lmax) write (200,'(6(4X,E18.12))')tension,r_s,req,req-r_s,adsorb,grandpot !,pressure_zero	
	
	!	if (dabs(delta/TL(1)).gt.1.0d0) stop
	else	  !plane
		write (200,'(4(4X,E18.12),5X,I2,4(4X,E18.12))')Tcurrent, local_tens(i),normal(i),tang(i),i,normal(i)*kT*1.0d1/vol_cell,tang(i)*kT*1.0d1/vol_cell
		if (i.eq.Lmax) write (200,'((8X,A7,7X),(3X,A21,2X),(3X,A15,2X),(2X,A18,2X),(2X,A10,2X))')'tension','radius_equimol(water)','adsorption','grandpot'
		if (i.eq.Lmax) write (200,'(5(4X,E18.12))')tension,req,adsorb,grandpot,-ghom(NK-1)/R(NK-1)*kT*1.0d1/vol_cell !,-ghom(NK-1)
	endif 
	enddo
	if (shape.eq.'sph') then
 	 write (200,'((4X,E18.12),A5,61X,(4X,E18.12),A9)')tension*kT/a/(6.02214076d0**0.3333333333d0),'mN/m',adsorb/a/6.02d3,'mole/m^2'
	else
	 write (200,'((4X,E18.12),A5,17X,(4X,E18.12),A9)')tension*kT/a/(6.02214076d0**0.3333333333d0),'mN/m',adsorb/a/6.02d3,'mole/m^2'
	 endif
	endif  !final
endif	 !L=Lmax
 !	TL(1)=TL(1)+delta
	
! write (200,'(12(4X,E18.12))')local(L)
!print *,'end'
return
end
!************************************************************************************************************
 SUBROUTINE surf_of_tension(delta,rs,tension,energy,normall)
 IMPLICIT NONE
  CHARACTER*3 shape
  COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
 common /pressure2/normal,tang
 real*8 tang(150),normal(150),normall(150),tension,rs,delta,sum0,sum_tang,LL,delta1,delta2,energy,Tcurrent,TL(150),TLmin,QT,RT,XlocL(150)
 REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
integer L,Lmax
common /layer/shape,TL,Lmax,TLmin
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
  REAL*8 ZZ,ZR,JX,IJ
integer NK,NG,k,Nmol,Nbonds

Tcurrent=0.0d0
delta1=0.0d0
delta2=0.0d0
sum_tang=0.0d0
sum0=0.0d0
do L=1,Lmax
Tcurrent=Tcurrent+TL(L)
LL=L
sum_tang=sum_tang+(tang(L)-normal(Lmax))*TL(L)/(tang(1)-normal(Lmax))
sum0=sum0+(tang(L)-normal(Lmax))/(tang(1)-normal(Lmax))*TL(L)*(Tcurrent-0.5d0*TL(L))
enddo

if ((sum_tang**2-2.0d0*sum0).gt.0.0d0) then
delta1=-sum_tang+dsqrt(sum_tang**2-2.0d0*sum0)
delta2=-sum_tang-dsqrt(sum_tang**2-2.0d0*sum0)
	if (( (sum0*delta1) .lt.0.0d0).and.( (sum0*delta2) .lt.0.0d0)) then
	print *,'2 delta for the first layer',delta1,delta2
	delta=min(dabs(delta1),dabs(delta2))
	if (sum0.gt.0.0d0) delta=-delta
	elseif ((sum0*delta1) .lt.0.0d0) then
	delta=delta1
	elseif ((sum0*delta2) .lt.0.0d0) then
	delta=delta2
	else
	print *,'2 wrong delta for the first layer',delta1,delta2
	delta=min(dabs(delta1),dabs(delta2))
	if (delta1.lt.0.0d0) delta=-delta
	endif
else
print *,'(sum_tang**2-2.0d0*sum0)<0'
delta=0.0d0
endif
if (delta.lt.(-TL(1))) print *,'the first layer is very small',TL(1)+delta

energy=0.0d0
Tcurrent=0.0d0
do L=1,Lmax
 DO k=1,NK
  XlocL(k)=Xloc(k,L)
  ENDDO
 call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output

Tcurrent=Tcurrent+TL(L)
LL=L
if (L.eq.1) energy=energy-(tang(L)-normal(Lmax))*RT*((Tcurrent+delta)**3)/3.0d0
if ((L.eq.Lmax).and.(delta.lt.0.5d0))	energy=energy-(tang(L)-normal(Lmax))*RT*((Tcurrent-delta)**3-(Tcurrent-TL(L)+delta)**3)/3.0d0
if ((L.ne.1).and.(L.ne.Lmax)) energy=energy-(tang(L)-normal(Lmax))*RT*((Tcurrent+delta)**3-(Tcurrent-TL(L)+delta)**3)/3.0d0 !*(3.0d0*(LL+delta)**2-3.0d0*(LL+delta)+1.0d0)/3.0d0

if (L.eq.1) normall(1)=tang(1)
if (L.ne.1)	normall(L)=normall(L-1)+2.0d0*(tang(L-1)-normall(L-1))*dlog(Tcurrent/(Tcurrent-TL(L)))
enddo
if ((energy/tang(1)).gt.0.0d0) then
rs=(6.0d0*energy/(tang(1)-normal(Lmax)))**(0.3333333333333333d0)			 !3*4pi/2pi=6
tension=3.0d0*energy/(rs**2)
else
print *,'expression for rs<0'
endif

return
end
!************************************************************************************************************
  SUBROUTINE adsorption(delta,rs,req,adsorb)
 IMPLICIT NONE
 CHARACTER*3 shape
 character*7 name_seg(150)
 real*8 Rout,req,TL(150),LL,delta,rs,sum_NK,sum_comp,adsorb,bulkL1,bulkLmax,RT,QT,XlocL(150),Tcurrent
 REAL*8 Xloc(150,150),XX_up(10,150),XX_dwn(10,150),XX_al(10,150),GAMloc(150,150),Xchem(150)
 REAL*8 ZZ,ZR,JX,IJ,teta_chem(150,150),segments,TLmin
integer NK,NG,k,L,Lmax,Nmol,chem(150),type_mol(150),Nbonds,rank_seg(150)
common /layer/shape,TL,Lmax,TLmin
COMMON /S_or/ Xloc,XX_up,XX_dwn,XX_al,GAMloc,Xchem
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds
 REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid 
segments=0.0d0
sum_NK=0.0d0
sum_comp=0.0d0
bulkL1=0.0d0
bulkLmax=0.0d0
do k=1,NK
		if ( (type_mol(k).eq.1).or.((Nmol.eq.0).and.(k.eq.1)) ) then
		bulkL1=bulkL1+Xloc(k,1)
		bulkLmax=bulkLmax+Xloc(k,Lmax)
		if  (type_mol(k).eq.1) then
			 segments=segments+1.0d0
		else
			 segments=1.0d0
		endif
		endif
enddo

if (shape.eq.'sph') then
	DO k=1,NK
		 XlocL(k)=Xloc(k,1)
		ENDDO
		call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
Tcurrent=TL(1)
sum_NK=sum_NK+(Xloc(NK,1)-Xloc(NK,Lmax))*((TL(1)+delta)**3-r_solid**3)/RT
		do k=1,NK
		if ( (type_mol(k).eq.1).or.((Nmol.eq.0).and.(k.eq.1)) ) then
		sum_comp=sum_comp+(Xloc(k,1)-Xloc(k,Lmax))*((TL(1)+delta)**3-r_solid**3)/RT
		endif
		enddo
	do L=2,Lmax
	Tcurrent=Tcurrent+TL(L)
	LL=L
	sum_NK=sum_NK+(Xloc(NK,L)-Xloc(NK,Lmax))*((Tcurrent+delta)**3-(Tcurrent-TL(L)+delta)**3)/RT
		do k=1,NK
		if ( (type_mol(k).eq.1).or.((Nmol.eq.0).and.(k.eq.1)) ) then
		sum_comp=sum_comp+(Xloc(k,L)-Xloc(k,Lmax))*((Tcurrent+delta)**3-(Tcurrent-TL(L)+delta)**3)/RT
		endif
		enddo
	enddo
LL=Lmax
if (r_solid.gt.1.0d-50) then  !YYY
	 adsorb=sum_comp/(3.0d0*segments*r_solid**2)
	  else

	if (dabs(Xloc(NK,1)-Xloc(NK,Lmax)).le.1.0d-50) then
	print *,'Requimol does not exist, for solvent Xbulk1=Xbulk2'
if (rs.gt.0.0d0) req=rs
	elseif ((  sum_NK/(Xloc(NK,1)-Xloc(NK,Lmax))  ).le.0.0d0) then
	print *,'Requimol <0 for solvent'
if (rs.gt.0.0d0) req=rs
	else
req=sum_NK/(Xloc(NK,1)-Xloc(NK,Lmax)) !(  (sum_NK-Xloc(NK,Lmax)*(LL+delta)**3)/(Xloc(NK,1)-Xloc(NK,Lmax))  )**0.3333333333d0
adsorb=(sum_comp+(bulkLmax-bulkL1)*req)/segments/3.0d0

LL=0.0d0
Tcurrent=0.0d0
L=0
do while (LL.lt.req)
L=L+1
 		DO k=1,NK
		 XlocL(k)=Xloc(k,L)
		ENDDO
		call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
Tcurrent=Tcurrent+TL(L)
if (L.eq.1)  LL=LL+((Tcurrent+delta)**3)/RT
if (L.ne.1) LL=LL+((Tcurrent+delta)**3-(Tcurrent-TL(L)+delta)**3)/RT
enddo
!req=TL(L)*( (LL-req)/((Tcurrent+delta)**3-(Tcurrent-TL(L)+delta)**3)*RT)**0.33333333d0+Tcurrent-TL(L)+delta  !in units of tau standart
req=((Tcurrent+delta)**3+(req-LL)*RT)**0.33333333d0
adsorb=adsorb/(req**2)
	endif
	endif
!adsorb=sum_comp+(bulkLmax-bulkL1)*req	!(sum_comp+(bulkLmax-bulkL1)*req**3-bulkLmax*(LL+delta)**3)/(3.0d0*segments*req**2)
else	   !plane
	do L=1,Lmax
		DO k=1,NK
		 XlocL(k)=Xloc(k,L)
		ENDDO
		call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
	sum_NK=sum_NK+(Xloc(NK,L)-Xloc(NK,Lmax))*TL(L)/RT
		do k=1,NK
		if ( (type_mol(k).eq.1).or.((Nmol.eq.0).and.(k.eq.1)) ) then
		sum_comp=sum_comp+(Xloc(k,L)-Xloc(k,Lmax))*TL(L)/RT
		endif
		enddo
	enddo  !L
LL=Lmax
if (r_solid.gt.1.0d-50) then  !YYY
	 adsorb=sum_comp/segments
	  else
	if (dabs(Xloc(NK,1)-Xloc(NK,Lmax)).le.1.0d-50) then
	print *,'Requimol does not exist, for solvent Xbulk1=Xbulk2'
	if (dabs(bulkLmax-bulkL1).le.1.0d-50) adsorb=sum_comp/segments !adsorb=	(sum_comp-bulkLmax*LL)/segments
	else
req=sum_NK/(Xloc(NK,1)-Xloc(NK,Lmax))  !(sum_NK-Xloc(NK,Lmax)*LL)/(Xloc(NK,1)-Xloc(NK,Lmax)) 
adsorb=	(sum_comp+(bulkLmax-bulkL1)*req)/segments !(sum_comp+(bulkLmax-bulkL1)*req-bulkLmax*LL)/segments
LL=0.0d0
Tcurrent=0.0d0
L=0
do while (LL.lt.req)
L=L+1
 		DO k=1,NK
		 XlocL(k)=Xloc(k,L)
		ENDDO
		call total(XlocL,QT,RT) !,NK,Q,R,QT,RT)		!QT, RT - output
Tcurrent=Tcurrent+TL(L)
LL=LL+TL(L)/RT
enddo
req=Tcurrent+(req-LL)*RT !TL(L)*((LL-req)/TL(L)*RT)+Tcurrent-TL(L)  !in units of tau standart
	endif
 endif
endif


 
return
end
!************************************************************************************************************

