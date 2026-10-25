subroutine indsurf(r_max) 
IMPLICIT NONE    
 	  REAL*8 T2,Pi,psi,k,e,Nav,r_max
      COMMON /A/ T2,Pi,psi,k,e,Nav
	REAL*8 X(150),X1(150),X2(150),XX(150),Y(150),GAM1(150),GAM2(150)
	integer NK,NG,Nmol,type_mol(150),num_seg(150),Nfix,chem(150),Nbonds,rank_seg(150) 
    integer M,N,I,J,s,ss,D,No_or(150),CHANGE(150),Nsolve,Nx,Nalf, rank_branch(150), n_branch, type_branch(150), type_branch2(150), seg_max_ch, seg_ch
	integer type_group(10)
	character*7 name_seg(150)
REAL*8 W(10,10),WT(10,10),WTT(10,10),REDT(10,10),chemm(150),comp(150,10),RL1(30),RL2(30),IG(30,30),IGT(30,30),IGTT(30,30)
REAL*8 QT,RT,T,TT,R(150),RS(150,10),Q(150),QS(150,10),vol_ref, q_ref
REAL*8 ZZ,ZR,JX,IJ,DPx(150),DPalf(150),EPS,EPF
REAL*8 BLIS(150,10),AL(10),ALI(10),ET(10,10),GP(150),teta_chem(150,150)
REAL*8 r_solid, ET_solid(10), W_solid(10), WT_solid(10),WTT_solid(10)
common /chem/chem,type_mol,name_seg,teta_chem,rank_seg
common /branch/rank_branch, n_branch, type_branch, type_branch2
COMMON /S/ X,X1,X2,XX,Y,GAM1,GAM2
COMMON /B/ W,WT,WTT,REDT                                 
COMMON /C/ T,TT,R,RS,Q,QS,vol_ref, q_ref
COMMON /D/ ZZ,ZR,NK,NG,JX,IJ,Nmol,Nbonds                                      
COMMON /H/ BLIS,AL,ALI,ET,GP                                           
common /without_or/ Nfix,No_or,CHANGE
common /forNEWT/Nx,Nalf,DPx,DPalf,EPS,EPF
COMMON /solid/ r_solid, ET_solid, W_solid, WT_solid,WTT_solid
      OPEN(UNIT=31,FILE='FILEgeom.FT')
      OPEN(UNIT=32,FILE='FILEenergy.FT')
	  OPEN(UNIT=5,FILE='input.dat')
	  	507  FORMAT (25X)
	   508  FORMAT (100X)
	   20 FORMAT (2(F13.5))
	   21 FORMAT (3(E15.8,1X))
	   22  FORMAT (3I4)
	   23 FORMAT (5X,20I4)
	   24 FORMAT (A8,20F4.1)
	   25 FORMAT  (4X,5I4)
	   do M=1,30
		READ (31,20)RL1(M),RL2(M)
			DO N=M,30 
			READ (32,21) IG(M,N),IGT(M,N),IGTT(M,N)
			IG(N,M)=IG(M,N)
			IGT(N,M)=IGT(M,N)
			IGTT(N,M)=IGTT(M,N)
			enddo
		enddo
	   READ (5,507)
	   READ (5,22) NK,NG,Nmol
	   !print *,NK,NG,Nmol 
	   READ (5,507)
	   READ (5,23)(type_group(J),J=1,NG)
	  ! print *,(type_group(J),J=1,NG)
	   READ (5,507)
	   do I=1,NK	
	   READ (5,24)name_seg(I),(COMP(I,J),J=1,NG)
	  ! print *,name_seg(I),(COMP(I,J),J=1,NG)
	   enddo
	   READ (5,507)
	   do I=1,NK
	   READ (5,25)chem(i),type_mol(i),rank_seg(i),type_branch(i),rank_branch(i)
	  ! print *,chem(i),type_mol(i),rank_seg(i),type_branch(i),rank_branch(i)
	   enddo
	TT=298.15d0	 !standart for param(T,TT)
	T=T2
	!Nmol=1   			!number of chains (molecules with chem bonds)
	!NK=11   !14     			!number of segments
	!NG=8 !4    ! 5   			! number of different nonchem groups 
	ZZ=10.0d0
	ZR=0.5d0*ZZ

    do I=1,NG
	W_solid(I)=0.0d0
	WT_solid(I)=0.0d0
	WTT_solid(I)=0.0d0
	do J=1,NG
	W(I,J)=0.0d0
	WT(I,J)=0.0d0
	WTT(I,J)=0.0d0
	end do
	end do
	W_solid(2)=-2.0d0 !-4.988d0   ! nucleophilic	  with H-water
	W_solid(4)=-2.0d0 !-5.744d0   ! nucleophilic	with H-BuoH
	W_solid(3)=0.2d0	! hydrophilic:  unfavorable  with CH3  BuOH
	 W_solid(6)=0.2d0	 !	 with CH2 BuOH
!	 W_solid(1)=-4.988d0 
!	 W_solid(5)=-5.744d0 
!===========================================================================================================
!  SETTING THE GEOMETRY PARAMETERS OF THE GROUPS:
	do 1300 I=1,NK
	 !  rank_branch(I) = 0
	 !  type_branch(I) = 0
	   type_branch2(I) = 0
	do 1300 J=1,NG
	    W(J,J)=0.0d0
	    RS(I,J)=0.0d0
 1300	BLIS(I,J)=1.0d0
!RS(comp,group) - volume
!BLIS(comp,group) - bulkiness factor
!chem(comp) - number of chemical bonds on comp
!type_mol(comp) - index of chain. =0 (segment without chem bonds), =1 or other (segment with chem bonds)
!rank_seg(comp) - ranking number of segment along the main chain. =0 (segment without chem bonds or segment of branch), =1 or other (segment of the main chain) 
!type_branch(comp) - number of segment from which the branching (with segment=comp) starts, =0 (segment=comp is not part of branching)
!rank_branch(comp)- ranking number of segment along the branching, =0 (segment=comp is not part of branching) 
!type_branch2(comp) - number of segment from which the second branching (with segment=comp) starts, =0 (segment=comp is not part of branching)

 	do  I=1,NK
	do  J=1,NG
	if (type_group(J).ne.0) RS(I,J)=COMP(I,J)*RL1(type_group(J))                                                                                             
	if (type_group(J).ne.0) BLIS(I,J)=COMP(I,J)*(RL2(type_group(J))-1.0d0)+1.0d0
	enddo
	write(5,*)'No	RS for groups'
	write (5,'(I4,20(E15.8,1X))')I,(RS(I,J),J=1,NG)
	write(5,*)'No	LS for groups'
	write (5,'(I4,20(E15.8,1X))')I,(BLIS(I,J),J=1,NG)
	enddo                                 
	do  I=1,NG
	do  J=I+1,NG
	if ((type_group(J).ne.0).and.(type_group(I).ne.0)) then
	WT(I,J)=IGT(type_group(I),type_group(J))
	WTT(I,J)=IGTT(type_group(I),type_group(J))
	W(I,J)=IG(type_group(I),type_group(J))
	endif
	enddo
	enddo

	n_branch = 0
	do I=1,NK
	if (chem(I) .eq. 3) then
	n_branch = n_branch + 1
	endif
	if (chem(I) .eq. 4) then
	n_branch = n_branch+2
	endif
	enddo

	do I=1,NK
	num_seg(I)=1
	enddo

   Nbonds=0
   do I=1,NK
   Nbonds=Nbonds+chem(I)
   enddo
   Nbonds=Nbonds/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  SETTING THE INTERACTION PARAMETERS OF THE GROUPS:			               (for an athermal solution  all W(I,J)=0.0)

		write(5,*)'s	Wst for groups'
	do  I=1,NG-1
	do  J=I+1,NG
		WT(J,I)=WT(I,J)
		WTT(J,I)=WTT(I,J)
 	W(J,I)=W(I,J)
	enddo
		write (5,'(I4,20(E15.8,1X))')I,(W(I,J),J=1,NG)
	enddo
		close (5)
	call param(T,TT)
	 r_max=0.0d0
      DO  I=1,NK
	R(I)=0.0d0
	Q(I)=0.0d0
	GP(I)=0.0d0
      DO  J=1,NG
	QS(I,J)=RS(I,J)-2.0d0*(BLIS(I,J)+RS(I,J)-1.0d0)/ZZ
	R(I)=R(I)+RS(I,J)
  Q(I)=Q(I)+QS(I,J)
  enddo
  if (R(I).gt.r_max) r_max=R(I)
  enddo


	DO  I=1,NK
	chemm(I)=chem(I)
	Q(I)=Q(I)+chemm(I)/ZZ
	ENDDO
  
  Nfix=0	
	DO J=1,NG-1
	ss=1
	DO  I=1,NK
	if ((Q(I)-QS(I,J)).le.1.0d-150) then
		do s=1,NK
		if ( (s.ne.I).and.(QS(s,J).gt.1.0d-150).and.((Q(s)-QS(s,J)).gt.1.0d-150) ) ss=0
		enddo
	!if (ss)	then
	if (ss /= 0) then
!	Nfix=Nfix+1
!	No_or(Nfix)=J		!number of group, Al_up=Al_dwn=Al_al=Al_loc always (monomer with one group in monolayer)
	endif

	endif
	ENDDO !I
	ENDDO  !J

	Nx=NK-1            
	Nalf=(NG-1)*2
	Nsolve=Nalf-Nfix*2
	 !if (Nsolve+Nbonds) then
	if ((Nsolve+Nbonds) /= 0) then

	SS=0 
	DO D=1,Nalf/2
	 	N=0
		do J=1,Nfix
	   if (D.eq.No_or(J)) N=1
	   enddo
	   if (N.eq.0) then
	   SS=SS+1 
	   CHANGE(SS)=Nx+D			    !up nonchem
	   endif
	enddo
	   !if (Nmol) then
		if (Nmol /= 0) then
	DO D=1+Nalf/2,Nalf/2+Nbonds
	   SS=SS+1
	   CHANGE(SS)=Nx+D				!up chem
	enddo
	   endif
	
	 DO D=1+Nalf/2+Nbonds,Nalf+Nbonds
		N=0
		do J=1,Nfix
	   if (D.eq.(No_or(J)+Nalf/2+Nbonds)) N=1
	   enddo
	   if (N.eq.0) then 
	    SS=SS+1
	   CHANGE(SS)=Nx+D			    !down nonchem
	   endif
	enddo
	endif
!	print *,'Rnk=',R(NK)
	return
	end
!************************************************************************************************************************************
