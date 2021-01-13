!##########################################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Mesh coarsening                                                         //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: January, 07, 2017                                                              //!
!// Developed by: M. Hashemi, Iran, Tehran, Mohammadhashemi@ut.ac.ir                     //!
!// Doc ID: MC2F101F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 program Coarsening
 implicit none
!===============================
 Integer,Parameter::Dim=120000
!===============================
 Integer::NC,NP,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,NR,NS
 Integer::I,J,K,a,bp1,bp2,bp3,bp4,p1,p2,P,MIS,COARSE,n,m,counter,q,counter1,counter2,counter3
 Integer,Dimension(1:100)::NFR,BC
 Real(8)::dx1,dx2,dy1,dy2,d1,d2,teta,C,subx,suby,pi,Beta,Teta0
 Real(8),Dimension(1:Dim)::X,Y,NNp,F0,F0_C,f1,DelRow,NICM,point1,point2,cp1,cp2,coarsepoint,bs,NOCP
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:10000,1:10000)::d
 real :: start, finish
!*********************************************************
 call cpu_time(start)
    
 !part 1:
 Open(1,File='Mesh.txt')
 Open(2,File='Eliminating nodes.txt')
     
 !part 2:
 C=2
 Beta=3.4
 Teta0=5
 
 !part 3:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
 
 !part 4:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
 !part 5:
 
 !part 6:
 Do I=1,NP
  NNp(I)=1000;
  f1(I)=1000; 
 End Do 
 
 !part 7:
 Do I=1,NP
  Do J=1,NP
   
    If (I==J) Then
   Cycle
    End If 
    
    !part 7.1:
    subx=X(I)-X(J)
    suby=Y(I)-Y(J)
    d(I,J)=sqrt(subx*subx+suby*suby)
    
    !part 7.2:
    If (d(I,J)<=NNp(I)) Then
       NNp(I)=d(I,J)
    End If
  End Do
  
  !part 8:
  F0(I)=0.5*NNp(I)
 End Do  
 
 !part 9:
 Do I=1,NP
  Do J=1,NP
  
   If (I==J) Then
    Cycle
   End If
   
   F0_C(I)=C*F0(J)+d(I,J)
   
   !part 9.1:
   If (F0_C(I)<=F1(I)) Then
    F1(I)=F0_C(I)
   End If
  End Do  
 End Do   
   
 a=0
 counter1=0;
 pi=3.141592654
 
 !part 10:
 Do I=NF2+1,NF
  bp1=IDS(3,I)
  bp2=IDS(4,I)
  
  Do J=I,NF
   bp3=IDS(3,J)
   bp4=IDS(4,J)
   
   !part 10.1:
   If (I/=J .AND. (bp1==bp3 .OR. bp1==bp4)) Then
    
    dx1=x(bp2)-x(bp1)
    dy1=y(bp2)-y(bp1)
    dx2=x(bp4)-x(bp3)
    dy2=y(bp4)-y(bp3)
    d1=sqrt(dx1*dx1+dy1*dy1)
    d2=sqrt(dx2*dx2+dy2*dy2)
    teta=acos((dx1*dx2+dy1*dy2)/(d1*d2))
    teta=(teta*180)/pi
    If(teta>=Teta0) Then
     a=a+1
     bs(a)=bp1
     
     counter1=counter1+1
     NOCP(counter1)=bp1
    End If 
   End If 
   
   !part 10.2:
   If (I/=J .AND. (bp2==bp3 .OR. bp2==bp4)) Then
    
    dx1=x(bp2)-x(bp1)
    dy1=y(bp2)-y(bp1)
    dx2=x(bp4)-x(bp3)
    dy2=y(bp4)-y(bp3)
    d1=sqrt(dx1*dx1+dy1*dy1)
    d2=sqrt(dx2*dx2+dy2*dy2)
    teta=acos((dx1*dx2+dy1*dy2)/(d1*d2))
    teta=(teta*180)/pi
    If(teta>=Teta0) Then
     a=a+1
     bs(a)=bp2
     
     counter1=counter1+1
     NOCP(counter1)=bp2
     
    End If 
   End If 
  End Do 
 End Do 
 
 !part 11:
 MIS=0
 COARSE=0
 Do I=NF2+1,NF
  counter=0
  p1=IDS(3,I)
  p2=IDS(4,I)
  
  !part 11.1:
  Do K=1,a
   If (p1==bs(K) .OR. p2==bs(K)) Then
    counter=counter+1
    Exit
   End If    
  End Do   
  
   If (counter/=0) Then
    Cycle
   End If  
   
  !part 11.2: 
  If(d(p1,p2)<((F1(p1)+F1(p2))/Beta)) Then
   MIS=MIS+1
   point1(MIS)=p1
   point2(MIS)=p2
   
  !part 11.3: 
  Else
   COARSE=COARSE+1
   cp1(COARSE)=p1
   cp2(COARSE)=p2
  End If   
 End Do  
 
 !part 12:
 n=0
 Do I=1,COARSE
  counter=0
  Do J=1,MIS
   If (cp1(I)==point1(J) .OR. cp1(I)==point2(J)) Then
    counter=counter+1
    Exit
   End If
  End Do
  
  If (counter==0) Then
   n=n+1
   coarsepoint(n)=cp1(i)
  End If
 End Do
 
 Do I=1,COARSE
  counter=0
  Do J=1,MIS
   If (cp2(I)==point1(J) .OR. cp2(I)==point2(J)) Then
    counter=counter+1
    Exit
   End If
  End Do
  
  If (counter==0) Then
   n=n+1
   coarsepoint(n)=cp2(i)
  End If
 End Do
 
 !part 13:
 Do I=1,n
  counter=0
  Do J=1,I
   If (I/=J .AND. coarsepoint(I)==coarsepoint(J)) Then
    counter=counter+1
    Exit
   End If 
  End Do
  
   If (counter==0) Then
   P=coarsepoint(I)
   
     counter1=counter1+1
     NOCP(counter1)=P
   
   End If
 End Do  
 
 !part 14:
 m=0
 Do I=1,MIS
  counter=0
  q=24000
  P=point1(I)
  
10   Do K=1,m
  
   !part 14.1:
   If (I==DelRow(K)) Then
    q=K
    Exit
   End If
   
   !part 14.2:
   If (P==NICM(K)) Then
    counter=counter+1
    P=point2(I)
    If (counter==1) Then
     goto 10
    End If
   End If
   
   !part 14.3
   If (counter==2) Then
    Exit
   End If
   
  End Do
  
  !part 14.4:
  If (counter==2 .OR. I==DelRow(q)) Then
   Cycle
  End If
  
  !part 15:
  Do J=1,MIS
   
   !part 15.1:
   If (P==point1(J)) Then
    m=m+1
    DelRow(m)=J
    NICM(m)=point2(J)
   End If
   
   !part 15.2:
   If (P==point2(J)) Then
    m=m+1
    DelRow(m)=J
    NICM(m)=point1(J)
   End If
  End Do
  
     counter1=counter1+1
     NOCP(counter1)=P
  
 End do
 

 !************************************************************
 
 !part 16:
 MIS=0
 COARSE=0
 Do I=NF1+1,NF2
  counter=0
  p1=IDS(3,I)
  p2=IDS(4,I)
  
    Do J=NF2+1,NF
     If (p1==IDS(3,J) .OR. p1==IDS(4,J) .OR. p2==IDS(3,J) .OR. p2==IDS(4,J)) Then
      counter=counter+1
      exit
     End If
    End Do
  
   If (counter/=0) Then
    Cycle
   End If  
   
  If(d(p1,p2)<((f1(p1)+f1(p2))/Beta)) Then
   MIS=MIS+1
   point1(MIS)=p1
   point2(MIS)=p2
   
  Else
   COARSE=COARSE+1
   cp1(COARSE)=p1
   cp2(COARSE)=p2
  End If   
 End Do  
 
 !part 17:
 n=0
 Do I=1,COARSE
  counter=0
  Do J=1,MIS
   If (cp1(I)==point1(J) .OR. cp1(I)==point2(J)) Then
    counter=counter+1
    Exit
   End If
  End Do
  
  If (counter==0) Then
   n=n+1
   coarsepoint(n)=cp1(i)
  End If
 End Do
 
 Do I=1,COARSE
  counter=0
  Do J=1,MIS
   If (cp2(I)==point1(J) .OR. cp2(I)==point2(J)) Then
    counter=counter+1
    Exit
   End If
  End Do
  
  If (counter==0) Then
   n=n+1
   coarsepoint(n)=cp2(i)
  End If
 End Do
 
 !part 18:
 Do I=1,n
  counter=0
  Do J=1,I
   If (I/=J .AND. coarsepoint(I)==coarsepoint(J)) Then
    counter=counter+1
    Exit
   End If 
  End Do
  
   If (counter==0) Then
   P=coarsepoint(I)
   
     counter1=counter1+1
     NOCP(counter1)=P
   
   End If
 End Do  
 
 !part 19:
 m=0
 Do I=1,MIS
  counter=0
  q=24000
  P=point1(I)
  
20  Do K=1,m
   
   If (I==DelRow(K)) Then
    q=K
    Exit
   End If
   
   If (P==NICM(K)) Then
    counter=counter+1
    P=point2(I)
    If (counter==1) Then
     goto 20
    End If
   End If
   
   If (counter==2) Then
    Exit
   End If
   
  End Do
  
  If (counter==2 .OR. I==DelRow(q)) Then
   Cycle
  End If
  
  !part 20:
  Do J=1,MIS
   
   If (P==point1(J)) Then
    m=m+1
    DelRow(m)=J
    NICM(m)=point2(J)
   End If
   
   If (P==point2(J)) Then
    m=m+1
    DelRow(m)=J
    NICM(m)=point1(J)
   End If
  End Do
  
     counter1=counter1+1
     NOCP(counter1)=P
  
 End do
 
 !part 21:
 Do I=1,Np
    counter2=0
    Do J=1,Counter1
        IF(I==NOCP(J))Then
            counter2=1
            exit
        End IF    
    End Do
    
    if(counter2==0)then
       write(2,*) I
    end if
 End Do
 

 
 Close(1)
 Close(2)
 call cpu_time(finish)
 print *, finish-start
 pause
 end program Coarsening
!*******************************************************************************************
 !*******************************************************************************************
 Subroutine Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim
 Intent(Out  )::NP,NC,NF,NR,NFR,BC,IDS,X,Y

 Integer::Dim,I,J,J1,JJ,NP,NC,NF,NR,SFace,FaceType,MeshDim
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
!*******************************************************************************************
!Part 1:
 Open(1,File='Mesh.Txt')

!Part 2:
 Read(1,*) MeshDim
 IF(MeshDim/=2)Print*,'Please Check the Mesh File. It is not a 2D Mesh'

!Part 3:
 Read(1,*) NP 

!Part 4:
 Read(1,*) NC

!Part 5:
 Read(1,*) NF

!Part 6:
 Read(1,*) NR

!Part 7:
 Read(1,*)
 Do J=1,NR
 Read(1,*) NFR(J) , BC(J)
 End Do 

!Part 8:
 Read(1,*)
 Do J=1,NF  
 Read(1,*) FaceType,IDS(1,J),IDS(2,J),IDS(3,J),IDS(4,J)
 End Do

!Part 9:
 Read(1,*)  
 Do J=1,NP
 Read(1,*) X(J),Y(J)
 End Do

 Close(1)
!*******************************************************************************************
 End
!###########################################################################################
!###########################################################################################
  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Renumbering Faces According to Boundary Conditions and      //!
!//     Determining the Index of first and last Faces of Each Boundary Conditions/!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F022F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NR,NF
 Intent(Out  )::NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFIF1,NFIF2
 Intent(InOut)::NFR,BC,IDS

 Integer::Dim,J,JJ,J1,I,SF,N,M,NR,NF,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,&
          NF1,NF2,NFIF1,NFIF2,NFN,NFW,NFF,NFI,NFO,NFS,NFIF
 Integer,Dimension(1:100)::NFR,TNFR,BC,TBC
 Integer,Dimension(1:4,1:Dim)::IDS,TIDS
!*******************************************************************************************
!Part 1:
 Do J=1,NF
	Do J1=1,4
       TIDS(J1,J) = IDS(J1,J)
    End do
 End Do

!Part 2:
 Do J=1,NR
    TNFR(J) = NFR(J)
	TBC(J)  = BC(J) 
 End do

!Part 3:
 N=0
 M=0
 Do JJ=1,10

   !Part 4:
    SF=0
    Do J=1,NR
       IF(TBC(J)==JJ)Then

	    Do I=SF+1,SF+TNFR(J)

           N=N+1
           Do J1=1,4
              IDS(J1,N) = TIDS(J1,I)
		   End do

        Enddo
		
	   !Part 5:		   
	    M=M+1
	    NFR(M) = TNFR(J)
		BC(M)  = TBC(J)  

	   Endif
	   SF=SF+TNFR(J)
    End Do

 End Do

!Part 6:
 NFN  = 0
 NFW  = 0
 NFF  = 0
 NFI  = 0
 NFO  = 0
 NFS  = 0
 NFIF = 0
 Do J=1,NR
    IF( BC(J)==1 ) NFN  = NFN  + NFR(J)
    IF( BC(J)==2 ) NFW  = NFW  + NFR(J)
    IF( BC(J)==3 ) NFF  = NFF  + NFR(J)
    IF( BC(J)==4 ) NFI  = NFI  + NFR(J)
    IF( BC(J)==5 ) NFO  = NFO  + NFR(J)
    IF( BC(J)==6 ) NFS  = NFS  + NFR(J)
    IF( BC(J)==7 ) NFIF = NFIF + NFR(J)
 End Do

!Part 7:
 NF1=0
 NF2=NF1+NFN

 NFW1=NF2
 NFW2=NFW1+NFW

 NFF1=NFW2
 NFF2=NFF1+NFF

 NFI1=NFF2
 NFI2=NFI1+NFI

 NFO1=NFI2
 NFO2=NFO1+NFO

 NFS1=NFO2
 NFS2=NFS1+NFS

 NFIF1=NFI2
 NFIF2=NFIF1+NFIF
!*******************************************************************************************
 End