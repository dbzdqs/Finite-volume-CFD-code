!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Interpolation between two Meshes by intersection of triangles           //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: January, 07, 2017                                                              //!
!// Developed by: M. Hashemi, Iran, Tehran, Mohammadhashemi@ut.ac.ir                     //!
!// Doc ID: MC2F100F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
Subroutine Interpolation(Dim,WNP1,RES,Error,NOL,FTC)
implicit none
!***********************************************************************************************************
Intent (IN)::Dim,NOL
Intent (InOut)::WNP1,RES,Error
!===============================
 Integer::Dim,I,J,K,NC,NP,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
          NR,NS,counter,NMF,NP1,NC1,NF_1,NF_2,NR1,NP2,NC2,NR2,P1,P2,ME,NE,h1,h2,counter10,counter5,counter6,num1,counter3,counter9,counter11,counter13,NP2p,NP1p,NOL,FTC
 Real(8)::F1,F2,F3,F_p1,F_p2,F_p3,F_s1,F_s2,F_s3,F_h,F_h1,F_h2,F_h3,collisionX,collisionY,M1,M2,Multi1,Multi2
 Integer,Dimension(1:100)::NFR,BC,NFR1,BC1,NFR2,BC2
 Integer,Dimension(1:3)::P_2,P_1
 Integer,Dimension(1:4,1:Dim)::IDS1,IDS2,IDS,corn
 Integer,Dimension(1:3,1:Dim)::Corn1,Corn2
 Real(8),Dimension(1:4,1:Dim)::WNP1,WC,Con,RES,Error,delta_s0,delta_s1,delta_s2,delta_sm,delta_sn,delta_sp,FUNC1,FUNC2,FUNC3,PFUNC1,PFUNC2,PFUNC3
 Real(8),Dimension(1:Dim)::X,Y,XC,YC,A,NX,NY,DA,DT,P,X2,Y2,A1,X1,Y1
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:3)::M_1,M_2
!************************************************************************************************************
!Part 1:
Call Read_2DMesh1(Dim,NP1,NC1,NF_1,NR1,NFR1,BC1,IDS1,X1,Y1,NOL,FTC)
Call Read_2DMesh2(Dim,NP2,NC2,NF_2,NR2,NFR2,BC2,IDS2,X2,Y2,NOL,FTC)

!Part 2:
NP1p=NP1
NP2p=NP2

!Part 3:
Do I=1,NC1
    FUNC1(1,I)=WNP1(1,I);FUNC1(2,I)=WNP1(2,I);FUNC1(3,I)=WNP1(3,I);FUNC1(4,I)=WNP1(4,I)
    FUNC2(1,I)=RES(1,I);FUNC2(2,I)=RES(2,I);FUNC2(3,I)=RES(3,I);FUNC2(4,I)=RES(4,I)
    FUNC3(1,I)=Error(1,I);FUNC3(2,I)=Error(2,I);FUNC3(3,I)=Error(3,I);FUNC3(4,I)=Error(4,I)
end Do


!Part 4:
Do I=1,NC2
    delta_s0(1,I)=0;delta_s0(2,I)=0;delta_s0(3,I)=0;delta_s0(4,I)=0
    delta_s1(1,I)=0;delta_s1(2,I)=0;delta_s1(3,I)=0;delta_s1(4,I)=0
    delta_s2(1,I)=0;delta_s2(2,I)=0;delta_s2(3,I)=0;delta_s2(4,I)=0
    delta_sm(1,I)=0;delta_sm(2,I)=0;delta_sm(3,I)=0;delta_sm(4,I)=0
    delta_sn(1,I)=0;delta_sn(2,I)=0;delta_sn(3,I)=0;delta_sn(4,I)=0
    delta_sp(1,I)=0;delta_sp(2,I)=0;delta_sp(3,I)=0;delta_sp(4,I)=0
end Do

!Part 5:
if(NOL==1 .AND. FTC==1)then
NMF=2
else if(NOL==1 .AND. FTC==0)then
NMF=1
else if(NOL==2 .AND. FTC==1)then
NMF=3
else if(NOL==2 .AND. FTC==0)then
NMF=2
end if

!Part 6:
Call Read_2DMeshMG(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)

!Part 7:
Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

!Part 8:
Call Edge_To_Cell(Dim,NF,NC,IDS,Corn)

!Part 9:
Do J=1,NC2
    Corn2(1,J)=Corn(1,J)
    Corn2(2,J)=Corn(2,J)
    Corn2(3,J)=Corn(3,J)
End Do


!Part 10:
Do I=1,NF_1
    
    !Part 11:
    counter5=0
    counter6=0
    counter13=0
    
    !Part 12:
    ME=IDS1(1,I)
    NE=IDS1(2,I)
    P1=IDS1(3,I)
    P2=IDS1(4,I)
   
    !Part 13:
10  Do J=1,NC2
        
        !Part 14:
        P_2(1)=Corn2(1,J)
        P_2(2)=Corn2(2,J)
        P_2(3)=Corn2(3,J)
        
    !Part 15:
    F1 = (X2(P_2(1))-X1(P1))*(Y2(P_2(2))-Y1(P1))-(X2(P_2(2))-X1(P1))*(Y2(P_2(1))-Y1(P1))
    F2 = (X2(P_2(2))-X1(P1))*(Y2(P_2(3))-Y1(P1))-(X2(P_2(3))-X1(P1))*(Y2(P_2(2))-Y1(P1))
    F3 = (X2(P_2(3))-X1(P1))*(Y2(P_2(1))-Y1(P1))-(X2(P_2(1))-X1(P1))*(Y2(P_2(3))-Y1(P1))  
    
    !Part 16:
    IF(abs(F1)<=1e-10)Then
        F1=0
    End IF
    IF(abs(F2)<=1e-10)Then
        F2=0
    End IF
    IF(abs(F3)<=1e-10)Then
        F3=0
    End IF 
    
    !Part 17:
    IF((F1>=0. .And. F2>=0. .And. F3>=0.))Then    
        
        counter5=counter5+1
        
        !Part 18:
        F_p1 = (X2(P_2(1))-X1(P2))*(Y2(P_2(2))-Y1(P2))-(X2(P_2(2))-X1(P2))*(Y2(P_2(1))-Y1(P2))
        F_p2 = (X2(P_2(2))-X1(P2))*(Y2(P_2(3))-Y1(P2))-(X2(P_2(3))-X1(P2))*(Y2(P_2(2))-Y1(P2))
        F_p3 = (X2(P_2(3))-X1(P2))*(Y2(P_2(1))-Y1(P2))-(X2(P_2(1))-X1(P2))*(Y2(P_2(3))-Y1(P2))
    
    !Part 19:
    IF(abs(F_p1)<=1e-16)Then
          F_p1=0
    End IF      
    IF(abs(F_p2)<=1e-16)Then
          F_p2=0
    End IF
    IF(abs(F_p3)<=1e-16)Then
          F_p3=0
    End IF
        
        !Part 20:
        IF((abs(F1)<=1e-10 .AND. abs(F_p1)<=1e-10 .AND. F_p2>=0 .AND. F_p3>=0) .OR. (abs(F2)<=1e-10 .AND. abs(F_p2)<=1e-10 .AND. F_p1>=0 .AND. F_p3>=0) .OR. (abs(F3)<=1e-10 .AND. abs(F_p3)<=1e-10 .AND. F_p1>=0 .AND. F_p2>=0))Then
            
            counter13=counter13+1
            Exit
            
        !Part 21:            
        else IF((F_p1>=0. .And. F_p2>=0. .And. F_p3>=0.))Then
            
            counter13=counter13+1
            
            If(NE==0)then
                
            delta_s0(1,J)= delta_s0(1,J) + 0.5*(FUNC1(1,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s0(2,J)= delta_s0(2,J) + 0.5*(FUNC1(2,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s0(3,J)= delta_s0(3,J) + 0.5*(FUNC1(3,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s0(4,J)= delta_s0(4,J) + 0.5*(FUNC1(4,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
                        
            delta_s1(1,J)= delta_s1(1,J) + 0.5*(FUNC2(1,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s1(2,J)= delta_s1(2,J) + 0.5*(FUNC2(2,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s1(3,J)= delta_s1(3,J) + 0.5*(FUNC2(3,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s1(4,J)= delta_s1(4,J) + 0.5*(FUNC2(4,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))

            delta_s2(1,J)= delta_s2(1,J) + 0.5*(FUNC3(1,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s2(2,J)= delta_s2(2,J) + 0.5*(FUNC3(2,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s2(3,J)= delta_s2(3,J) + 0.5*(FUNC3(3,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s2(4,J)= delta_s2(4,J) + 0.5*(FUNC3(4,ME))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))

            else
                
            delta_s0(1,J)= delta_s0(1,J) + 0.5*(FUNC1(1,ME)-FUNC1(1,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s0(2,J)= delta_s0(2,J) + 0.5*(FUNC1(2,ME)-FUNC1(2,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s0(3,J)= delta_s0(3,J) + 0.5*(FUNC1(3,ME)-FUNC1(3,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s0(4,J)= delta_s0(4,J) + 0.5*(FUNC1(4,ME)-FUNC1(4,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
                       
            delta_s1(1,J)= delta_s1(1,J) + 0.5*(FUNC2(1,ME)-FUNC2(1,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s1(2,J)= delta_s1(2,J) + 0.5*(FUNC2(2,ME)-FUNC2(2,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s1(3,J)= delta_s1(3,J) + 0.5*(FUNC2(3,ME)-FUNC2(3,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s1(4,J)= delta_s1(4,J) + 0.5*(FUNC2(4,ME)-FUNC2(4,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            
            delta_s2(1,J)= delta_s2(1,J) + 0.5*(FUNC3(1,ME)-FUNC3(1,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s2(2,J)= delta_s2(2,J) + 0.5*(FUNC3(2,ME)-FUNC3(2,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s2(3,J)= delta_s2(3,J) + 0.5*(FUNC3(3,ME)-FUNC3(3,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))
            delta_s2(4,J)= delta_s2(4,J) + 0.5*(FUNC3(4,ME)-FUNC3(4,NE))*(X1(P1)*Y1(P2)-X1(P2)*Y1(P1))


            end if
                        
            Exit
            
        !Part 22:   
        else    
            
            !Part 23:
            counter10=0
            if(X1(P2)==X1(P1))Then
                
                collisionX=X1(P2)
                counter10=counter10+1
                
            else
                
                M1=(Y1(P2)-Y1(P1))/(X1(P2)-X1(P1))
                
            End IF
            
            !Part 24:
            Do K=1,3
                
                !Part 25:
                h1=K;h2=K+1
                if(K==3)Then
                    h2=1
                    h1=3  
                End if
                
                !Part 26:
                counter11=0
                if(X2(P_2(h2))==X2(P_2(h1)))Then
                
                    collisionX=X2(P_2(h2))
                    counter11=counter11+1
                    
                else
            
                    M_2(K)=(Y2(P_2(h2))-Y2(P_2(h1)))/(X2(P_2(h2))-X2(P_2(h1)))
                    
                End if
                
                !Part 27:
                
                !Part 27.1:
                if(counter10==1 .AND. counter11==0)then
                    
                    collisionY=Y2(P_2(h1))+M_2(K)*(collisionX-X2(P_2(h1)))
                    
                else if(counter10==0 .AND. counter11==1)then
                    
                    collisionY=Y1(P1)+M1*(collisionX-X1(P1))
                    
                end if
                
                !Part 27.2:
                If(counter10==1 .AND. counter11==1)Then
                
                    cycle
                    
                else If(M1==M_2(K))Then
                
                    cycle
            
                !Part 27.3:
                else if(counter10==0 .AND. counter11==0)Then
                    
                    collisionX=(M1*X1(P1)-M_2(K)*X2(P_2(K))+Y2(P_2(K))-Y1(P1))/(M1-M_2(K))
                    collisionY=M1*(collisionX-X1(P1))+Y1(P1)
                    
                End If
                
                
           !Part 28:
           !Part 28.1:
           F_s1 = (X2(P_2(1))-collisionX)*(Y2(P_2(2))-collisionY)-(X2(P_2(2))-collisionX)*(Y2(P_2(1))-collisionY)
           F_s2 = (X2(P_2(2))-collisionX)*(Y2(P_2(3))-collisionY)-(X2(P_2(3))-collisionX)*(Y2(P_2(2))-collisionY)
           F_s3 = (X2(P_2(3))-collisionX)*(Y2(P_2(1))-collisionY)-(X2(P_2(1))-collisionX)*(Y2(P_2(3))-collisionY)

           !Part 28.2:
           IF(abs(F_s1)<=1e-10)Then
               F_s1=0
           end if
           IF(abs(F_s2)<=1e-10)Then
               F_s2=0
           end if
           IF(abs(F_s3)<=1e-10)Then
               F_s3=0
           End IF
           
 
           !Part 28.3:
           Multi1=(Y1(P2)-collisionY)*(collisionY-Y1(P1))          
           Multi2=(X1(P2)-collisionX)*(collisionX-X1(P1))  
           
           !Part 28.4:
           IF(abs(Multi1)<=1e-10 .AND. abs(Multi2)<=1e-10)Then
                Multi1=0
                Multi2=0
           else IF(abs(Multi1)<=1e-10)Then
                Multi1=0
           else IF(abs(Multi2)<=1e-10)Then
                Multi2=0
           End IF
           
           !Part 29:
           If((F_s1>=0. .And. F_s2>=0. .And. F_s3>=0. .And. Multi1>0) .OR. (F_s1>=0. .And. F_s2>=0. .And. F_s3>=0. .And. Multi2>0) .OR. (F_s1<=0. .And. F_s2<=0. .And. F_s3<=0. .And. Multi1>0) .OR. (F_s1<=0. .And. F_s2<=0. .And. F_s3<=0. .And. Multi2>0))Then
            
            !Part 29.1:
            IF((abs(F1)<=1e-10 .AND. abs(F_s1)<=1e-10) .OR. (abs(F2)<=1e-10 .AND. abs(F_s2)<=1e-10) .OR. (abs(F3)<=1e-10 .AND. abs(F_s3)<=1e-10))Then
                   
                   NP1p=NP1p+1
                   P1=NP1p
                   X1(P1)=collisionX
                   Y1(P1)=collisionY
            
                   goto 10
            !Part 29.2:                   
            else
                
                IF(counter6==1)Then
                    
                     If(NE==0)then
                    
                        delta_s0(1,J)= delta_s0(1,J) + 0.5*(FUNC1(1,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s0(2,J)= delta_s0(2,J) + 0.5*(FUNC1(2,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s0(3,J)= delta_s0(3,J) + 0.5*(FUNC1(3,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s0(4,J)= delta_s0(4,J) + 0.5*(FUNC1(4,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        
                        delta_s1(1,J)= delta_s1(1,J) + 0.5*(FUNC2(1,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s1(2,J)= delta_s1(2,J) + 0.5*(FUNC2(2,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s1(3,J)= delta_s1(3,J) + 0.5*(FUNC2(3,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s1(4,J)= delta_s1(4,J) + 0.5*(FUNC2(4,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))

                        delta_s2(1,J)= delta_s2(1,J) + 0.5*(FUNC3(1,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s2(2,J)= delta_s2(2,J) + 0.5*(FUNC3(2,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s2(3,J)= delta_s2(3,J) + 0.5*(FUNC3(3,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s2(4,J)= delta_s2(4,J) + 0.5*(FUNC3(4,ME))*(-X1(P1)*collisionY+collisionX*Y1(P1))

                    else
                                   
                        delta_s0(1,J)= delta_s0(1,J) + 0.5*(FUNC1(1,ME)-FUNC1(1,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s0(2,J)= delta_s0(2,J) + 0.5*(FUNC1(2,ME)-FUNC1(2,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s0(3,J)= delta_s0(3,J) + 0.5*(FUNC1(3,ME)-FUNC1(3,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s0(4,J)= delta_s0(4,J) + 0.5*(FUNC1(4,ME)-FUNC1(4,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        
                        delta_s1(1,J)= delta_s1(1,J) + 0.5*(FUNC2(1,ME)-FUNC2(1,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s1(2,J)= delta_s1(2,J) + 0.5*(FUNC2(2,ME)-FUNC2(2,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s1(3,J)= delta_s1(3,J) + 0.5*(FUNC2(3,ME)-FUNC2(3,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s1(4,J)= delta_s1(4,J) + 0.5*(FUNC2(4,ME)-FUNC2(4,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))

                        delta_s2(1,J)= delta_s2(1,J) + 0.5*(FUNC3(1,ME)-FUNC3(1,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s2(2,J)= delta_s2(2,J) + 0.5*(FUNC3(2,ME)-FUNC3(2,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s2(3,J)= delta_s2(3,J) + 0.5*(FUNC3(3,ME)-FUNC3(3,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))
                        delta_s2(4,J)= delta_s2(4,J) + 0.5*(FUNC3(4,ME)-FUNC3(4,NE))*(-X1(P1)*collisionY+collisionX*Y1(P1))

                        
                    end if
                    
                else
                    
                     If(NE==0)then
                    
                        delta_s0(1,J)= delta_s0(1,J) + 0.5*(FUNC1(1,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s0(2,J)= delta_s0(2,J) + 0.5*(FUNC1(2,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s0(3,J)= delta_s0(3,J) + 0.5*(FUNC1(3,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s0(4,J)= delta_s0(4,J) + 0.5*(FUNC1(4,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                                       
                        delta_s1(1,J)= delta_s1(1,J) + 0.5*(FUNC2(1,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s1(2,J)= delta_s1(2,J) + 0.5*(FUNC2(2,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s1(3,J)= delta_s1(3,J) + 0.5*(FUNC2(3,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s1(4,J)= delta_s1(4,J) + 0.5*(FUNC2(4,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        
                        delta_s2(1,J)= delta_s2(1,J) + 0.5*(FUNC3(1,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s2(2,J)= delta_s2(2,J) + 0.5*(FUNC3(2,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s2(3,J)= delta_s2(3,J) + 0.5*(FUNC3(3,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s2(4,J)= delta_s2(4,J) + 0.5*(FUNC3(4,ME))*(X1(P1)*collisionY-collisionX*Y1(P1))


                    else
                                   
                        delta_s0(1,J)= delta_s0(1,J) + 0.5*(FUNC1(1,ME)-FUNC1(1,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s0(2,J)= delta_s0(2,J) + 0.5*(FUNC1(2,ME)-FUNC1(2,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s0(3,J)= delta_s0(3,J) + 0.5*(FUNC1(3,ME)-FUNC1(3,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s0(4,J)= delta_s0(4,J) + 0.5*(FUNC1(4,ME)-FUNC1(4,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                    
                        delta_s1(1,J)= delta_s1(1,J) + 0.5*(FUNC2(1,ME)-FUNC2(1,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s1(2,J)= delta_s1(2,J) + 0.5*(FUNC2(2,ME)-FUNC2(2,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s1(3,J)= delta_s1(3,J) + 0.5*(FUNC2(3,ME)-FUNC2(3,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s1(4,J)= delta_s1(4,J) + 0.5*(FUNC2(4,ME)-FUNC2(4,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))

                        delta_s2(1,J)= delta_s2(1,J) + 0.5*(FUNC3(1,ME)-FUNC3(1,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s2(2,J)= delta_s2(2,J) + 0.5*(FUNC3(2,ME)-FUNC3(2,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s2(3,J)= delta_s2(3,J) + 0.5*(FUNC3(3,ME)-FUNC3(3,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))
                        delta_s2(4,J)= delta_s2(4,J) + 0.5*(FUNC3(4,ME)-FUNC3(4,NE))*(X1(P1)*collisionY-collisionX*Y1(P1))

                    end if
                    
                End IF
                
            
            NP1p=NP1p+1
            P1=NP1p
            X1(P1)=collisionX
            Y1(P1)=collisionY
            
            goto 10
            
            End If
            
           End IF
           
            
        End Do ! Do K
            
            
        End IF ! for intersection between edge and triangle of mesh 2
        
    !Part 30:    
    else
        
      Cycle  
      
    End IF
          
    End Do ! Do J

       !Part 31:
       if((counter5==0 .AND. counter6==0) .OR. (counter13==0 .AND. counter5/=0))Then
          
          num1=P1
          P1=P2
          P2=num1
          counter6=1
          
          !Part 31.1:
          if(counter5/=0)then
              counter13=1
          end if
          
          goto 10
          
      end if
            
End Do ! Do I

!********************************************************************************************************************************************
!Part 32:
if(NOL==1 .AND. FTC==1)then
NMF=1
else if(NOL==1 .AND. FTC==0)then
NMF=2
else if(NOL==2 .AND. FTC==1)then
NMF=2
else if(NOL==2 .AND. FTC==0)then
NMF=3
end if

!Part 33:
Call Read_2DMeshMG(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF) 
!Part 34:
Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
!Part 35:
Call Edge_To_Cell(Dim,NF,NC,IDS,Corn)

!Part 36:
Do J=1,NC1
    Corn1(1,J)=Corn(1,J)
    Corn1(2,J)=Corn(2,J)
    Corn1(3,J)=Corn(3,J)
End Do

!Part 37:
Do I=1,NF_2
    
    !Part 38:
    counter3=0
    counter5=0
    counter6=0
    counter13=0
    
    !Part 39:
    ME=IDS2(1,I)
    NE=IDS2(2,I)
    P1=IDS2(3,I)
    P2=IDS2(4,I)
     
    !Part 40:
20  Do J=1,NC1
    
        !Part 41:
        P_1(1)=Corn1(1,J)
        P_1(2)=Corn1(2,J)
        P_1(3)=Corn1(3,J)
     
    !Part 42:   
    F1 = (X1(P_1(1))-X2(P1))*(Y1(P_1(2))-Y2(P1))-(X1(P_1(2))-X2(P1))*(Y1(P_1(1))-Y2(P1))
    F2 = (X1(P_1(2))-X2(P1))*(Y1(P_1(3))-Y2(P1))-(X1(P_1(3))-X2(P1))*(Y1(P_1(2))-Y2(P1))
    F3 = (X1(P_1(3))-X2(P1))*(Y1(P_1(1))-Y2(P1))-(X1(P_1(1))-X2(P1))*(Y1(P_1(3))-Y2(P1))
    
    !Part 43:    
    IF(abs(F1)<=1e-10)Then
          F1=0
    end if
    IF(abs(F2)<=1e-10)Then
          F2=0
    end if
    IF(abs(F3)<=1e-10)Then
          F3=0
     End IF 
     
    !Part 44:
    IF((F1>=0. .And. F2>=0. .And. F3>=0.))Then    
        
        counter5=counter5+1
    
    !Part 45:    
    F_p1 = (X1(P_1(1))-X2(P2))*(Y1(P_1(2))-Y2(P2))-(X1(P_1(2))-X2(P2))*(Y1(P_1(1))-Y2(P2))
    F_p2 = (X1(P_1(2))-X2(P2))*(Y1(P_1(3))-Y2(P2))-(X1(P_1(3))-X2(P2))*(Y1(P_1(2))-Y2(P2))
    F_p3 = (X1(P_1(3))-X2(P2))*(Y1(P_1(1))-Y2(P2))-(X1(P_1(1))-X2(P2))*(Y1(P_1(3))-Y2(P2))
    
     !Part 46:
     IF(abs(F_p1)<=1e-16)Then
          F_p1=0
     else IF(abs(F_p2)<=1e-16)Then
          F_p2=0
     else IF(abs(F_p3)<=1e-16)Then
          F_p3=0
     End IF

        !Part 47:    
        IF((abs(F1)<=1e-10 .AND. abs(F_p1)<=1e-10 .AND. F_p2>=0 .AND. F_p3>=0) .OR. (abs(F2)<=1e-10 .AND. abs(F_p2)<=1e-10 .AND. F_p1>=0 .AND. F_p3>=0) .OR. (abs(F3)<=1e-10 .AND. abs(F_p3)<=1e-10 .AND. F_p1>=0 .AND. F_p2>=0))Then
            
            counter3=counter3+1
            counter13=counter13+1
          
            !Part 47.1:
            F_h1=(X2(P1)-X1(P_1(1)))*(Y2(P2)-Y1(P_1(1)))-(X2(P2)-X1(P_1(1)))*(Y2(P1)-Y1(P_1(1)))
            F_h2=(X2(P1)-X1(P_1(2)))*(Y2(P2)-Y1(P_1(2)))-(X2(P2)-X1(P_1(2)))*(Y2(P1)-Y1(P_1(2)))
            F_h3=(X2(P1)-X1(P_1(3)))*(Y2(P2)-Y1(P_1(3)))-(X2(P2)-X1(P_1(3)))*(Y2(P1)-Y1(P_1(3)))
            
     !Part 47.2:      
     IF(abs(F_h1)<=1e-16)Then
          F_h1=0
     End if
     IF(abs(F_h2)<=1e-16)Then
          F_h2=0
     End if
     IF(abs(F_h3)<=1e-16)Then
          F_h3=0
     End IF
            
            !Part 47.3:
            if((X1(P_1(1))/=X2(P1) .AND. X1(P_1(1))/=X2(P2) .AND. F_h1/=0) .OR. (Y1(P_1(1))/=Y2(P1) .AND. Y1(P_1(1))/=Y2(P2)) .AND. F_h1/=0)then
                
                F_h=F_h1
                
            else if((X1(P_1(2))/=X2(P1) .AND. X1(P_1(2))/=X2(P2) .AND. F_h2/=0) .OR. (Y1(P_1(2))/=Y2(P1) .AND. Y1(P_1(2))/=Y2(P2)) .AND. F_h2/=0)then
                
                F_h=F_h2
                
            else if((X1(P_1(3))/=X2(P1) .AND. X1(P_1(3))/=X2(P2) .AND. F_h3/=0) .OR. (Y1(P_1(3))/=Y2(P1) .AND. Y1(P_1(3))/=Y2(P2)) .AND. F_h3/=0)then
                
                F_h=F_h3
                
            end if
            
            
            if(F_h>0)then
                
                delta_sm(1,ME)= delta_sm(1,ME) + 0.5*(FUNC1(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sm(2,ME)= delta_sm(2,ME) + 0.5*(FUNC1(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sm(3,ME)= delta_sm(3,ME) + 0.5*(FUNC1(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sm(4,ME)= delta_sm(4,ME) + 0.5*(FUNC1(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))

                delta_sn(1,ME)= delta_sn(1,ME) + 0.5*(FUNC2(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sn(2,ME)= delta_sn(2,ME) + 0.5*(FUNC2(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sn(3,ME)= delta_sn(3,ME) + 0.5*(FUNC2(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sn(4,ME)= delta_sn(4,ME) + 0.5*(FUNC2(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))

                delta_sp(1,ME)= delta_sp(1,ME) + 0.5*(FUNC3(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sp(2,ME)= delta_sp(2,ME) + 0.5*(FUNC3(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sp(3,ME)= delta_sp(3,ME) + 0.5*(FUNC3(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sp(4,ME)= delta_sp(4,ME) + 0.5*(FUNC3(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))


            else if( F_h<0 .AND. NE/=0)then
                
                delta_sm(1,NE)= delta_sm(1,NE) - 0.5*(FUNC1(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sm(2,NE)= delta_sm(2,NE) - 0.5*(FUNC1(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sm(3,NE)= delta_sm(3,NE) - 0.5*(FUNC1(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sm(4,NE)= delta_sm(4,NE) - 0.5*(FUNC1(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                
                delta_sn(1,NE)= delta_sn(1,NE) - 0.5*(FUNC2(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sn(2,NE)= delta_sn(2,NE) - 0.5*(FUNC2(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sn(3,NE)= delta_sn(3,NE) - 0.5*(FUNC2(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sn(4,NE)= delta_sn(4,NE) - 0.5*(FUNC2(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                
                delta_sp(1,NE)= delta_sp(1,NE) - 0.5*(FUNC3(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sp(2,NE)= delta_sp(2,NE) - 0.5*(FUNC3(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sp(3,NE)= delta_sp(3,NE) - 0.5*(FUNC3(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
                delta_sp(4,NE)= delta_sp(4,NE) - 0.5*(FUNC3(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))

            end if
            
        !Part 47.4:    
        if(counter3==1 .AND. NE/=0)then
            cycle
        else
            exit
        end if        
         
        !Part 48:
        else IF((F_p1>=0. .And. F_p2>=0. .And. F_p3>=0.))Then
            
            counter13=counter13+1
            
            If(NE==0)then
                
            delta_sm(1,ME)= delta_sm(1,ME) + 0.5*(FUNC1(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(2,ME)= delta_sm(2,ME) + 0.5*(FUNC1(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(3,ME)= delta_sm(3,ME) + 0.5*(FUNC1(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(4,ME)= delta_sm(4,ME) + 0.5*(FUNC1(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))

            delta_sn(1,ME)= delta_sn(1,ME) + 0.5*(FUNC2(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(2,ME)= delta_sn(2,ME) + 0.5*(FUNC2(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(3,ME)= delta_sn(3,ME) + 0.5*(FUNC2(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(4,ME)= delta_sn(4,ME) + 0.5*(FUNC2(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            
            delta_sp(1,ME)= delta_sp(1,ME) + 0.5*(FUNC3(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(2,ME)= delta_sp(2,ME) + 0.5*(FUNC3(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(3,ME)= delta_sp(3,ME) + 0.5*(FUNC3(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(4,ME)= delta_sp(4,ME) + 0.5*(FUNC3(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))


            else          
  
            delta_sm(1,ME)= delta_sm(1,ME) + 0.5*(FUNC1(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(2,ME)= delta_sm(2,ME) + 0.5*(FUNC1(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(3,ME)= delta_sm(3,ME) + 0.5*(FUNC1(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(4,ME)= delta_sm(4,ME) + 0.5*(FUNC1(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            
            delta_sm(1,NE)= delta_sm(1,NE) - 0.5*(FUNC1(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(2,NE)= delta_sm(2,NE) - 0.5*(FUNC1(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(3,NE)= delta_sm(3,NE) - 0.5*(FUNC1(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sm(4,NE)= delta_sm(4,NE) - 0.5*(FUNC1(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))


            delta_sn(1,ME)= delta_sn(1,ME) + 0.5*(FUNC2(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(2,ME)= delta_sn(2,ME) + 0.5*(FUNC2(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(3,ME)= delta_sn(3,ME) + 0.5*(FUNC2(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(4,ME)= delta_sn(4,ME) + 0.5*(FUNC2(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            
            delta_sn(1,NE)= delta_sn(1,NE) - 0.5*(FUNC2(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(2,NE)= delta_sn(2,NE) - 0.5*(FUNC2(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(3,NE)= delta_sn(3,NE) - 0.5*(FUNC2(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sn(4,NE)= delta_sn(4,NE) - 0.5*(FUNC2(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))


            delta_sp(1,ME)= delta_sp(1,ME) + 0.5*(FUNC3(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(2,ME)= delta_sp(2,ME) + 0.5*(FUNC3(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(3,ME)= delta_sp(3,ME) + 0.5*(FUNC3(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(4,ME)= delta_sp(4,ME) + 0.5*(FUNC3(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            
            delta_sp(1,NE)= delta_sp(1,NE) - 0.5*(FUNC3(1,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(2,NE)= delta_sp(2,NE) - 0.5*(FUNC3(2,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(3,NE)= delta_sp(3,NE) - 0.5*(FUNC3(3,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))
            delta_sp(4,NE)= delta_sp(4,NE) - 0.5*(FUNC3(4,J))*(X2(P1)*Y2(P2)-X2(P2)*Y2(P1))

            end if
                        
            Exit
         
        !Part 49:
        else    
            
            counter10=0
            if(X2(P2)==X2(P1))Then
                
                collisionX=X2(P2)
                counter10=counter10+1
                
            else
                
                M2=(Y2(P2)-Y2(P1))/(X2(P2)-X2(P1))
                
            End IF
            
            Do K=1,3
                
                h1=K;h2=K+1
                
                if(K==3)Then
                    h2=1
                    h1=3  
                End if
                
                
                if(X1(P_1(h2))==X1(P_1(h1)))Then
                
                    collisionX=X1(P_1(h2))
                    counter10=counter10+1
                    
                else
            
                    M_1(K)=(Y1(P_1(h2))-Y1(P_1(h1)))/(X1(P_1(h2))-X1(P_1(h1)))
                    if(counter10==1)then
                        
                        collisionY=Y1(P_1(h1))+M_1(K)*(collisionX-X1(P_1(h1)))
                    end if
                    
                End if
            
                If(counter10==2)Then
                
                    counter10=1
                    cycle
                else If(M2==M_1(K))Then
                
                    cycle
            
                else if(counter10==0)Then
                    
                    collisionX=(M2*X2(P1)-M_1(K)*X1(P_1(K))+Y1(P_1(K))-Y2(P1))/(M2-M_1(K))
                    collisionY=M2*(collisionX-X2(P1))+Y2(P1)
                    
                End If
                
              
                        
               
           F_s1 = (X1(P_1(1))-collisionX)*(Y1(P_1(2))-collisionY)-(X1(P_1(2))-collisionX)*(Y1(P_1(1))-collisionY)
           F_s2 = (X1(P_1(2))-collisionX)*(Y1(P_1(3))-collisionY)-(X1(P_1(3))-collisionX)*(Y1(P_1(2))-collisionY)
           F_s3 = (X1(P_1(3))-collisionX)*(Y1(P_1(1))-collisionY)-(X1(P_1(1))-collisionX)*(Y1(P_1(3))-collisionY)


           IF(abs(F_s1)<=1e-10)Then
               F_s1=0
           end if
           IF(abs(F_s2)<=1e-10)Then
               F_s2=0
           end if
           IF(abs(F_s3)<=1e-10)Then
               F_s3=0
           End IF
           
 
           Multi1=(Y2(P2)-collisionY)*(collisionY-Y2(P1))        !check            
           Multi2=(X2(P2)-collisionX)*(collisionX-X2(P1))  
           
           IF(abs(Multi1)<=1e-10 .AND. abs(Multi2)<=1e-10)Then
                Multi1=0
                Multi2=0
           else IF(abs(Multi1)<=1e-10)Then
                Multi1=0
           else IF(abs(Multi2)<=1e-10)Then
                Multi2=0
           End IF
           
           If((F_s1>=0. .And. F_s2>=0. .And. F_s3>=0. .And. Multi1>0) .OR. (F_s1>=0. .And. F_s2>=0. .And. F_s3>=0. .And. Multi2>0) .OR. (F_s1<=0. .And. F_s2<=0. .And. F_s3<=0. .And. Multi1>0) .OR. (F_s1<=0. .And. F_s2<=0. .And. F_s3<=0. .And. Multi2>0))Then
                           
                IF(counter6==1)Then
                    
                    If(NE==0)then
                    
                        delta_sm(1,ME)= delta_sm(1,ME) + 0.5*(FUNC1(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(2,ME)= delta_sm(2,ME) + 0.5*(FUNC1(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(3,ME)= delta_sm(3,ME) + 0.5*(FUNC1(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(4,ME)= delta_sm(4,ME) + 0.5*(FUNC1(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        
                        delta_sn(1,ME)= delta_sn(1,ME) + 0.5*(FUNC2(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(2,ME)= delta_sn(2,ME) + 0.5*(FUNC2(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(3,ME)= delta_sn(3,ME) + 0.5*(FUNC2(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(4,ME)= delta_sn(4,ME) + 0.5*(FUNC2(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))

                        delta_sp(1,ME)= delta_sp(1,ME) + 0.5*(FUNC3(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(2,ME)= delta_sp(2,ME) + 0.5*(FUNC3(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(3,ME)= delta_sp(3,ME) + 0.5*(FUNC3(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(4,ME)= delta_sp(4,ME) + 0.5*(FUNC3(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))

                
                    else
                                   
                        delta_sm(1,ME)= delta_sm(1,ME) + 0.5*(FUNC1(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(2,ME)= delta_sm(2,ME) + 0.5*(FUNC1(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(3,ME)= delta_sm(3,ME) + 0.5*(FUNC1(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(4,ME)= delta_sm(4,ME) + 0.5*(FUNC1(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                
                        delta_sm(1,NE)= delta_sm(1,NE) - 0.5*(FUNC1(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(2,NE)= delta_sm(2,NE) - 0.5*(FUNC1(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(3,NE)= delta_sm(3,NE) - 0.5*(FUNC1(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sm(4,NE)= delta_sm(4,NE) - 0.5*(FUNC1(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        
                        
                        delta_sn(1,ME)= delta_sn(1,ME) + 0.5*(FUNC2(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(2,ME)= delta_sn(2,ME) + 0.5*(FUNC2(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(3,ME)= delta_sn(3,ME) + 0.5*(FUNC2(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(4,ME)= delta_sn(4,ME) + 0.5*(FUNC2(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        
                        delta_sn(1,NE)= delta_sn(1,NE) - 0.5*(FUNC2(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(2,NE)= delta_sn(2,NE) - 0.5*(FUNC2(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(3,NE)= delta_sn(3,NE) - 0.5*(FUNC2(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sn(4,NE)= delta_sn(4,NE) - 0.5*(FUNC2(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))

                
                        delta_sp(1,ME)= delta_sp(1,ME) + 0.5*(FUNC3(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(2,ME)= delta_sp(2,ME) + 0.5*(FUNC3(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(3,ME)= delta_sp(3,ME) + 0.5*(FUNC3(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(4,ME)= delta_sp(4,ME) + 0.5*(FUNC3(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                    
                        delta_sp(1,NE)= delta_sp(1,NE) - 0.5*(FUNC3(1,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(2,NE)= delta_sp(2,NE) - 0.5*(FUNC3(2,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(3,NE)= delta_sp(3,NE) - 0.5*(FUNC3(3,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))
                        delta_sp(4,NE)= delta_sp(4,NE) - 0.5*(FUNC3(4,J))*(-X2(P1)*collisionY+collisionX*Y2(P1))

                    end if
                    
                else
                    
                    If(NE==0)then
                    
                        delta_sm(1,ME)= delta_sm(1,ME) + 0.5*(FUNC1(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(2,ME)= delta_sm(2,ME) + 0.5*(FUNC1(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(3,ME)= delta_sm(3,ME) + 0.5*(FUNC1(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(4,ME)= delta_sm(4,ME) + 0.5*(FUNC1(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                         
                        delta_sn(1,ME)= delta_sn(1,ME) + 0.5*(FUNC2(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(2,ME)= delta_sn(2,ME) + 0.5*(FUNC2(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(3,ME)= delta_sn(3,ME) + 0.5*(FUNC2(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(4,ME)= delta_sn(4,ME) + 0.5*(FUNC2(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))

                        delta_sp(1,ME)= delta_sp(1,ME) + 0.5*(FUNC3(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(2,ME)= delta_sp(2,ME) + 0.5*(FUNC3(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(3,ME)= delta_sp(3,ME) + 0.5*(FUNC3(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(4,ME)= delta_sp(4,ME) + 0.5*(FUNC3(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))

                        
                    else
                                   
                        delta_sm(1,ME)= delta_sm(1,ME) + 0.5*(FUNC1(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(2,ME)= delta_sm(2,ME) + 0.5*(FUNC1(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(3,ME)= delta_sm(3,ME) + 0.5*(FUNC1(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(4,ME)= delta_sm(4,ME) + 0.5*(FUNC1(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                
                        delta_sm(1,NE)= delta_sm(1,NE) - 0.5*(FUNC1(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(2,NE)= delta_sm(2,NE) - 0.5*(FUNC1(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(3,NE)= delta_sm(3,NE) - 0.5*(FUNC1(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sm(4,NE)= delta_sm(4,NE) - 0.5*(FUNC1(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        
                      
                        delta_sn(1,ME)= delta_sn(1,ME) + 0.5*(FUNC2(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(2,ME)= delta_sn(2,ME) + 0.5*(FUNC2(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(3,ME)= delta_sn(3,ME) + 0.5*(FUNC2(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(4,ME)= delta_sn(4,ME) + 0.5*(FUNC2(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        
                        delta_sn(1,NE)= delta_sn(1,NE) - 0.5*(FUNC2(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(2,NE)= delta_sn(2,NE) - 0.5*(FUNC2(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(3,NE)= delta_sn(3,NE) - 0.5*(FUNC2(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sn(4,NE)= delta_sn(4,NE) - 0.5*(FUNC2(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))

                        
                        delta_sp(1,ME)= delta_sp(1,ME) + 0.5*(FUNC3(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(2,ME)= delta_sp(2,ME) + 0.5*(FUNC3(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(3,ME)= delta_sp(3,ME) + 0.5*(FUNC3(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(4,ME)= delta_sp(4,ME) + 0.5*(FUNC3(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        
                        delta_sp(1,NE)= delta_sp(1,NE) - 0.5*(FUNC3(1,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(2,NE)= delta_sp(2,NE) - 0.5*(FUNC3(2,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(3,NE)= delta_sp(3,NE) - 0.5*(FUNC3(3,J))*(X2(P1)*collisionY-collisionX*Y2(P1))
                        delta_sp(4,NE)= delta_sp(4,NE) - 0.5*(FUNC3(4,J))*(X2(P1)*collisionY-collisionX*Y2(P1))

                end if
                    
                End IF
                 
                
            NP2p=NP2p+1    
            P1=NP2p
            X2(P1)=collisionX
            Y2(P1)=collisionY
            
            goto 20
            
           End IF
           
            
        End Do ! Do K
            
        End IF ! for intersection between edge and triangle of fine mesh
        
     
    !Part 50:
    else
        
      Cycle  
      
    End IF
        
    
    End Do ! Do J

    !Part 51:
    if((counter5==0 .AND. counter6==0) .OR. (counter13==0 .AND. counter5/=0))Then
          
          num1=P1
          P1=P2
          P2=num1
          counter6=1
          
          if(counter5/=0)then
              counter13=1
          end if
          
          goto 20
          
      end if
     
End Do ! Do I

!Part 52:
Do I=1,NC1
    
    PFUNC1(1,I)=FUNC1(1,I)
    PFUNC1(2,I)=FUNC1(2,I)
    PFUNC1(3,I)=FUNC1(3,I)
    PFUNC1(4,I)=FUNC1(4,I)
    
    PFUNC2(1,I)=FUNC2(1,I)
    PFUNC2(2,I)=FUNC2(2,I)
    PFUNC2(3,I)=FUNC2(3,I)
    PFUNC2(4,I)=FUNC2(4,I)
    
    PFUNC3(1,I)=FUNC3(1,I)
    PFUNC3(2,I)=FUNC3(2,I)
    PFUNC3(3,I)=FUNC3(3,I)
    PFUNC3(4,I)=FUNC3(4,I)
    
End Do

!Part 53:
if(NOL==1 .AND. FTC==1)then
NMF=2
else if(NOL==1 .AND. FTC==0)then
NMF=1
else if(NOL==2 .AND. FTC==1)then
NMF=3
else if(NOL==2 .AND. FTC==0)then
NMF=2
end if

!Part 54:
Call Read_2DMeshMG(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)
!Part 55:
Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
!Part 56:
Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)

!Part 57:
Do I=1,NC2
    
    FUNC1(1,I)=(delta_sm(1,I)+delta_s0(1,I))/A(I)
    FUNC1(2,I)=(delta_sm(2,I)+delta_s0(2,I))/A(I)
    FUNC1(3,I)=(delta_sm(3,I)+delta_s0(3,I))/A(I)
    FUNC1(4,I)=(delta_sm(4,I)+delta_s0(4,I))/A(I)

    FUNC2(1,I)=(delta_sn(1,I)+delta_s1(1,I))/A(I)
    FUNC2(2,I)=(delta_sn(2,I)+delta_s1(2,I))/A(I)
    FUNC2(3,I)=(delta_sn(3,I)+delta_s1(3,I))/A(I)
    FUNC2(4,I)=(delta_sn(4,I)+delta_s1(4,I))/A(I)

    FUNC3(1,I)=(delta_sp(1,I)+delta_s2(1,I))/A(I)
    FUNC3(2,I)=(delta_sp(2,I)+delta_s2(2,I))/A(I)
    FUNC3(3,I)=(delta_sp(3,I)+delta_s2(3,I))/A(I)
    FUNC3(4,I)=(delta_sp(4,I)+delta_s2(4,I))/A(I)

End Do

!Part 58:
Do I=1,NC2
    counter9=0
    P_2(1)=Corn2(1,I)
    P_2(2)=Corn2(2,I)
    P_2(3)=Corn2(3,I)
    
    Do J=1,NC1
        
    P_1(1)=Corn1(1,J)
    P_1(2)=Corn1(2,J)
    P_1(3)=Corn1(3,J)
        
    
    if((X2(P_2(1))==X1(P_1(1)) .AND. X2(P_2(2))==X1(P_1(2)) .AND. X2(P_2(3))==X1(P_1(3))) .OR. (X2(P_2(1))==X1(P_1(1)) .AND. X2(P_2(2))==X1(P_1(3)) .AND. X2(P_2(3))==X1(P_1(2))))then
        counter9=counter9+1
    Else if((X2(P_2(1))==X1(P_1(2)) .AND. X2(P_2(2))==X1(P_1(1)) .AND. X2(P_2(3))==X1(P_1(3))) .OR. (X2(P_2(1))==X1(P_1(2)) .AND. X2(P_2(2))==X1(P_1(3)) .AND. X2(P_2(3))==X1(P_1(1))))then
        counter9=counter9+1
    Else if((X2(P_2(1))==X1(P_1(3)) .AND. X2(P_2(2))==X1(P_1(1)) .AND. X2(P_2(3))==X1(P_1(2))) .OR. (X2(P_2(1))==X1(P_1(3)) .AND. X2(P_2(2))==X1(P_1(2)) .AND. X2(P_2(3))==X1(P_1(1))))then
        counter9=counter9+1
    end if
    
    
    if(counter9/=0)then
        
            FUNC1(1,I)=PFUNC1(1,J)
            FUNC1(2,I)=PFUNC1(2,J)
            FUNC1(3,I)=PFUNC1(3,J)
            FUNC1(4,I)=PFUNC1(4,J)

            FUNC2(1,I)=PFUNC2(1,J)
            FUNC2(2,I)=PFUNC2(2,J)
            FUNC2(3,I)=PFUNC2(3,J)
            FUNC2(4,I)=PFUNC2(4,J)

            FUNC3(1,I)=PFUNC3(1,J)
            FUNC3(2,I)=PFUNC3(2,J)
            FUNC3(3,I)=PFUNC3(3,J)
            FUNC3(4,I)=PFUNC3(4,J)

            Exit    
    end if
        
    End Do
End Do


Do I=1,NC2
    WNP1(1,I)=FUNC1(1,I);WNP1(2,I)=FUNC1(2,I);WNP1(3,I)=FUNC1(3,I);WNP1(4,I)=FUNC1(4,I)
    RES(1,I)=FUNC2(1,I);RES(2,I)=FUNC2(2,I);RES(3,I)=FUNC2(3,I);RES(4,I)=FUNC2(4,I)
    Error(1,I)=FUNC3(1,I);Error(2,I)=FUNC3(2,I);Error(3,I)=FUNC3(3,I);Error(4,I)=FUNC3(4,I)
end Do



end
!##########################################################################################################

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Read Edge Based Mesh From 'Mesh.Txt' File                            //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: November, 6, 2015                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Doc ID: MC5F003F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
Subroutine Read_2DMeshMG(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NMF
 Intent(Out  )::NP,NC,NF,NR,NFR,BC,IDS,X,Y

 Integer::Dim,I,J,J1,JJ,NP,NC,NF,NR,SFace,FaceType,MeshDim,NMF
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
!*******************************************************************************************
!Part 1:
 IF(NMF==1)Then
 Open(3,File='Mesh1.Txt')
Else IF(NMF==2)Then
 Open(3,File='Mesh2.Txt')
Else IF(NMF==3)Then
 Open(3,File='Mesh3.Txt')
End IF

!Part 2:
 Read(3,*) MeshDim
 IF(MeshDim/=2)Print*,'Please Check the Mesh File. It is not a 2D Mesh'

!Part 3:
 Read(3,*) NP    

!Part 4:
 Read(3,*) NC

!Part 5:
 Read(3,*) NF

!Part 6:
 Read(3,*) NR

!Part 7:
 Read(3,*)
 Do J=1,NR
    Read(3,*) NFR(J) , BC(J)
 End Do 

!Part 8:
 Read(3,*)
 Do J=1,NF  
    Read(3,*) FaceType,IDS(1,J),IDS(2,J),IDS(3,J),IDS(4,J)
 End Do

!Part 9:
 Read(3,*)        
 Do J=1,NP
    Read(3,*) X(J),Y(J)
 End Do

 Close(3)
!*******************************************************************************************
    End
!###########################################################################################

 
    


!*******************************************************************************************
 Subroutine EbasedToCbased(Dim,NC,NF,IDS,Corn)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NF,IDS 
 Intent(Out  )::Corn

 Integer::Dim,J,I,NC,NF,ME,NE,P1,P2 
 Integer,Dimension(1:4,1:Dim)::IDS,Corn
 Integer,Dimension(1:Dim)::NCorn
!*******************************************************************************************

!Part 6:
 Do J=1,NC
	Corn(1,J)=0
	Corn(2,J)=0
	Corn(3,J)=0
	Corn(4,J)=0
    NCorn(J)=0
 End Do

!Part 6:
 Do J=1,NF
    ME = IDS(1,J)
    NE = IDS(2,J)
    P1 = IDS(3,J)
	P2 = IDS(4,J)

	Do I=1,NCorn(ME)
	   IF( P1==Corn(I,ME) ) goto 1
    End Do
    NCorn(ME) = NCorn(ME) + 1
	Corn( NCorn(ME),ME ) = P1

1	Do I=1,NCorn(ME)
	   IF( P2==Corn(I,ME) ) goto 2
    End Do
    NCorn(ME) = NCorn(ME) + 1
	Corn( NCorn(ME), ME) = P2

2   IF(NE==0)goto 4

	Do I=1,NCorn(NE)
	   IF( P1==Corn(I,NE) ) goto 3
    End Do
    NCorn(NE) = NCorn(NE) + 1
	Corn(NCorn(NE), NE ) = P1

3	Do I=1,NCorn(NE)
	   IF( P2==Corn(I,NE) ) goto 4
    End Do
    NCorn(NE) = NCorn(NE) + 1
	Corn(NCorn(NE), NE ) = P2  
 
4 End Do
!*******************************************************************************************
 End
!###########################################################################################
 Subroutine Read_2DMesh1(Dim,NP1,NC1,NF_1,NR1,NFR1,BC1,IDS1,X1,Y1,NOL,FTC)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NOL,FTC
 Intent(Out  )::NP1,NC1,NF_1,NR1,NFR1,BC1,IDS1,X1,Y1

 Integer::Dim,I,J,J1,JJ,NP1,NC1,NF_1,NR1,SFace,FaceType,MeshDim,NOL,FTC
 Integer,Dimension(1:100)::NFR1,BC1
 Integer,Dimension(1:4,1:Dim)::IDS1
 Real(8),Dimension(1:Dim)::X1,Y1
!*******************************************************************************************
!Part 1:
 IF(NOL==1 .AND. FTC==1)then
   Open(1,File='Mesh1.Txt')
 else IF(NOL==1 .AND. FTC==0)then
   Open(1,File='Mesh2.Txt')
 else if(NOL==2 .AND. FTC==1)then
   Open(1,File='Mesh2.Txt')
 else if(NOL==2 .AND. FTC==0)then
   Open(1,File='Mesh3.Txt')
 end if

!Part 2:
 Read(1,*) MeshDim
 IF(MeshDim/=2)Print*,'Please Check the Mesh File. It is not a 2D Mesh'

!Part 3:
 Read(1,*) NP1    

!Part 4:
 Read(1,*) NC1

!Part 5:
 Read(1,*) NF_1

!Part 6:
 Read(1,*) NR1

!Part 7:
 Read(1,*)
 Do J=1,NR1
    Read(1,*) NFR1(J) , BC1(J)
 End Do 

!Part 8:
 Read(1,*)
 Do J=1,NF_1  
    Read(1,*) FaceType,IDS1(1,J),IDS1(2,J),IDS1(3,J),IDS1(4,J)
 End Do

!Part 9:
 Read(1,*)        
 Do J=1,NP1
    Read(1,*) X1(J),Y1(J)
 End Do

 Close(1)
!*******************************************************************************************
End
!###########################################################################################
Subroutine Read_2DMesh2(Dim,NP2,NC2,NF_2,NR2,NFR2,BC2,IDS2,X2,Y2,NOL,FTC)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NOL,FTC
 Intent(Out  )::NP2,NC2,NF_2,NR2,NFR2,BC2,IDS2,X2,Y2

 Integer::Dim,I,J,J1,JJ,NP2,NC2,NF_2,NR2,SFace,FaceType,MeshDim,NOL,FTC
 Integer,Dimension(1:100)::NFR2,BC2
 Integer,Dimension(1:4,1:Dim)::IDS2
 Real(8),Dimension(1:Dim)::X2,Y2
!*******************************************************************************************
!Part 1:
 IF(NOL==1 .AND. FTC==1)then
   Open(2,File='Mesh2.Txt')
 else IF(NOL==1 .AND. FTC==0)then
   Open(2,File='Mesh1.Txt')
 else if(NOL==2 .AND. FTC==1)then
   Open(2,File='Mesh3.Txt')
 else if(NOL==2 .AND. FTC==0)then
   Open(2,File='Mesh2.Txt')
 end if

!Part 2:
 Read(2,*) MeshDim
 IF(MeshDim/=2)Print*,'Please Check the Mesh File. It is not a 2D Mesh'

!Part 3:
 Read(2,*) NP2    

!Part 4:
 Read(2,*) NC2

!Part 5:
 Read(2,*) NF_2

!Part 6:
 Read(2,*) NR2

!Part 7:
 Read(2,*)
 Do J=1,NR2
    Read(2,*) NFR2(J) , BC2(J)
 End Do 

!Part 8:
 Read(2,*)
 Do J=1,NF_2  
    Read(2,*) FaceType,IDS2(1,J),IDS2(2,J),IDS2(3,J),IDS2(4,J)
 End Do

!Part 9:
 Read(2,*)        
 Do J=1,NP2
    Read(2,*) X2(J),Y2(J)
 End Do

 Close(2)
!*******************************************************************************************
End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Find Index of Points Constructing the cell                           //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Doc ID: MC5F088F1                                                                    //!
!//                                                                                      //!
!// This Program is Available Through the Website: www.MarketCode.ir                     //!
!// It May be Copied, Modified, and Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Edge_To_Cell(Dim,NF,NC,IDS,Corn)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NF,NC,IDS
 Intent(Inout)::Corn
 
 Integer::Dim,NF,NC,ME,NE,I,J,J1,J2,E,E1,E2,E3,P1_E2,P2_E1,P
 Integer,Dimension(1:4, 1:Dim)::IDS
 Integer,Dimension(1:4, 1:Dim)::CELL_EDGE,Corn
 Integer,Dimension(1:Dim)::NCELL_EDGE
!*******************************************************************************************
!Part 1:
 Do I=1,NC
    NCELL_EDGE(I) = 0   
	Corn(1,I)     = 0  
	Corn(2,I)     = 0  
	Corn(3,I)     = 0  
	Corn(4,I)     = 0
 End Do
 
!Part 2:
 Do I=1,NF

   !Part 3:
    ME = IDS(1, I)
    NE = IDS(2, I)
	
   !Part 4:
	NCELL_EDGE(ME) = NCELL_EDGE(ME) + 1
    CELL_EDGE(NCELL_EDGE(ME),ME)=I

   !Part 5:
    IF(NE/=0)Then
	 NCELL_EDGE(NE) = NCELL_EDGE(NE) + 1
     CELL_EDGE(NCELL_EDGE(NE),NE)=-I
    EndIF

 End Do

!Part 6:
 Do I=1,NC

   !Part 7:
    Do J1=1,NCELL_EDGE(I)
       E1 = CELL_EDGE(J1,I)

      !Part 8:
       IF( E1>0 )Then
        P2_E1 = IDS(4,E1)
       Else
        P2_E1 = IDS(3,-E1)
	   EndIF

      !Part 9:
	   Do J2=J1+1,NCELL_EDGE(I)
          E2 = CELL_EDGE(J2,I)

         !Part 10:
          IF( E2>0 )Then
           P1_E2 = IDS(3,E2)
          Else
           P1_E2 = IDS(4,-E2)
          EndIF

         !Part 11:
          IF( P2_E1==P1_E2 )Then
	       E                 = CELL_EDGE(J1+1,I)
	       CELL_EDGE(J1+1,I) = CELL_EDGE(J2,I)
           CELL_EDGE(J2,I)   = E
          EndIF

       End Do

    End Do

 End Do 

!Part 12:
 Do I=1,NC
    Do J=1,NCELL_EDGE(I)

       E = CELL_EDGE(J,I)

       IF( E>0 )Then
        P = IDS(3,E)
       Else
        P = IDS(4,-E)
       EndIF

       Corn(J,I) = P

    End Do
 End Do
!*******************************************************************************************
 End
!###########################################################################################

!###########################################################################################
 Subroutine Mesh_Info(Dim,IDS,X,Y,CornC,NeibC,NC,NF)
Implicit None
!*******************************************************************************************
Intent(In   )::Dim,NC,NF,IDS
Intent(Out  )::CornC,NeibC

Integer::Dim  
Integer::I,J,ME,NE,P1,P2,counter,NC,NF
Integer,Dimension(1:3,1:Dim)::CornC,NeibC
Real(8),Dimension(1:Dim)::X,Y
Integer,Dimension(1:4,1:Dim)::IDS
!*******************************************************************************************    
Do I=1,NC
    counter=0
    Do J=1,NF
        ME=IDS(1,J)
        NE=IDS(2,J)
        P1=IDS(3,J)
        P2=IDS(4,J)
        IF (I==ME .AND. counter==0) Then
            CornC(1,I)=P1
            CornC(2,I)=P2
            NeibC(3,I)=NE
            counter=counter+1
            cycle
        End IF
            
        IF (I==NE .AND. counter==0) Then
            CornC(1,I)=P1
            CornC(2,I)=P2
            NeibC(3,I)=ME
            counter=counter+1
            cycle
        End IF
            
        IF (I==ME) Then
                
            IF((CornC(3,I)==P1 .AND. CornC(2,I)==P2) .OR. (CornC(3,I)==P2 .AND. CornC(2,I)==P1)) Then
                counter=counter+1
                NeibC(1,I)=NE
                cycle
            End IF
                
            IF((CornC(3,I)==P1 .AND. CornC(1,I)==P2) .OR. (CornC(3,I)==P2 .AND. CornC(1,I)==P1)) Then
                counter=counter+1
                NeibC(2,I)=NE
                cycle
            End IF
                
            IF (CornC(1,I)==P1) Then
                counter=counter+1;
                CornC(3,I)=P2;
                NeibC(2,I)=NE;
                Cycle;
            End IF
                
            IF (CornC(1,I)==P2) Then
                counter=counter+1;
                CornC(3,I)=P1;
                NeibC(2,I)=NE;
                Cycle;
            End IF
                
            IF (CornC(2,I)==P1) Then
                counter=counter+1;
                CornC(3,I)=P2;
                NeibC(1,I)=NE;
                Cycle;
            End IF
                
            IF (CornC(2,I)==P2) Then
                counter=counter+1;
                CornC(3,I)=P1;
                NeibC(1,I)=NE;
                Cycle;
            End IF
                
        End IF
            
            
        IF (I==NE) Then
                
            IF((CornC(3,I)==P1 .AND. CornC(2,I)==P2) .OR. (CornC(3,I)==P2 .AND. CornC(2,I)==P1)) Then
                counter=counter+1
                NeibC(1,I)=ME
                cycle
            End IF
                
            IF((CornC(3,I)==P1 .AND. CornC(1,I)==P2) .OR. (CornC(3,I)==P2 .AND. CornC(1,I)==P1)) Then
                counter=counter+1
                NeibC(2,I)=ME
                cycle
            End IF
                
            IF (CornC(1,I)==P1) Then
                counter=counter+1;
                CornC(3,I)=P2;
                NeibC(2,I)=ME;
                Cycle;
            End IF
                
            IF (CornC(1,I)==P2) Then
                counter=counter+1;
                CornC(3,I)=P1;
                NeibC(2,I)=ME;
                Cycle;
            End IF
                
            IF (CornC(2,I)==P1) Then
                counter=counter+1;
                CornC(3,I)=P2;
                NeibC(1,I)=ME;
                Cycle;
            End IF
                
            IF (CornC(2,I)==P2) Then
                counter=counter+1;
                CornC(3,I)=P1;
                NeibC(1,I)=ME;
                Cycle;
            End IF
        End IF
                    
    End Do
End Do
           
!Do I=1,NC
!    print*,CornC(1,I),CornC(2,I),CornC(3,I),NeibC(1,I),NeibC(2,I),NeibC(3,I)
!End Do

End