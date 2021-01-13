!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//       /////////////       ////////////////    ////////     //////    ////////////////  //!
!//       /////////////      ////////////////    //////////   //////    ////////////////   //!
!//      /////    /////     //////    //////    //////////// //////    /////               //!
!//     /////    //////    ////////////////    ///////////////////    ////////////////     //!
!//    /////    //////    ////////////////    ////// ////////////               /////      //!
!//   ///////////////    //////    //////    //////   //////////    ////////////////       //!
!// ///////////////     //////    //////    //////     ////////    ////////////////        //!
!//    Developer            Assistant    in      Numerical             Sciences            //!
!//----------------------------------------------------------------------------------------//!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace eng. Amirkabir University of Technology      //!
!// Date: May., 15, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: M. Valadkhani, Mechanical Eng., Amirkabir University of Technology       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Defin_BoundPoint_Displac(Dim,NFW1,NFW2,IDS,X,Y,Time,dT,Test_Case,Minf,alfa,Omega,Delx,Dely)
 Implicit None
!********************************************************************************************* 
 Intent(In ):: Dim,NFW1,NFW2,IDS,X,Y,Time,dT,Test_Case,Minf
 Intent(Out )::omega,Delx,Dely
 Intent(InOut )::alfa
 
 Integer::Dim,I,J,NFW1,NFW2,P1,Test_Case
 Real(8)::PI,Amplitude,Omega,Time,dT,k,Xo,Yo,X_Vec,Y_Vec,X_Rot,Y_Rot,&
          Minf,alfa,Del_Teta,TetaMax,Alpha_Dot,DeltaY
 Real(8),Dimension(1:Dim)::Delx,Dely,X,Y
 Integer,Dimension(1:4,1:Dim)::IDS
!********************************************************************************************* 
!Part 1:
 Delx(1:Dim)=0.0
 Dely(1:Dim)=0.0
 
!Part 2:
 PI=4.*Atan(1.0)
 
!Part 3:
!Select Test Case
If (Test_Case==1) Then
    !===========================Pitching==========================NACA 0012========
    !Minf=0.6 and Alfa_Mean=2.89  difined at setting.txt file    (Re=4.8E+06)    T=64.802 s
     
    !Part 4:
    !Note:Center of rotation for each mesh may vary
    Xo = 0.25
    Yo = 4.00 ! 0.0
    
    !Part 5:
    TetaMax = 2.41 * (PI/180.0)
    k = 0.0808
    omega = 2.0*k*Minf
    
    !Part 6:
    Del_Teta = -1 * TetaMax * (sin(omega*Time)-sin(omega*(Time-dT)))
    
    !Part 7:
    DO J=NFW1+1,NFW2
        !Part 8:
        P1  = IDS(3,J)
        
        !Part 9:
        X_Vec = X(P1) - Xo
        Y_Vec = Y(P1) - Yo
        
        !Part 10:
        X_Rot = Xo + X_Vec*cos(Del_Teta) - Y_Vec*sin(Del_Teta)
        Y_Rot = Yo + X_Vec*sin(Del_Teta) + Y_Vec*cos(Del_Teta)
        
        !Part 11:
        Delx(P1)= X_Rot - X(P1)
        Dely(P1)= Y_Rot - Y(P1)
	 
    End Do
    
    !Part 12:
    Alfa=alfa - Del_Teta * (180./PI)
    
    !Part 13:
    If ( Time > 130.0 ) stop
    !=============================================================
!Part 3:
Else If (Test_Case==2) Then
    !===========================Pitching============================NACA 0012========
    !Minf=0.755 and Alfa_Mean=0.016  difined at setting.txt file    (Re=5.5E+06)    T=51.1185 s
     
    !Part 4:
    !Note:Center of rotation for each mesh may vary
    Xo = 0.25
    Yo = 4.00 ! 0.0
    
    !Part 5:
    TetaMax = 2.51 * (PI/180.0)
    k = 0.0814
    omega = 2.0*k*Minf
    
    !Part 6:
    Del_Teta = -1 * TetaMax * (sin(omega*Time)-sin(omega*(Time-dT)))
    
    !Part 7:
    DO J=NFW1+1,NFW2
        !Part 8:
        P1  = IDS(3,J)
        
        !Part 9:
        X_Vec = X(P1) - Xo
        Y_Vec = Y(P1) - Yo
        
        !Part 10:
        X_Rot = Xo + X_Vec*cos(Del_Teta) - Y_Vec*sin(Del_Teta)
        Y_Rot = Yo + X_Vec*sin(Del_Teta) + Y_Vec*cos(Del_Teta)
        
        !Part 11:
        Delx(P1)= X_Rot - X(P1)
        Dely(P1)= Y_Rot - Y(P1)
	 
    End Do
    
    !Part 12:
    Alfa=alfa - Del_Teta * (180./PI)
    
    !Part 13:
    If ( Time > 120.0 ) stop
    !=============================================================
!Part 3:
Else If (Test_Case==3) Then
    !===========================Pitching==========================NACA 64A010========
    !Minf=0.80 and Alfa_Mean=0.0  difined at setting.txt file    (Re=12.5E+06)  T=19.4405 s
    
    !Part 4:
    !Note:Center of rotation for each mesh may vary
    Xo = 0.25
    Yo = 4.00 ! 0.0
    
    !Part 5:
    TetaMax = 1.01 * (PI/180.0)
    k = 0.2020
    omega = 2.0*k*Minf
    
    !Part 6:
    Del_Teta = -1 * TetaMax * (sin(omega*Time)-sin(omega*(Time-dT)))
    
    !Part 7:
    DO J=NFW1+1,NFW2
        !Part 8:
        P1  = IDS(3,J)
        
        !Part 9:
        X_Vec = X(P1) - Xo
        Y_Vec = Y(P1) - Yo
        
        !Part 10:
        X_Rot = Xo + X_Vec*cos(Del_Teta) - Y_Vec*sin(Del_Teta)
        Y_Rot = Yo + X_Vec*sin(Del_Teta) + Y_Vec*cos(Del_Teta)
        
        !Part 11:
        Delx(P1)= X_Rot - X(P1)
        Dely(P1)= Y_Rot - Y(P1)
        
    End Do
    
    !Part 12:
    Alfa=alfa - Del_Teta * (180./PI)
    
    !Part 13:
    If ( Time > 60.0 ) stop
    !=============================================================
!Part 3:
Else If (Test_Case==4) Then
    !===========================Ramp_Motion==========================NACA 0012========
    !Minf=0.30 and Alfa_initial=-0.03  and  Alfa_Final=15.54  difined at setting.txt file    (Re=12.5E+06)   t=15.57/1280.=0.01216 s
    
    !Part 4:
    !Note:Center of rotation for each mesh may vary
    Xo = 0.25
    Yo = 4.00 ! 0.0
    
    !Part 5:
    TetaMax = 15.54 * (PI/180.0)
    Alpha_Dot= 0.43745 * (PI/180.0)  !Alfa_Dot_Star=ALfa_Dot*(cord_Length/speed_of_Sound)=1280*(0.1016/297.284)=0.43745    t=35.5926 s
    omega = 0.3530   ! omega=2*PI/T=0.3530     T= 17.8 
    
    !Part 6:
    Del_Teta = -1 * (Alpha_Dot * dT) 
    
    !Part 7:
    DO J=NFW1+1,NFW2
        !Part 8:
        P1  = IDS(3,J)
        
        !Part 9:
        X_Vec = X(P1) - Xo
        Y_Vec = Y(P1) - Yo
        
        !Part 10:
        X_Rot = Xo + X_Vec*cos(Del_Teta) - Y_Vec*sin(Del_Teta)
        Y_Rot = Yo + X_Vec*sin(Del_Teta) + Y_Vec*cos(Del_Teta)
        
        !Part 11:
        Delx(P1)= X_Rot - X(P1)
        Dely(P1)= Y_Rot - Y(P1)
        
    End Do
    
    !Part 12:
    Alfa=alfa - Del_Teta * (180./PI)
    
    !Part 13:
    If( Time > 35.6 )  stop
    !=============================================================
!Part 3:
Else If (Test_Case==5) Then
    !===========================Plunging==========================NACA 0012========
    !Minf=0.10 and Alfa_Mean=0.00  difined at setting.txt file    (Re=4.E+06)  T=7.854 s
    
    !Part 5:
    Amplitude=0.0125
    k=4.0
    omega=2.0*k*Minf
    
    !Part 6:
    DeltaY=Amplitude*(sin(omega*Time)-sin(omega*(Time-dT)))
    
    !Part 7:
    Do J=NFW1+1,NFW2
        !Part 8:
        P1  = IDS(3,J)
        
        !Part 11:
       !Delx(P1)=0.0
        Dely(P1)=DeltaY
    End Do
    
    !Part 12:
    Alfa = - Atan(Amplitude*Omega*Cos(omega*Time)/Minf)*(180./PI)  
    
    !Part 13:
    If ( Time > 45.0 ) stop
    !=============================================================
!Part 3:
Else If (Test_Case==6) Then
    !===========================Plunging==========================NACA 0012========
    !Minf=0.05 and Alfa_Mean=0.00  difined at setting.txt file    (Re=2.0E+04)  T=15.0708 s
    
    !Part 5:
    Amplitude=0.0125
    k=4.0
    omega=2.0*k*Minf
    
    !Part 6:
    DeltaY=Amplitude*(sin(omega*Time)-sin(omega*(Time-dT)))
    
    !Part 7:
    Do J=NFW1+1,NFW2
        !Part 7:
        P1  = IDS(3,J)
        
        !Part 11:
       !Delx(P1)=0.0
        Dely(P1)=DeltaY
    End Do
    
    !Part 12:
    Alfa= - Atan(Amplitude*Omega*Cos(omega*Time)/Minf)*(180./PI) 
    
    !Part 13:
    If ( Time > 60.0 ) stop
    !=============================================================
!Part 3:
Else If (Test_Case==7) Then
    !===========================Plunging==========================NACA 0012========
    !Minf=0.30 and Alfa_Mean=0.00  difined at setting.txt file    (Re=3.0E+06)  T=6.98 s
    
    !Part 5:
    Amplitude=0.10
    k=1.5
    omega=2.0*k*Minf
    
    !Part 6:
    DeltaY=-1*Amplitude*(cos(omega*Time)-cos(omega*(Time-dT)))
    
    !Part 7:
    Do J=NFW1+1,NFW2
        !Part 8:
        P1  = IDS(3,J)
        
        !Part 11:
       !Delx(P1)=0.0
        Dely(P1)=DeltaY
    End Do
    
    !Part 12:
    Alfa= - Atan(Amplitude*Omega*Sin(omega*Time)/Minf)*(180./PI)
    
    !Part 13:
    If ( Time > 69.82 ) stop
    !============================================================= 

End If

!********************************************************************************************* 
 End 
!###########################################################################################
