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
!// Date: Feb., 10, 2018                                                                   //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!// Developed by: M. Valadkhani, Mechanical Eng., Amirkabir University of Technology       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine GridFlux_3D(Dim,NF,IDS,X,Xo,Y,Yo,Z,Zo,FaceType,delT,GF,Face_Velocity)
 Implicit None
!********************************************************************************************* 
Intent(IN)::Dim,NF,IDS,X,Xo,Y,Yo,Z,Zo,FaceType,delT
Intent(OUT)::GF,Face_Velocity

Integer::Dim,NF,I,P1,P2,P3,P4
Real(8)::Xa,Ya,Za,Xb,Yb,Zb,Xc,Yc,Zc,Xd,Yd,Zd,NI1,NJ1,NK1,NI2,NJ2,NK2,delT
Real(8)::Xba_n,Yba_n,Zba_n,Xba_o,Yba_o,Zba_o,Xca_n,Yca_n,Zca_n,Xca_o,Yca_o,Zca_o,Xda_n,Yda_n,Zda_n,Xda_o,Yda_o,Zda_o
Real(8),Dimension(1:Dim)::X,Xo,Y,Yo,Z,Zo,DT,GF
Integer,Dimension(1:6,1:Dim)::IDS
Integer,Dimension(1:Dim)::FaceType
Real(8),Dimension(1:3,1:Dim)::Face_Velocity
!********************************************************************************************* 
!Part 1:
 DO I=1,NF
     If(FaceType(I)==3) Then
         P1 = IDS(3,I)
         P2 = IDS(4,I)
         P3 = IDS(5,I)
         
         Xa=X(P1)-Xo(P1) ; Ya=Y(P1)-Yo(P1) ; Za=Z(P1)-Zo(P1)
         Xb=X(P2)-Xo(P2) ; Yb=Y(P2)-Yo(P2) ; Zb=Z(P2)-Zo(P2)
         Xc=X(P3)-Xo(P3) ; Yc=Y(P3)-Yo(P3) ; Zc=Z(P3)-Zo(P3)
         
         Face_Velocity(1,I) = (Xa+Xb+Xc) / (3.*DelT)
         Face_Velocity(2,I) = (Ya+Yb+Yc) / (3.*DelT)
         Face_Velocity(3,I) = (Za+Zb+Zc) / (3.*DelT)
         
         Xba_n = X(P2)-X(P1)   ; Yba_n = Y(P2)-Y(P1)   ; Zba_n = Z(P2)-Z(P1)
         Xba_o = Xo(P2)-Xo(P1) ; Yba_o = Yo(P2)-Yo(P1) ; Zba_o = Zo(P2)-Zo(P1)
         
         Xca_n = X(P3)-X(P1)   ; Yca_n = Y(P3)-Y(P1)   ; Zca_n = Z(P3)-Z(P1)
         Xca_o = Xo(P3)-Xo(P1) ; Yca_o = Yo(P3)-Yo(P1) ; Zca_o = Zo(P3)-Zo(P1)
         
         NI1=(1./3.)*(Yba_n*Zca_n-Zba_n*Yca_n) + (1./6.)*(Yba_n*Zca_o+Yba_o*Zca_n-Zba_n*Yca_o-Zba_o*Yca_n) + (1./3.)*(Yba_o*Zca_o-Zba_o*Yca_o)
         NJ1=(1./3.)*(Zba_n*Xca_n-Xba_n*Zca_n) + (1./6.)*(Zba_n*Xca_o+Zba_o*Xca_n-Xba_n*Zca_o-Xba_o*Zca_n) + (1./3.)*(Zba_o*Xca_o-Xba_o*Zca_o)
         NK1=(1./3.)*(Xba_n*Yca_n-Yba_n*Xca_n) + (1./6.)*(Xba_n*Yca_o+Xba_o*Yca_n-Yba_n*Xca_o-Yba_o*Xca_n) + (1./3.)*(Xba_o*Yca_o-Yba_o*Xca_o)
         
         GF(I)=( Face_Velocity(1,I)*NI1 +  Face_Velocity(2,I)*NJ1 +  Face_Velocity(3,I)*NK1 )/2.0
         
     Else If(FaceType(I)==4) Then
         P1 = IDS(3,I)
         P2 = IDS(4,I)
         P3 = IDS(5,I)
         P4 = IDS(6,I)
         
         Xa=X(P1)-Xo(P1) ; Ya=Y(P1)-Yo(P1) ; Za=Z(P1)-Zo(P1)
         Xb=X(P2)-Xo(P2) ; Yb=Y(P2)-Yo(P2) ; Zb=Z(P2)-Zo(P2)
         Xc=X(P3)-Xo(P3) ; Yc=Y(P3)-Yo(P3) ; Zc=Z(P3)-Zo(P3)
         Xd=X(P4)-Xo(P4) ; Yd=Y(P4)-Yo(P4) ; Zd=Z(P4)-Zo(P4)
         
         Face_Velocity(1,I) = (Xa+Xb+Xc+Xd) / (4.*DelT)
         Face_Velocity(2,I) = (Ya+Yb+Yc+Yd) / (4.*DelT)
         Face_Velocity(3,I) = (Za+Zb+Zc+Zd) / (4.*DelT)
         
         Xba_n = X(P2)-X(P1)   ; Yba_n = Y(P2)-Y(P1)   ; Zba_n = Z(P2)-Z(P1)
         Xba_o = Xo(P2)-Xo(P1) ; Yba_o = Yo(P2)-Yo(P1) ; Zba_o = Zo(P2)-Zo(P1)
         
         Xca_n = X(P3)-X(P1)   ; Yca_n = Y(P3)-Y(P1)   ; Zca_n = Z(P3)-Z(P1)
         Xca_o = Xo(P3)-Xo(P1) ; Yca_o = Yo(P3)-Yo(P1) ; Zca_o = Zo(P3)-Zo(P1)
         
         NI1=(1./3.)*(Yba_n*Zca_n-Zba_n*Yca_n) + (1./6.)*(Yba_n*Zca_o+Yba_o*Zca_n-Zba_n*Yca_o-Zba_o*Yca_n) + (1./3.)*(Yba_o*Zca_o-Zba_o*Yca_o)
         NJ1=(1./3.)*(Zba_n*Xca_n-Xba_n*Zca_n) + (1./6.)*(Zba_n*Xca_o+Zba_o*Xca_n-Xba_n*Zca_o-Xba_o*Zca_n) + (1./3.)*(Zba_o*Xca_o-Xba_o*Zca_o)
         NK1=(1./3.)*(Xba_n*Yca_n-Yba_n*Xca_n) + (1./6.)*(Xba_n*Yca_o+Xba_o*Yca_n-Yba_n*Xca_o-Yba_o*Xca_n) + (1./3.)*(Xba_o*Yca_o-Yba_o*Xca_o)
  
         GF(I)=( Face_Velocity(1,I)*NI1 +  Face_Velocity(2,I)*NJ1 +  Face_Velocity(3,I)*NK1 )/1.0
     
     End If
 End Do
!*********************************************************************************************
 End
!########################################################################################### 
