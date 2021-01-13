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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!// Developed by: M. Valadkhani, Mechanical Eng., Amirkabir University of Technology       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine GridFlux(Dim,NF1,NF2,NFW1,NFW2,IDS,X,Xo,Y,Yo,delT,GF,Face_Velocity)
 Implicit None
!********************************************************************************************* 
Intent(IN)::Dim,NF1,NF2,NFW1,NFW2,IDS,X,Xo,Y,Yo,delT
Intent(OUT)::GF,Face_Velocity

Integer::Dim,NF1,NF2,NFW1,NFW2,I,P1,P2
Real(8)::Xa,Ya,Xb,Yb,delT
Real(8),Dimension(1:Dim)::X,Xo,Y,Yo,DT,GF
Integer,Dimension(1:4,1:Dim)::IDS
Real(8),Dimension(1:2,1:Dim)::Face_Velocity
!********************************************************************************************* 

!Part 1:
 DO I=NFW1+1,NFW2

   !part 2
    P1 = IDS(3,I)
    P2 = IDS(4,I)

    !part 3
    Xa = ((X(P2)+X(P1))/2.0)-((Xo(P2)+Xo(P1))/2.0)
    Ya = ((Y(P2)+Y(P1))/2.0)-((Yo(P2)+Yo(P1))/2.0)
    
    !part 4
    Xb = ((X(P2)+Xo(P2))/2.0)-((X(P1)+Xo(P1))/2.0)
    Yb = ((Y(P2)+Yo(P2))/2.0)-((Y(P1)+Yo(P1))/2.0)
    
    !part 5
    GF(I)=((Xa*Yb)-(Xb*Ya))/delT

    !part 6
    Face_Velocity(1,I)=Xa/DelT
    Face_Velocity(2,I)=Ya/DelT
 End Do


!Part 12:
 DO I=NF1+1,NF2
 
    !part 13
    P1 = IDS(3,I)
    P2 = IDS(4,I)

   !part 14
    Xa = ((X(P2)+X(P1))/2.0)-((Xo(P2)+Xo(P1))/2.0)
    Ya = ((Y(P2)+Y(P1))/2.0)-((Yo(P2)+Yo(P1))/2.0)
    
    !part 15
    Xb = ((X(P2)+Xo(P2))/2.0)-((X(P1)+Xo(P1))/2.0)
    Yb = ((Y(P2)+Yo(P2))/2.0)-((Y(P1)+Yo(P1))/2.0)
    
    !part 16
    GF(I)=((Xa*Yb)-(Xb*Ya))/delT
    
    Face_Velocity(1,I)=Xa/DelT
    Face_Velocity(2,I)=Ya/DelT
   
 End Do

!*********************************************************************************************
 End
!########################################################################################### 
