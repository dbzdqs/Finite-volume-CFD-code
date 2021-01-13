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
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Update_Mesh(Dim,DimU,DimL,NUP,IBP_UP,Ynew_Up,Y_Up,NLP,IBP_Lw,Ynew_Lw,Y_Lw,NBP,NP,NC,NF,&
                        NF1,NF2,IBP,X,Y,IDS,DelX,DelY,Xc,Yc,NX,NY,DA,Vol)
 Implicit None
!************************************************************************************* 
 Intent(In    )Dim,DimU,DimL,NF1,NF2,NBP,NP,IBP,IDS,NUP,IBP_UP,Ynew_Up,Y_Up,NLP,IBP_Lw,Ynew_Lw,Y_Lw
 Intent(InOut )DelX,DelY
 Intent(Out   )X,Y,Xc,Yc,NX,NY,DA,Vol

 Integer Dim,DimU,DimL,NF1,NF2,I,NBP,NP,J,Move,NC,NF,NUP,NLP,P1,P2
 Real(8),Dimension(Dim)::X,Y,DelX,DelY,Xc,Yc,NX,NY,DA,Vol
 Integer,Dimension(1:Dim)::IBP
 Integer,Dimension(1:4,1:Dim)::IDS 

 Integer,Dimension(1:DimL)::IBP_Lw
 Integer,Dimension(1:DimU)::IBP_Up
 Real(8),Dimension(1:DimL)::Y_Lw,Ynew_Lw
 Real(8),Dimension(1:DimU)::Y_Up,Ynew_Up
!************************************************************************************* 
!Part 1:
 Do I=1,NP
    DelX(I) = 0.0 
    DelY(I) = 0.0 
 End Do

!Part 2:
 Do I=1,NUP
    P1 = IBP_UP(I)
    DelY(P1) = Ynew_Up(I) - Y_Up(I) 
 End Do

 Do I=1,NLP
	P2 = IBP_Lw(I)
    DelY(P2) = Ynew_Lw(I) - Y_Lw(I) 
 End Do

!Part 3:
 Call RBF_Moving_Mesh(Dim,NBP,NP,IBP,X,Y,DelX,DelY)

!Part 4:
 Do J=1,NP
    X(J)=X(J)+DelX(J)
    Y(J)=Y(J)+DelY(J)
 End Do

!Part 5:
 Call MoveCheckEbased2D(Dim,NF,NC,IDS,X,Y,Move)
 IF( Move==-1 )Then
  Print*,'Can not Move the Mesh'
  pause
 EndIF

!Part 6:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,Vol)
!************************************************************************************
 End
!####################################################################################