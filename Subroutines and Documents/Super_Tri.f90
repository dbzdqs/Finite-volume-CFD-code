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
 Subroutine Super_Tri(Dim,NP,X,Y,NC,Corn,Neib)
 Implicit none
!*********************************************************************************************
 Intent(In   )::Dim
 Intent(Out  )::NC,Corn,Neib
 Intent(Inout)::NP,X,Y

 Integer::Dim,NC,NP
 Real(8)::Xmid,Ymid,R,Xmax,Xmin,Ymax,Ymin
 Integer,Dimension(1:Dim,1:4)::Corn,Neib
 Real(8),Dimension(1:Dim)::X,Y
!*********************************************************************************************
!part1 :
 Xmax=Maxval(X)
 Xmin=Minval(X)

 Ymax=Maxval(Y)
 Ymin=Minval(Y)

!part2 : 
 Xmid = 0.5*(Xmax+Xmin)
 Ymid = 0.5*(Ymax+Ymin)
 R    = 0.5*Max( abs(Xmax-Xmin) , abs(Ymax-Ymin) )

!part3 : 
 X(NP+1)=Xmid 
 Y(NP+1)=Ymid + 2*R

 X(NP+2)=Xmid - 4*R
 Y(NP+2)=Ymid - 2*R

 X(NP+3)=Xmid + 4*R
 Y(NP+3)=Ymid - 2*R

!part4 : 
 Corn(1,1)=NP+1
 Corn(1,2)=NP+2
 Corn(1,3)=NP+3

 Neib(1,1)=0
 Neib(1,2)=0
 Neib(1,3)=0

!part5 : 
 NC=1
 NP=NP+3
!*********************************************************************************************
 End 
!###########################################################################################
