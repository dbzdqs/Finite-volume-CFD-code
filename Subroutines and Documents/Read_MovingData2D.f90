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
!// Chief Developer: M. Namvar, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. A. Jahangirian, Aerospace eng. Amirkabir University of Technology      //!
!// Date: Feb., 10, 2018                                                                   //!
!// Developed by: M. Namvar, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Read_MovingData2D(Dim,NR,Nstp,Vel_X,Vel_Y,AnulVelz,Xo1,Yo1)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NR
 Intent(Out  )::Nstp,Vel_X,Vel_Y,AnulVelz,Xo1,Yo1

 Integer::Dim,i,j,K,NStp,NBP
 Integer,Dimension(1:Dim)::IBP
 Real(8)::Dt,Pi,xx , yy, ww
 Integer::NR
 Real(8),Dimension(1:100,1:100)::DX,DY,DZ,LinVel,Dtetx,Dtety,Dtetz,AnulVelx,AnulVely, AnulVelz,Vel_X,Vel_Y,Vel_W
 Real(8),Dimension(1:100):: Xo1,Yo1,Zo1,Xo2,Yo2,Zo2,MoveTime
!*********************************************************************************************
!Part 1:
 Open(8,File='MovigData.Txt')

 Pi=4.*atan(1.)
 
 Read(8,*)
 Read(8,*) NStp

 Do K=1,NR
    Read(8,*) Xo1(k),Yo1(k),Zo1(k) , Xo2(k),Yo2(k),Zo2(k)
 End do

 Read(8,*)

 Do I=1,NStp
    Do K=1,NR
       Read(8,*) DX(i,k) , DY(i,k) , DZ(i,k) , LinVel(i,k) , Dtetx(i,k) , Dtety(i,k) , Dtetz(i,k) , AnulVelx(i,k) , AnulVely(i,k) , AnulVelz(i,k)
    End do
 End do

 Do I=1,NStp
    Do K=1,NR
       Vel_X(i,k) = LinVel(i,k) * ( DX(i,k)/sqrt(DX(i,k)*DX(i,k) + DY(i,k)*DY(i,k)+0.00000001 ) ) 
       Vel_Y(i,k) = LinVel(i,k) * ( DY(i,k)/sqrt(DX(i,k)*DX(i,k) + DY(i,k)*DY(i,k)+0.00000001 ) )
    End do
 End do

 Close(8)
!*********************************************************************************************
 End
!#############################################################################################
