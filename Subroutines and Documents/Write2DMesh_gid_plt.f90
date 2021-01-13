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
Subroutine Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,FileName)
Implicit None
!*********************************************************************************************
Intent(In   )::Dim,NP,NF,IDS,X,Y
Character *9 FileName
Integer::Dim,J,NP,NF
Integer,Dimension(1:4,1:Dim)::IDS
Real(8),Dimension(1:Dim)::X,Y
!*********************************************************************************************
!Part 1:
Open(189,File=FileName)

!Part 2:
Write(189,*) ' TITLE = "Title" '
Write(189,*) ' VARIABLES  = X , Y '
Write(189,*) ' ZONE T=" Title',  '", N= ', NP , ' , E= ', NF ,', ET=LINESEG, F=FEBLOCK'

!Part 3:
Do J=1,NP
    Write(189,*) X(J)
End Do
Do J=1,NP
    Write(189,*) Y(J)
End Do

!Part 1:
Do J=1,NF
    Write(189,*) IDS(3,J),IDS(4,J)
End do
!*********************************************************************************************
 END 
!###########################################################################################