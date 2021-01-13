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
 Subroutine Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,FileName)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,NR,NFR,IDS,X,Y

 Integer::Dim,J,NP,NF,NR,I,sumNF 
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
 Integer,Dimension(1:100)::NFR
Character *9 FileName
!*********************************************************************************************
!Part 1:
Open(16,File=FileName)

 sumNF=0
 Do I=1,NR

!Part 2:
 Write(16,*) ' TITLE = "Title" ' 
 Write(16,*) ' VARIABLES  = X , Y '
 Write(16,*) ' ZONE T=" Title',  '", N= ', NP , ' , E= ', NFR(I) ,', ET=LINESEG, F=FEBLOCK'

!Part 3:
 Do J=1,NP
	Write(16,*) X(J)
 End Do
 Do J=1,NP
	Write(16,*) Y(J)
 End Do

!Part 1:
 Do J=sumNF+1,sumNF+NFR(I)        
    Write(16,*) IDS(3,J),IDS(4,J)
 End do

sumNF=sumNF+NFR(I) 
 End do
 
 close(16)
!*********************************************************************************************
 End
!###########################################################################################

