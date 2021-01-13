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
 Subroutine Write3DMeshSepRgn_gid_plt(Dim,NP,NF,NR,NFR,IDS,X,Y,Z,FaceType,FileName)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,NF,NR,NFR,IDS,X,Y,Z,FaceType

 Integer::Dim,I,J,JJ,NP,NC,NF,NR,SFace
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:Dim)::X,Y,Z
Character *9 FileName
!*********************************************************************************************

 Open(2,File=FileName)

 SFace=0
 Do I=1,NR

    Write(2,*) ' TITLE = "Title" '            
    Write(2,*) ' VARIABLES  = X , Y , Z'
    Write(2,*) ' ZONE T="Title", N= ', NP , ' , E= ', NFR(I) ,', ET=QUADRILATERAL, F=FEBLOCK'

    Do J=1,NP
	   Write(2,*) X(J)
    End Do
    Do J=1,NP
	   Write(2,*) Y(J)
    End Do
    Do J=1,NP
	   Write(2,*) Z(J)
    End Do

    Do J=SFace+1,SFace+NFR(I)
	   Do JJ=3,FaceType(J)+2            
          Write(2,'(I15)',Advance='No') IDS(JJ,J)         
       End Do
       
	   IF(FaceType(J)==3) Write(2,'(I15)',Advance='No') IDS(FaceType(J)+2,J )
       Write(2,*)
    End Do
    SFace=SFace+NFR(I)

 End do
 
 Close(2) 
!*********************************************************************************************
 End
!###########################################################################################

