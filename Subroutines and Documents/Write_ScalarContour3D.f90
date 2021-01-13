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
 Subroutine Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,Scalar,FileName)
 Implicit None
!**********************************************************************************************
 Intent(In)::Dim,NC,NP,X,Y,Z,Corn,Scalar

 Integer::Dim,I,J,NP,NC
 Integer,Dimension(1:8,1:Dim)::Corn
 INTEGER::P1,P2,P3,P4,P5,P6,P7,P8
 Real(8),Dimension(1:Dim)::X,Y,Z,Scalar
 Character *10 FileName
 Real(8)::GM
 Real(8),Dimension(1:4,1:Dim)::WNP1
!**********************************************************************************************	
!Part 1:
 Open(105,File=FileName)

!Part 2:
 Write(105,*) 'Variables="X","Y","Z","Scalar","no" '
 Write(105,*) 'ZONE N=' ,   NP , ' E=' ,  NC  
 Write(105,*) ' ZONETYPE=FEBRICK DATAPACKING=BLOCK VARLOCATION=([4-5]=CELLCENTERED)'


!Part 3:
 Do J=1,NP
	Write(105,*) X(J)
 End Do
 Do J=1,NP
	Write(105,*) Y(J) 
 End Do
 Do J=1,NP
	Write(105,*) Z(J) 
 End Do

!Part 4:
 Do J=1,NC
	Write(105,*) Scalar(J)
 End Do
 Do J=1,NC
	Write(105,*) '0.0'
 End Do
 
 !Part 5:
 DO I=1,NC
    P1 = Corn(1,I) 
    P2 = Corn(2,I) 
    P3 = Corn(3,I)
    P4 = Corn(4,I)
    P5 = Corn(5,I) 
    P6 = Corn(6,I) 
    P7 = Corn(7,I)
    P8 = Corn(8,I)

	Write(105,*) P1,P2,P3,P4,P5,P6,P7,P8
 End Do
 

 Close(105)
!**********************************************************************************************	
 END 
!##############################################################################################
 
 
 
 
