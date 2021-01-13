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
 Subroutine Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)
 Implicit None
!*********************************************************************************************
 Intent(In)::Dim,NC,NP,Corn,X,Y,GM,WNP1,P

 Integer::Dim,I,J,P1,P2,P3,P4,ME,NP,NC,NFW1,NFW2,NF
 Real(8)::Minf,GM,CP,CCP,Xm,U,V,CC,Mach
 Integer,Dimension(1:4,1:Dim)::Corn,IDS
 Real(8),Dimension(1:Dim)::X,Y,P,Xc,Yc
 Real(8),Dimension(1:4,1:Dim)::WNP1
!*********************************************************************************************	
!Part 1:
 Open(1,File='Contours.Plt')


!Part 3:
 Write(1,*) 'Variables="X","Y","Ro","U","V","P","Mach"'
 Write(1,*) 'ZONE N=' ,   NP , ' E=' ,  NC  
 Write(1,*) ' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK VARLOCATION=([3-7]=CELLCENTERED)'

 Do J=1,NP
	Write(1,*) X(J)
 End Do
 Do J=1,NP
	Write(1,*) Y(J) 
 End Do

 Do J=1,NC
	Write(1,*) WNP1(1,J)
 End Do
 Do J=1,NC
	Write(1,*) WNP1(2,J)/WNP1(1,J)
 End Do
 Do J=1,NC
	Write(1,*) WNP1(3,J)/WNP1(1,J)
 End Do
 Do J=1,NC
	Write(1,*) P(J)
 End Do   
 Do J=1,NC
    U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    CC = GM*P(J)/WNP1(1,J)
    Mach = SQRT((U*U+V*V)/CC)
	Write(1,*) Mach
 End Do


 Do I=1,NC
    P1 = Corn(1,I) 
    P2 = Corn(2,I) 
    P3 = Corn(3,I)
    P4 = Corn(4,I)
    if(P4==0) P4=P3

	Write(1,*) P1,P2,P3,P4
 End Do

!Part 10:
 Close(1)
!**********************************************************************************************
 End
!##############################################################################################

