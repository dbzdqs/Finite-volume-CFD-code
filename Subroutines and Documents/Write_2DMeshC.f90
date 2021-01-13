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
Subroutine Write_2DMeshC(Dim,NP,NC,NBC,NEC,BEP,Corn,Neib,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,NC,NBC,NEC,BEP,Corn,Neib,X,Y

 Integer::Dim,I,J,J1,NP,NBC,NC,NBE,Dumy
 Integer,Dimension(1:Dim,1:4)::Corn,Neib,BEP
 Real(8),Dimension(1:Dim)::X,Y
 Integer,Dimension(1:10)::NEC
!*********************************************************************************************
!Part 1:
 Open(1,File='2DMeshC.Txt')
 Open(2,File='2DMeshC.Plt')

!Part 2: 
 Write(1,*) NP  ,' Number of Points'
 
!Part 3: 
 Write(1,*) NC  ,' Number of Cells'

!Part 4:
 Write(1,*) NBC ,' Number of Boundary Curves' 

!Part 5:
 Do I=1,NBC
    Write(1,*) NEC(I) , ' Number of Edges Belong to Each Curves' 
 End Do

!Part 6: 
 NBE=0
 Do J=1,NBC
    Do J1=NBE+1,NBE+NEC(J)
       Write(1,*) BEP(J1,1) , BEP(J1,2)
    End Do
	NBE = NBE + NEC(J)
 End Do

!Part 7: 
 Do I=1,NC
	Write(1,'(3x,I5,2x,I5,2x,I5,2x,I5, 5x,I5)') Corn(I,1),Corn(I,2),Corn(I,3),Corn(I,4) , I
	Write(1,'(3x,I5,2x,I5,2x,I5,2x,I5, 5x,I5)') Neib(I,1),Neib(I,2),Neib(I,3),Neib(I,4) , I   
 End do
 
!Part 8:
 Do I=1,NP
	Write(1,*) X(I),Y(I)  
 End do


!Part 9:
if(NC/=0)then
 Write(2,*) 'Variables="X","Y"'
 Write(2,*) 'Zone T="Grid"'
 Write(2,*) ' N=  ', NP, ',E= ' , NC, ',F=FEPOINT ET=QUADRILATERAL'
 Do I=1,NP
	Write(2,*) X(I),Y(I)
 End Do
 Do I=1,NC
    Dumy=Corn(I,4)
    if(Dumy==0) Dumy=Corn(I,3)
	Write(2,*) Corn(I,1),Corn(I,2),Corn(I,3),Dumy
 End Do
Endif

!Part 10:
 NBE=0
 Do I=1,NBC
    
    Write(2,*) 'Variables="X","Y"'
    Write(2,*) 'Zone N=',NP,' E=',NEC(I),',Datapacking=Point, Zonetype=Fetriangle'
    Do J=1,NP
	   Write(2,*) X(J),Y(J)
    End Do
    Do J=NBE+1,NBE+NEC(I)
  	   Write(2,*) BEP(J,1),BEP(J,2),BEP(J,2)
    End Do

	NBE = NBE + NEC(I)
 End Do

 Close(1)
 Close(2)
!*********************************************************************************************
 End
!*********************************************************************************************
