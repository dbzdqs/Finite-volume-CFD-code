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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Read_2DMeshC_TriToQuad(Dim,NP,NC,NBC,NEC,BEP,Corn,Neib,X,Y)
Implicit None
!===========================================================================================
Intent(In   )::Dim
Intent(Out  )::NP,NC,NBC,NEC,BEP,Corn,Neib,X,Y

Integer::Dim,I,J,J1,NP,NBC,NC,NBE,Dumy
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Real(8),Dimension(1:Dim)::X,Y
Integer,Dimension(1:10)::NEC
Integer,Dimension(1:Dim,1:2)::BEP
!===========================================================================================
!Part 1: 
Open(1,File='MeshIn.cgid')
Open(2,File='MeshIn.Plt')

!Part 2: 
Read(1,*) NP 

!Part 3: 
Read(1,*) NC  

!Part 4:
Read(1,*) NBC  

!Part 5:
Do I=1,NBC
	Read(1,*) NEC(I) 
End Do

!Part 6: 
NBE=0
Do J=1,NBC
	Do J1=NBE+1,NBE+NEC(J)
	   Read(1,*) BEP(J1,1) , BEP(J1,2)
	End Do
	NBE = NBE + NEC(J)
End Do

!Part 7: 
Do I=1,NC
	Read(1,*) Corn(I,1),Corn(I,2),Corn(I,3),Corn(I,4)
	Read(1,*) Neib(I,1),Neib(I,2),Neib(I,3),Neib(I,4)  
End do

!Part 8:
Do I=1,NP
Read(1,*) X(I),Y(I)  
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
!===========================================================================================
End Subroutine Read_2DMeshC_TriToQuad
!*********************************************************************************************
