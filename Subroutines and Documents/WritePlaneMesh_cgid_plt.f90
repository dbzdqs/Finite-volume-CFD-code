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
!// Date: Nov., 15, 2014                                                                   //!
!// Developed by: *//*-+/                       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine WritePlaneMesh_cgid_plt(Dim,NP,NC,Corn,X,Y,Z)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,NC,Corn,X,Y,Z

 Integer::Dim,I,NP,NC,Dumy
 Integer,Dimension(1:4,1:Dim)::Corn
 Real(8),Dimension(1:Dim)::X,Y,Z
!*********************************************************************************************
 !part 1:
 Open(1,File='PlaneMesh.plt')

 Write(1,*) 'Variables="X","Y","Z"'
 Write(1,*) 'Zone T="Grid"'
 Write(1,*) ' N=  ', NP, ',E= ' , NC, ',F=FEPOINT ET=QUADRILATERAL'
 Do I=1,NP
	Write(1,*) X(I),Y(I),Z(I)
 End Do
 
 !part 2:
 Do I=1,NC
    Dumy=Corn(4,I)
    if(Dumy==0) Dumy=Corn(3,I)
	Write(1,*) Corn(1,I),Corn(2,I),Corn(3,I),Dumy
 End Do

 Close(1)
!*********************************************************************************************
 End
!###########################################################################################