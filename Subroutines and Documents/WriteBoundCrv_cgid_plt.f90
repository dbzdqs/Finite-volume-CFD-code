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
!// Date: Oct., 05, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine WriteBoundCrv_cgid_plt(Dim,NP,NBoundCrv,NFacCrv,BFacPt,X,Y,Z)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,NBoundCrv,NFacCrv,BFacPt,X,Y,Z

 Integer::Dim,I,J,SumBedg,NP,NBoundCrv
 Integer,Dimension(1:100)::NFacCrv
 Integer,Dimension(1:Dim,1:2)::BFacPt
 Real(8),Dimension(1:Dim)::X,Y,Z
!**********************************************************************************************	

 Open(1,File='BoundCrv.plt')

 SumBedg=0
 Do I=1,NBoundCrv
    
    Write(1,*) ' TITLE = "Title" ' 
    Write(1,*) ' VARIABLES  = X , Y , Z '
    Write(1,*) ' ZONE T=" Title',  '", N= ', NP , ' , E= ', NFacCrv(I) ,', ET=LINESEG, F=FEBLOCK'

    Do J=1,NP
       Write(1,*) X(J)
    End Do
    Do J=1,NP
       Write(1,*) Y(J)
    End Do
    Do J=1,NP
       Write(1,*) Z(J)
    End Do

    Do J=SumBedg+1,SumBedg+NFacCrv(I)
  	   Write(1,*) BFacPt(J,1),BFacPt(J,2)
    End Do

	SumBedg = SumBedg + NFacCrv(I)
 End Do
 
 close(1)
!*********************************************************************************************
 End
!###########################################################################################

