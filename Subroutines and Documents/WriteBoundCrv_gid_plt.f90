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
 Subroutine WriteBoundCrv_gid_plt(Dim,NP,NR,NFR,IDS,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,NR,NFR,IDS,X,Y

 Integer::Dim,I,J,SumBedg,NP,NR
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
!**********************************************************************************************	

 Open(1,File='BoundCrv.plt')

 SumBedg=0
 Do I=1,NR
   
    IF( IDS(2,SumBedg+1)==0 )Then
    
     Write(1,*) ' TITLE = "Title" ' 
     Write(1,*) ' VARIABLES  = X , Y '
     Write(1,*) ' ZONE T=" Title',  '", N= ', NP , ' , E= ', NFR(I) ,', ET=LINESEG, F=FEBLOCK'

     Do J=1,NP
        Write(1,*) X(J)
     End Do
     Do J=1,NP
        Write(1,*) Y(J)
     End Do
     
     Do J=SumBedg+1,SumBedg+NFR(I)
  	    Write(1,*) IDS(3,J),IDS(4,J)
     End Do

    EndIf
    
	SumBedg = SumBedg + NFR(I)
 End Do
 
 !close(1)
!*********************************************************************************************
 End
!###########################################################################################

