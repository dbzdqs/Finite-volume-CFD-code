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
 Subroutine RegionsPointData(Dim,NR,NFR,IDS,PtRegion,PtIndx)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NR,NFR,IDS
 Intent(Out  )::PtRegion,PtIndx
 
 Integer::Dim,I,J,JJ,Sum,P,Cnt,FacTyp
 Integer::NR
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:100,1:2)::PtRegion
 Integer,Dimension(1:Dim)::IP
 Integer,Dimension(1:Dim)::PtIndx
!*********************************************************************************************
 FacTyp=2

 Cnt=0
 Sum=0
 Do I=1,NR
  
    IP(:)=0

    PtRegion(I,1)=Cnt+1
    DO J=Sum+1,Sum+NFR(I)

       Do JJ=3,FacTyp+2
	      P = IDS(JJ,J)	
    
          If( IP(P)==0 )Then
           IP(P)=1
         
           Cnt=Cnt+1
           PtIndx(Cnt) = P
          End If
       End Do

    END DO
    PtRegion(I,2)=Cnt

    Sum = Sum + NFR(I)		
 END DO

!*********************************************************************************************
 END
!###########################################################################################
