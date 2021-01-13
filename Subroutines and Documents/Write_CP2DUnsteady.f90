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
 Subroutine Write_CP2DUnsteady(Dim,NFW1,NFW2,IDS,X,Y,Minf,P,GM,Naverage,SumCP)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,IDS,X,Y,Minf,P,GM,Naverage
 Intent(InOut)::SumCP

 Integer::Dim,I,J,ME,NFW1,NFW2,Naverage,P1,P2
 Real(8)::Minf,GM,CCP,CP,Xm,Ym
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,P,SumCP
!*********************************************************************************************	
!Part 1:
 Open(101,File='UnstsyCP.Plt')
 Open(4,File='AveragedCP.Plt')

 WRITE(101,*)'VARIABLES="X","UnstsyCP"'  
 WRITE(101,*)'ZONE'
 
 WRITE(4,*)'VARIABLES="X","AveragedCP"'  
 WRITE(4,*)'ZONE'
 
!Part 2: 
 CCP = 0.5*GM*Minf*Minf
 
 !Part 3:
 DO I=NFW1+1,NFW2

    !Part 4:
    Xm = 0.5*( X(P1)+X(P2) )
    
    !Part 5:
     ME = IDS(1,I)
     CP = (GM*P(ME)-1.0)/CCP
     Write(101,*) Xm,  CP

     !Part 6:
     SumCP(I) = SumCP(I) + CP
     CP = SumCP(I)/Naverage
     
     !Part 7:
	 Write(4,*) Xm,CP
    
 End do
 
 
!Part 10:
 Close(4)
!**********************************************************************************************
    End
!##############################################################################################
  