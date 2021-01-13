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
 Subroutine Write_CF3DUnsteady(Dim,NFW1,NFW2,IDS,FaceType,X,Y,Z,Minf,Rinf,Mu,DUY,Naverage,SumCF)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,IDS,FaceType,X,Y,Z,Minf,Rinf,Mu,DUY,Naverage
 Intent(InOut)::SumCF

 Integer::Dim,I,J,ME,NFW1,NFW2,Naverage,Pi
 Real(8)::Minf,Rinf,CFC,CF,Xm,Ym,Zm
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:Dim)::X,Y,Z,Mu,DUY,SumCF
!*********************************************************************************************	
!Part 1:
 Open(102,File='UnstsyCF.Plt')
 Open(4,File='AveragedCF.Plt')

 WRITE(102,*)'VARIABLES="X","UnstsyCF"'  
 WRITE(102,*)'ZONE'
 
 WRITE(4,*)'VARIABLES="X","AveragedCF"'  
 WRITE(4,*)'ZONE'
 
!Part 2: 
 CFC = 2/(Rinf*Minf)
 
 !Part 3:
 DO I=NFW1+1,NFW2


!Part 4:
    Zm = 0.0
    Xm = 0.0
    Do J=1,FaceType(I)
       Pi = IDS(J+2,I)
       Zm = Zm + Z(Pi)
       Xm = Xm + X(Pi)
    End Do
    Zm = Zm / FaceType(I) 
    Xm = Xm / FaceType(I)

!Part 5:    
    IF( Zm>0.0 .and. Zm<0.6 )Then
        
     ME = IDS(1,I)
     
     CF = CFC * Mu(ME)*DUY(I)
     Write(102,*) Xm,  CF

!Part 6:
     SumCF(I) = SumCF(I) + CF
    
     CF = SumCF(I)/Naverage
     
 !Part 7:    
	 Write(4,*) Xm,CF
    End IF
    
 End do
 
 
!Part 10:
 Close(4)
!**********************************************************************************************
 End
!##############################################################################################