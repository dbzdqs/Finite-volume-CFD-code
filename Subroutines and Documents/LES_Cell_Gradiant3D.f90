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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,PHI,DPHIX,DPHIY,DPHIZ)
 Implicit None
!********************************************************************************************* 
 Intent(IN)::Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,PHI
 Intent(InOut)::DPHIX,DPHIY,DPHIZ

 Integer::Dim,NC,NF1,NF2,NF,ME,NE,I
 Real(8)::PHI_NX,PHI_NY,PHI_NZ
 Real(8),Dimension(1:Dim)::Vol,NX,NY,NZ,PHI,DPHIX,DPHIY,DPHIZ
 Integer,Dimension(1:6,1:Dim)::IDS
!*********************************************************************************************
 DPHIX(:)=0.0
 DPHIY(:)=0.0 
 DPHIZ(:)=0.0
    
!Part 1:
 DO I=NF1+1,NF2

   !Part 2:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 3:
    PHI_NX = PHI(I)*NX(I)
    PHI_NY = PHI(I)*NY(I)
    PHI_NZ = PHI(I)*NZ(I)

   !Part 4:
    DPHIX(ME)=DPHIX(ME) + PHI_NX
    DPHIY(ME)=DPHIY(ME) + PHI_NY  
    DPHIZ(ME)=DPHIZ(ME) + PHI_NZ         

    DPHIX(NE)=DPHIX(NE) - PHI_NX
    DPHIY(NE)=DPHIY(NE) - PHI_NY 
    DPHIZ(NE)=DPHIZ(NE) - PHI_NZ 
 
 End Do

!Part 5:
 DO I=NF2+1,NF

    ME = IDS(1,I)

    PHI_NX = PHI(I)*NX(I)
    PHI_NY = PHI(I)*NY(I)
    PHI_NZ = PHI(I)*NZ(I)

    DPHIX(ME)=DPHIX(ME) + PHI_NX
    DPHIY(ME)=DPHIY(ME) + PHI_NY 
    DPHIZ(ME)=DPHIZ(ME) + PHI_NZ         
 
 End Do

!Part 6:
 DO I=1,NC    
    DPHIX(I)=DPHIX(I)/Vol(I)
    DPHIY(I)=DPHIY(I)/Vol(I)  
    DPHIZ(I)=DPHIZ(I)/Vol(I)
 End Do
!*********************************************************************************************   
 End
!###########################################################################################