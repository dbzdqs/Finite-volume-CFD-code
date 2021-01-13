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
!// Chief Developer: N. msnkre,Aerospace eng. Amirkabir University of Technology           //!
!// Supervisor: Dr. h. hdhrnuidn,Aerospace eng. Amirkabir University of Technology         //!
!// Date: May.,04,2018                                                                     //!
!// Developed by: N. msnkre,Aerospace Eng.,Amirkabir University of Technology              //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied,Modified and Redistributed for Non-Commercial Use.                    //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
 Subroutine Write_CFwithwallFunction(Dim,NFW1,NFW2,IDS,Minf,Rinf,X,TauWall)
 Implicit None
!**********************************************************************************************
INTEGER                      ,INTENT(IN) :: DIM
INTEGER                      ,INTENT(IN) :: NFW1
INTEGER                      ,INTENT(IN) :: NFW2
INTEGER,DIMENSION(1:4,1:DIM) ,INTENT(IN) :: IDS
REAL(8)                      ,INTENT(IN) :: MINF
REAL(8)                      ,INTENT(IN) :: RINF
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN) :: X
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN) :: TAUWALL

INTEGER                                  ::I,P1,P2,ME
REAL(8)                                  ::xm,CFC
!**********************************************************************************************	
!Part 1:
 Open(21,File='CF.Plt')
 
!Part 2:
 CFC = 2/(Rinf*Minf)

!Part 3:
 DO I=NFW1+1,NFW2

   !Part 4:
    ME = IDS(1,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
   
   !Part 5:
    Xm = 0.5*( X(P1)+X(P2) )

   !Part 6:
    Write(21,*) Xm, CFC * TauWall(I)

 End do
 
 Close(21)
!**********************************************************************************************
 End
!##############################################################################################
 
 
 
 
 
