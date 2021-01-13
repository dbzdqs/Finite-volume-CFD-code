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
!// Date: June, 10, 2017                                                                   //!
!// Developed by: R. Amery, Mathmatical, Amirkabir university of Technology                //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine GenBLayerLastPt2D(Dim,NPtCurvs,NBL,BLThick,EdgInvSlop,NPBL,BLPt,XBL,YBL)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NPtCurvs,NBL,BLThick,EdgInvSlop
 Intent(InOut)::NPBL,BLPt,XBL,YBL
 
 Integer::Dim,I,J,P1,P2,Pt
 Real(8)::X1,Y1,Xn,Yn,SlopX,SlopY,Thicknes,Len1
 Integer::NPtCurvs
 Integer::NBL
 Integer::NPBL
 Integer,Dimension(1:25,1:Dim)::BLPt
 Real(8),Dimension(1:Dim)::BLThick
 Real(8),Dimension(1:2,1:Dim)::EdgInvSlop
 Real(8),Dimension(1:Dim)::XBL,YBL
!*********************************************************************************************
!Part 1:
 Do I=1,NPtCurvs 	
    
   !Part 2:
    X1 = XBL(I)
	Y1 = YBL(I)
    
   !Part 3:
    SlopX = EdgInvSlop(1,I)
    SlopY = EdgInvSlop(2,I)
    
   !Part 4:
    Xn = X1 + SlopX*BLThick(I)
	Yn = Y1 + SlopY*BLThick(I)
    
   !Part 5:
    NPBL = NPBL + 1
	XBL(NPBL) = Xn
	YBL(NPBL) = Yn
    
   !Part 6:
    BLPt(NBL+1,I) = NPBL

 End do
!*********************************************************************************************
 END
!###########################################################################################