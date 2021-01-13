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
 Subroutine Linear_Spring2D(Dim,NConectPoints,IConectPoints,NP,X,Y,DelX,DelY)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NConectPoints,IConectPoints,NP,X,Y
 Intent(InOut)::DelX,DelY

 Integer::Dim,I,J,J1,JJ,NP
 Real(8)::Res,SumX_Old,SumX_New,SumX,SumY_Old,SumY_New,SumY,Stiff,Sum_Stiff
 Real(8),Dimension(1:Dim)::X,Y,DelX,DelY
 Integer,Dimension(1:10,1:Dim)::IConectPoints
 Integer,Dimension(1:Dim)::NConectPoints
!*********************************************************************************************

!Part 1:
 Res=100000.0
 Do While( Dabs(Res)>0.00000001 )

   !Part 2:
    SumX_Old = 0.0
    SumY_Old = 0.0

    SumX_New = 0.0
    SumY_New = 0.0

   !Part 3:
    Do J=1,NP
	
      !Part 4:
       If(NConectPoints(J)==0)Cycle
       
      !Part 5:
	   SumX_Old = SumX_Old + DelX(J)
	   SumY_Old = SumY_Old + DelY(J)
    
      !Part 6:
       SumX      = 0.0
       SumY      = 0.0
	   Sum_Stiff = 0.0
    
      !Part 7:
	   Do J1=1,NConectPoints(J)
		  JJ=IConectPoints(J1,J)

		  Stiff=1/Dsqrt( (X(Jj)-X(J))**2 + (Y(Jj)-Y(J))**2 )

          SumX = SumX + Stiff*DelX(JJ)
		  SumY = SumY + Stiff*DelY(JJ)

		  Sum_Stiff = Sum_Stiff + Stiff
	   End Do
	    
      !Part 8:
       SumX = SumX/Sum_Stiff
       SumY = SumY/Sum_Stiff

      !Part 9:
       DelX(J) = DelX(J) + 1.4 * ( SumX - DelX(J) )
       DelY(J) = DelY(J) + 1.4 * ( SumY - DelY(J) )

      !Part 10:
	   SumX_New = SumX_New + DelX(J)
	   SumY_New = SumY_New + DelY(J)

    End Do

   !Part 11:
    Res = Abs( ((SumX_New+SumY_New) - (SumX_Old+SumY_Old)) / (NP*(SumX_Old+SumY_Old)) )
  
 End Do
!*********************************************************************************************
 End 
!###########################################################################################
