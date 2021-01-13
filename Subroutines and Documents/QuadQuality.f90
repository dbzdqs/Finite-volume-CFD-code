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
!// Date: Apr., 25, 2013                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function QuadQuality(Dim,X,Y,A,B,C,D) !------>> Implementation Based on Lee (1994) <<-------
Implicit None
!=========================================================================================== 
Intent(In)::Dim,X,Y,A,B,C,D

Integer::Dim,A,B,C,D,I,J
Real(8)::TriangleQuality,QuadQuality,key
Real(8),Dimension(1:4)::Values
Real(8),Dimension(1:Dim)::X,Y
!=========================================================================================== 
!Part 1:
!-------------------->>> Notice: Suppose ABCD is the Quad then: <<<-------------------------

Values(1) = TriangleQuality(Dim,X,Y,A,B,C)
Values(2) = TriangleQuality(Dim,X,Y,A,C,D)
Values(3) = TriangleQuality(Dim,X,Y,B,C,D)
Values(4) = TriangleQuality(Dim,X,Y,A,B,D)

!------------ Sorting alpha values using 'INSERSTION SORT' (descending order) --------------
do I=2,4

	key = Values(I)
	J = I-1

	do
        
        if(J < 1) then
            exit    
        elseif(Values(J)<key) then
            exit    
        endif

		Values(J+1) = Values(J)
		J = J-1

	end do

	Values(J+1) = key

end do

QuadQuality = (Values(3)*Values(4))/(Values(1)*Values(2))

!=========================================================================================== 
End Function QuadQuality
!*********************************************************************************************
