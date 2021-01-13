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
 Subroutine Tri_Contain_Point(Dim,NC,Corn,X,Y,Xn,Yn,Itri)
 Implicit None            
!*********************************************************************************************
 Intent(In   )::Dim,NC,Corn,X,Y,Xn,Yn
 Intent(Out  )::Itri

 Integer::Dim,J,Jp1,NC,I,E,Itri
 Real(8)::Xn,Yn,Dxn,Dyn,Dot,Dsx,Dsy
 Integer,Dimension(1:Dim,1:4)::Corn
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1:3)::Xt,Yt
!*********************************************************************************************
!part 1:
 Do I=1,NC

   !part 2: 
    Do J=1,3
	   Xt(J)=X(Corn(I,J))
	   Yt(J)=Y(Corn(I,J))
	End Do

   !part 3: 
    E=0
	Do J=1,3

	   Jp1=J+1
	   If(J==3)Jp1=1
      
	  !part 4:
	   Dsy =- ( Xt(Jp1) - Xt(J) )
	   Dsx =  ( Yt(Jp1) - Yt(J) )

      !part 5: 
	   Dxn = Xn - (Xt(Jp1)+Xt(J))/2
	   Dyn = Yn - (Yt(Jp1)+Yt(J))/2
      
	  !part 6: 
	   Dot = Dsx*Dxn + Dsy*Dyn
      
	  !part 7: 
	   If(Dot<=0.0) E=E+1

    End Do
  
   !part 8: 
	If(E==3)then
	Itri=I
	Exit
	else
    Itri=-1
	End If

 End Do
!*********************************************************************************************
 End 
!*********************************************************************************************
