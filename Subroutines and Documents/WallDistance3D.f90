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
 Subroutine WallDistance3D(Dim,NC,NFW1,NFW2,FaceType,IDS,X,Y,Z,Xc,Yc,Zc,INW,DW)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,FaceType,IDS,X,Y,Z,Xc,Yc,Zc
 Intent(Out  )::INW,DW

 Integer::Dim,J,I,II,JJ,Pi,NC,NFW1,NFW2,FacTyp
 Real(8)::Dmin,Dis,Xj,Yj,Zj,Xi,Yi,Zi,DX,DY,DZ
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW,FaceType
 Real(8),Dimension(1:Dim)::X,Y,Z,Xc,Yc,Zc,DW
!*********************************************************************************************	
!Part 1:
 Do J=1,NC

   !Part 2:
	Xj = Xc(J)
	Yj = Yc(J)
	Zj = Zc(J)
   
   !Part 3:    
    Dmin=1000000.0

   !Part 4:
    Do I=NFW1+1,NFW2
      
	  !Part 5:
       Xi = 0.0
       Yi = 0.0
       Zi = 0.0
       FacTyp = FaceType(I)
       Do JJ=1,FacTyp
          Pi = IDS(JJ+2,I)
          
          Xi = Xi + X(Pi)
          Yi = Yi + Y(Pi)
          Zi = Zi + Z(Pi)
       End Do
       Xi = Xi/FacTyp
       Yi = Yi/FacTyp
       Zi = Zi/FacTyp
       
       
	  !Part 6:  
	   DX = Xj-Xi
	   DY = Yj-Yi
	   DZ = Zj-Zi
       Dis = Dsqrt(DX*DX+DY*DY+DZ*DZ) 
      
	  !Part 7:
	   If(Dis<Dmin)then
	    Dmin = Dis
	    II   = I
	   Endif

	End do
      
   !Part 8: 
    DW(j)  = Dmin
    INW(j) = II

 End do
!*********************************************************************************************
 End
!###########################################################################################