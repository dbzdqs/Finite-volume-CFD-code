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
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine VolumeCheck3D(Dim,NF,NC,NP,IDS,X,Y,Z,FaceType,NConectCell_Loc,IConectCell_Loc,NFace_Cell,IFace_Cell,NegativeVol)
 Implicit None
!*********************************************************************************************
 Intent(In   )::NF,NP,NC,IDS,X,Y,Z,FaceType,NConectCell_Loc,IConectCell_Loc,NFace_Cell,IFace_Cell
 Intent(Out  )::NegativeVol

 Integer::Dim,I,J,K,NF,P1,P2,P3,ME,NE,NFacePnt,IFace,NC,NP,NConectCell_Loc,Cell,Face,NConectFaces
 Real(8)::a,b,c,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,X21,X32,Y21,Y32,Z21,Z32,DV,Vol
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:Dim)::X,Y,Z
 Integer,Dimension(1:100)::IConectCell_Loc
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::IConectFaces
 Logical::NegativeVol
!*********************************************************************************************	 
!Part 1:
 NegativeVol = .FALSE.

!Part 2:
 Do K=1,NConectCell_Loc
    Cell = IConectCell_Loc(K)
    Vol  = 0.0
   
   !Part 3:
    Do J=1,NFace_Cell(Cell)
       Face     = IFace_Cell(J,Cell)
       NFacePnt = FaceType(Face) + 2
       DV       = 0.0 
    
      !Part 4:
       Do I=4, NFacePnt - 1
          P1 = IDS(3  ,Face)
	      P2 = IDS(I  ,Face)
	      P3 = IDS(I+1,Face)

	      X1 = X(P1) ; Y1 = Y(P1) ; Z1 = Z(P1)
	      X2 = X(P2) ; Y2 = Y(P2) ; Z2 = Z(P2)
	      X3 = X(P3) ; Y3 = Y(P3) ; Z3 = Z(P3)

          x21 = x2 - x1
          x32 = x3 - x2

          y21 = y2 - y1
          y32 = y3 - y2

          z21 = z2 - z1
          z32 = z3 - z2

          a = y21*z32 - z21 * y32
          b = z21*x32 - x21 * z32
          c = x21*y32 - y21 * x32

	      DV = DV + a*( x1+ x2+ x3 ) + b*( y1+ y2+ y3 ) + c*( z1+ z2+ z3 )
       End Do
        
       DV = DV / 18
        
      !Part 5:
       If(IDS(1,Face)==Cell)Then
        Vol =  Vol + DV
       Else
        Vol =  Vol - DV
       Endif
        
    EndDo
    
   !Part 6:
    !If(Vol < 0.0 .AND. abs(Vol)>0.00000000001)Then
    If(Vol < 0.0 .AND. abs(Vol)>0.00000000001)Then
     NegativeVol = .TRUE.
     exit
    Endif
    
 EndDo
!*********************************************************************************************
 End
!###########################################################################################