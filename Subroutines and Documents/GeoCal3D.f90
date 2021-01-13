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
!// Date: April, 01, 2017                                                                  //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine GeoCal3D(Dim,NF,NC,IDS,X,Y,Z,FaceType,Vol,DA,Nx,Ny,Nz,XC,YC,ZC)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF,NC,IDS,X,Y,Z,FaceType
 Intent(Out  )::Vol,DA,Nx,Ny,Nz,XC,YC,ZC

 Integer::Dim,I,J,K,NF,P1,P2,P3,ME,NE,NFacePnt,IFace,NC
 Real(8)::a,b,c,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,X21,X32,Y21,Y32,Z21,Z32,DV,SumX,SumY,SumZ,aX,bX,cX,aY,bY,cY,aZ,bZ,cZ,vol_tt,VOLUME,Xcc,Ycc,Zcc
 Integer,Dimension(1:6,1:Dim)::IDS
  Integer,Dimension(1:4)::Pf
 Integer,Dimension(1:Dim)::FaceType,NPC1
 Real(8),Dimension(1:Dim)::X,Y,Z,Vol,Nx,Ny,Nz,DA,XC,YC,ZC
!*********************************************************************************************	
!Part 1:
 DO IFace=1,NF

    DV =0.0 
    Xcc=0.0
    Ycc=0.0
    Zcc=0.0
    
    NFacePnt = FaceType(IFace)+2
    Do I=4,NFacePnt-1

       P1= IDS(3  ,IFace)
	   P2= IDS(I  ,IFace)
	   P3= IDS(I+1,IFace)

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

       DA(IFace) = DA(IFace) + 0.5*Dsqrt( a*a + b*b + c*c )
       
       Xcc=Xcc+a*((X1+X2)**2+(X2+X3)**2+(X1+X3)**2) 
       Ycc=Ycc+b*((Y1+Y2)**2+(Y2+Y3)**2+(Y1+Y3)**2) 
       Zcc=Zcc+c*((Z1+Z2)**2+(Z2+Z3)**2+(Z1+Z3)**2) 
   

    End Do 
    DV = DV /18
     
    ME = IDS(1,IFace)
    NE = IDS(2,IFace)
 
    NX(IFace) = a / Dsqrt( a*a + b*b + c*c )
	NY(IFace) = b / Dsqrt( a*a + b*b + c*c )
	NZ(IFace) = c / Dsqrt( a*a + b*b + c*c )

    NX(IFace) = NX(IFace) * DA(IFace)
	NY(IFace) = NY(IFace) * DA(IFace)
	NZ(IFace) = NZ(IFace) * DA(IFace)

    XC(ME)= XC(ME)+Xcc
    YC(ME)= YC(ME)+Ycc
    ZC(ME)= ZC(ME)+Zcc
    IF(NE/=0) THEN
    XC(NE)= XC(NE)-Xcc
    YC(NE)= YC(NE)-Ycc
    ZC(NE)= ZC(NE)-Zcc    
    END IF
    
    Vol(ME) =  Vol(ME) + DV
	IF(NE/=0)Vol(NE) =  Vol(NE) - DV

 End do
 
 DO I=1,NC
    XC(I)= XC(I)/(48*vol(I))
    YC(I)= YC(I)/(48*vol(I))
    ZC(I)= ZC(I)/(48*vol(I)) 
 End do
!*********************************************************************************************
 End
!###########################################################################################


