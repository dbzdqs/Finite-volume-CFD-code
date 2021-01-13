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
 Subroutine First_Elem_Orient(Dim,X,Y,Z,NC,Corn,Neib,FirstOE)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,X,Y,Z,NC
 Intent (InOut)::Neib,Corn
 Intent (Out  )::FirstOE

 Integer::Dim,I,N1,N2,N3,NC,Dumy,Nearst,FirstOE
 Real(8)::Xout,Yout,Zout,Xc,Yc,Zc,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X21,Y21,X31,Y31,Ax,Ay,Az,Min,&
          Dis,Volum
 Real(8),Dimension(1:Dim)::X,Y,Z
 Integer,Dimension(1:4,1:Dim)::Corn,Neib
!*********************************************************************************************
!Part 1:
 If( Dabs( Maxval(Z)-Minval(Z) )<0.0000001 )Then

 !Part 2:
  X21 = X(2) - X(1)
  Y21 = Y(2) - Y(1)

  X31 = X(3) - X(1)
  Y31 = Y(3) - Y(1)

  Volum = (X21*Y31) - (Y21*X31)

 !Part 3:
  If( Volum<0. )Then
   Dumy      = Corn(2,1)
   Corn(2,1) = Corn(3,1)
   Corn(3,1) = Dumy

   Dumy       = Neib(2,1)
   Neib(2,1)  = Neib(3,1)
   Neib(3,1)  = Dumy
  Endif

 !Part 4:
  FirstOE =  1
  Goto 100

 Endif


!Part 5:
 Xout = 2*Maxval(X) 
 Yout = 2*Maxval(Y) 
 Zout = 2*Maxval(Z) 

!Part 6:
 Min=100000000.0
 Do I=1,NC
 
    N1 = Corn(1,I)
	N2 = Corn(2,I)
	N3 = Corn(3,I)

    Xc = ( X(N1)+X(N2)+X(N3) ) / 3
	Yc = ( Y(N1)+Y(N2)+Y(N3) ) / 3
	Zc = ( Z(N1)+Z(N2)+Z(N3) ) / 3

    Dis = (Xc-Xout)**2 + (Yc-Yout)**2 + (Zc-Zout)**2
	
	If(Dis<Min)Then
	 Min = Dis
	 Nearst = I
	Endif

 End Do

!Part 7:
 N1 = Corn(1,Nearst)
 N2 = Corn(2,Nearst)
 N3 = Corn(3,Nearst)

!Part 8:
 X1= X(N1)-Xout  ;  Y1= Y(N1)-Yout  ;  Z1= Z(N1)-Zout 
 X2= X(N2)-Xout  ;  Y2= Y(N2)-Yout  ;  Z2= Z(N2)-Zout 
 X3= X(N3)-Xout  ;  Y3= Y(N3)-Yout  ;  Z3= Z(N3)-Zout 

 Ax =  Y1*Z3 - Z1*Y3
 Ay = -X1*Z3 + Z1*X3
 Az =  X1*Y3 - Y1*X3

 Volum = X2*Ax + Y2*Ay + Z2*Az 

!Part 9:
 FirstOE =  Nearst

!Part 10:
 If( Volum>0.0 )Then
  Dumy            = Corn(2,FirstOE)
  Corn(2,FirstOE) = Corn(3,FirstOE)
  Corn(3,FirstOE) = Dumy

  Dumy            = Neib(2,FirstOE)
  Neib(2,FirstOE) = Neib(3,FirstOE)
  Neib(3,FirstOE) = Dumy
 Endif
!*********************************************************************************************
 100 End 
!*********************************************************************************************
