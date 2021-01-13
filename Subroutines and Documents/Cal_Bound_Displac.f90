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
!// Date: May., 15, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Cal_Bound_Displac(Dim,Vel_X,Vel_Y,Vel_W,Xo,Yo,IDS,IStp,NR,NFR,BC,DT,X,Y,DelX,DelY)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Vel_X,Vel_Y,Vel_W,IDS,IStp,NR,NFR,BC,DT,X,Y
 Intent(Out  )::DelX,DelY
 Intent(Inout)::XO,YO

 Integer::Dim,IStp,NR,I,J,K,Sum_Fac,ii
 Real(8)::DT,Omega,U,V,Move_RotatX,Move_RotatY,Move_LinearX,Move_LinearY
 Real(8),Dimension(1:Dim)::X,Y,DelX,DelY
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:100)::NFR,BC
 Real(8),Dimension(1:100,1:100)::Vel_X,Vel_Y,Vel_W
 Real(8),Dimension(1:100)::Xo,Yo
!*********************************************************************************************

 Sum_Fac=0
 Do k=1,NR

    Omega = Vel_W(IStp,K) 
    U     = Vel_X(IStp,K)  
    V     = Vel_Y(IStp,K) 

	if( Dabs(Omega)<0.00000001 .and. Dabs(U)<0.00000001 .and. Dabs(V)<0.00000001 ) goto 10 

    Do J=Sum_Fac+1,Sum_Fac+NFR(k)

       I = IDS(3,J)
	   ii=ii+1
      
       Move_LinearX = U*dt 
       Move_LinearY = V*dt

       Move_RotatX = ((x(i)-Xo(k))*(dcos(Omega*dt)-1.0)+(y(i)-Yo(k))*dsin(Omega*dt))
	   Move_RotatY = ((y(i)-Yo(k))*(dcos(Omega*dt)-1.0)-(x(i)-Xo(k))*dsin(Omega*dt))

       DelX(i) = Move_RotatX + Move_LinearX
       DelY(i) = Move_RotatY + Move_LinearY

	end do

    xo(k)=xo(k)+Move_LinearX
    yo(k)=yo(k)+Move_LinearY

10  Sum_Fac = Sum_Fac + NFR(k)
 END DO
!*********************************************************************************************
 End
!###########################################################################################
