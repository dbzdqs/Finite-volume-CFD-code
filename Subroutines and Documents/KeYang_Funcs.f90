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
!// Developed by: M. H. Saadat, Aerospace Eng., Amirkabir University of Technology         //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeYang_Funcs(Dim,NC,IDS,X,Y,Xc,Yc,NX,NY,A,DW,INW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,&
                         DDVY,DDKX,DDKY,DDEPSX,DDEPSY,DD2UY,Mu,Mut,f1,f2,fmu,Lk,Le)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,IDS,X,Y,Xc,Yc,NX,NY,A,DW,INW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,DDKX,&
                DDKY,DDEPSX,DDEPSY,DD2UY,Mu,Mut
 Intent(Out  )::f1,f2,fmu,Lk,Le
 
 Integer::Dim,I,II,NC,ME,P1,P2
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut,X,Y,Xc,Yc,f1,f2,fmu,Lk,Le
 Real(8),Dimension(1:2,1:Dim)::Wntp1
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:Dim)::DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDEPSX,DDEPSY,DD2UY,NX,NY
 Real(8)::k,Epsilon,Rho,Yn,DX,DY,DX2,DY2,DL,Area,MR,U,DUX,DUY,tauwall,ustar,yplus,RY,Ret,DDUN
!***************************************************************************************************  
!Part 1:
 Do I=1,NC
        
    Rho     = Wnp1(1,I)
    k       = Wntp1(1,I)/Rho
    Epsilon = Wntp1(2,I)/Rho
        
    II = INW(I)  
    Yn = DW(I)  
    ME = IDS(1,II)
    P1 = IDS(3,II)
    P2 = IDS(4,II)
    DX = X(P2)-X(P1)
    DY = Y(P2)-Y(P1)
    DX2 = Xc(ME)-X(P1)
	DY2 = Yc(ME)-Y(P1)
    DL=Dsqrt(DX*DX + DY*DY)
    AREA = 0.5*Dabs( DX*DY2 - DY*DX2 )
    U = Wnp1(2,ME)/Wnp1(1,ME)
        
    DUX = -0.5*U*DY/AREA
    DUY =  0.5*U*DX/AREA
    DDUN=(DUX*DY-DUY*DX)/DL

    tauwall = (MR)*Mu(ME)*DDUN 
        
    ustar = Dsqrt(abs(tauwall/Wnp1(1,ME)))
    yplus = (1.0/MR)*Rho*ustar*Yn/Mu(I)
        
    Ret =  (Rho*k*k)/(MR*Mu(I)*Epsilon)
    Ry  =  Rho*Yn*sqrt(abs(k))/(MR*Mu(I))
    f1(I)  =  1.0 / (1.0 + 1.0/dsqrt(Ret) )
    f2(I)  =  1.0 / (1.0 + 1.0/dsqrt(Ret) )
    fmu(I) =  (1.0 + 1.0/sqrt(Ret) )*dsqrt(dabs(1.0 - exp( (-1.5d-4)*Ry - (5.0d-7)*Ry*Ry*Ry - (1.0d-10)*Ry*Ry*Ry*Ry*Ry )) )
        
    Lk(I)  =  0.0 
    Le(I)  =  MR*MR*Mu(I)*Mut(I)*DD2UY(I)*DD2UY(I)
      
 END Do
!*********************************************************************************************
 End
!###########################################################################################

