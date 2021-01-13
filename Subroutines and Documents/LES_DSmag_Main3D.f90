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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine LES_DSmag_Main(Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,Mut,Taukk)
 Implicit None
!********************************************************************************************* 
Intent(In   )::Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ
Intent(InOut)::Taukk,Mut

Integer::Dim,I,NC,NF1,NF2,NF,Allocatestatus,DeAllocatestatus,ME,NE
Real(8)::MR,R1,U1,V1,W1,R2,U2,V2,W2,Rc,Uc,Vc,Wc

Real(8),Dimension(1:Dim)::NX,NY,NZ,Vol,Taukk,Mut
Integer,Dimension(1:6,1:Dim)::IDS
Real(8),Dimension(1:5,1:Dim)::WNP1
Real(8),Dimension(1:6,1:Dim)::WB
Real(8),Dimension(1:Dim) ::U,V,W,&
                           DUX,DUY,DUZ ,DVX,DVY,DVZ ,DWX,DWY,DWZ,&
                           Uhat,Vhat,What,&
                           DUXhat,DUYhat,DUZhat ,DVXhat,DVYhat,DVZhat ,DWXhat,DWYhat,DWZhat,&
                           Rho    ,RhoU    ,RhoV    ,RhoW    ,RhoUU    ,RhoVV    ,RhoWW    ,RhoUV    ,RhoUW    ,RhoWV,&
		                   Rhohat ,RhoUhat ,RhoVhat ,RhoWhat ,RhoUUhat ,RhoVVhat ,RhoWWhat ,RhoUVhat ,RhoUWhat ,RhoWVhat,&
                           Sxx    ,Syy     ,Szz     ,SXY     ,SXZ     ,SZY     ,Skk      ,Sabs     ,&
		                   Sxxhat ,Syyhat  ,Szzhat  ,Sxyhat  ,Sxzhat  ,Szyhat  ,Skkhat   ,Sabshat  ,&
                           Lxx    ,Lyy     ,Lzz     ,Lxy     ,Lxz     ,Lzy     ,Lkk      ,&
		                   VF     ,UF      ,WF      ,VhatF   ,UhatF   ,WhatF  
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    Rc = WNP1(1,I)
    Uc = WNP1(2,I)/Rc
    Vc = WNP1(3,I)/Rc
    Wc = WNP1(4,I)/Rc

    Rho(I) = Rc
    U(I)   = Uc
    V(I)   = Vc
    W(I)   = Wc

    RhoU(I) = WNP1(2,I)
    RhoV(I) = WNP1(3,I)
    RhoW(I) = WNP1(4,I)

    RhoUU(I) = Rc*Uc*Uc
    RhoVV(I) = Rc*Vc*Vc
    RhoWW(I) = Rc*Wc*Wc
    
    RhoUV(I) = Rc*Uc*Vc
    RhoUW(I) = Rc*Uc*Wc
    RhoWV(I) = Rc*Wc*Vc
 End Do

!Part 2:
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Rho,Vol,Rhohat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,U,Vol,Uhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,V,Vol,Vhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,W,Vol,What)

 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoU,Vol,RhoUhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoV,Vol,RhoVhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoW,Vol,RhoWhat)

 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoUU,Vol,RhoUUhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoVV,Vol,RhoVVhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoWW,Vol,RhoWWhat)
 
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoUV,Vol,RhoUVhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoUW,Vol,RhoUWhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,RhoWV,Vol,RhoWVhat)

!Part 3:
 Do I=NF1+1,NF2

    ME = IDS(1,I)
    NE = IDS(2,I)

    UF(I) = 0.5*(U(ME)+U(NE))
    VF(I) = 0.5*(V(ME)+V(NE))
    WF(I) = 0.5*(W(ME)+W(NE))
   
    UhatF(I) = 0.5*(Uhat(ME)+Uhat(NE))
    VhatF(I) = 0.5*(Vhat(ME)+Vhat(NE))
    WhatF(I) = 0.5*(What(ME)+What(NE))

 End Do

 Do I=NF2+1,NF

    R1 = WB(1,I)
    U1 = WB(2,I)/R1
    V1 = WB(3,I)/R1
    W1 = WB(4,I)/R1

    UF(I) = U1
    VF(I) = V1
    WF(I) = W1

    UhatF(I) = U1
    VhatF(I) = V1
    WhatF(I) = W1

 End Do


!Part 4:
 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,UF,DUX,DUY,DUZ)
 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,VF,DVX,DVY,DVZ)
 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,VF,DWX,DWY,DWZ)

 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,UhatF,DUXhat,DUYhat,DUZhat)
 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,VhatF,DVXhat,DVYhat,DVZhat)
 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,WhatF,DWXhat,DWYhat,DWZhat)
    
!Part 5:  
 Do I=1,NC
       
   !Part 6: 
    Sxx(I) = 0.5 * (DUX(I)+DUX(I))
    Syy(I) = 0.5 * (DVY(I)+DVY(I))
    Szz(I) = 0.5 * (DWZ(I)+DWZ(I))
    
    Sxy(I) = 0.5 * (DUY(I)+DVX(I))
    Sxz(I) = 0.5 * (DUZ(I)+DWX(I))
    Szy(I) = 0.5 * (DWY(I)+DVZ(I))

   !Part 7:
    Skk(I)  = Sxx(I) + Syy(I) + SZZ(I)
    Sabs(I) = Dsqrt( 2*( Sxx(I)*Sxx(I) + Syy(I)*Syy(I) + Szz(I)*Szz(I) + 2*Sxy(I)*Sxy(I) + 2*Sxz(I)*Sxz(I) + 2*Szy(I)*Szy(I) ) )
        
   !Part 8:
    Sxxhat(I) = 0.5 * (DUXhat(I)+DUXhat(I))
    Syyhat(I) = 0.5 * (DVYhat(I)+DVYhat(I))
    Szzhat(I) = 0.5 * (DWZhat(I)+DWZhat(I))
    
    Sxyhat(I) = 0.5 * (DUYhat(I)+DVXhat(I))
    Sxzhat(I) = 0.5 * (DUZhat(I)+DWXhat(I))
    Szyhat(I) = 0.5 * (DWYhat(I)+DVZhat(I))
       
   !Part 9:
    Skkhat(I) = Sxxhat(I) + Syyhat(I) + Szzhat(I)
    Sabshat(I) = Dsqrt( 2*(  Sxxhat(I)*Sxxhat(I) +  Syyhat(I)*Syyhat(I) +   Szzhat(I)*Szzhat(I) + &
                           2*Sxyhat(I)*Sxyhat(I) +2*Sxzhat(I)*Sxzhat(I) + 2*Szyhat(I)*Szyhat(I)  )   )
       
   !Part 10: 
    Lxx(I) = RhoUUhat(I) - RhoUhat(I)*RhoUhat(I)/Rhohat(I)
    Lyy(I) = RhoVVhat(I) - RhoVhat(I)*RhoVhat(I)/Rhohat(I)
    Lzz(I) = RhoWWhat(I) - RhoWhat(I)*RhoWhat(I)/Rhohat(I)
    
    Lxy(I) = RhoUVhat(I) - RhoUhat(I)*RhoVhat(I)/Rhohat(I)
    Lxz(I) = RhoUWhat(I) - RhoUhat(I)*RhoWhat(I)/Rhohat(I)
    Lzy(I) = RhoWVhat(I) - RhoWhat(I)*RhoVhat(I)/Rhohat(I)
    
    Lkk(I) = Lxx(I)+Lyy(I)+Lzz(I)
    
 End Do

!Part 11:
 Call LES_DSmag_Eddy3D(Dim,NC,NF1,NF2,IDS,Vol,MR,Rho,Rhohat,Skk,Skkhat,Sabs,Sabshat,Lkk,Mut,&
                              Sxx,Syy,Szz,Sxy,Sxz,Szy,&
                              Sxxhat,Syyhat,Szzhat,Sxyhat,Sxzhat,Szyhat,&
                              Lxx,Lyy,Lzz,Lxy,Lxz,Lzy)
!Part 12:
 Call LES_DSmag_IsoSGS3D(Dim,NC,NF1,NF2,IDS,MR,Vol,Rho,Rhohat,Sabs,Sabshat,Lkk,Taukk)
!********************************************************************************************* 
 End 
!########################################################################################### 