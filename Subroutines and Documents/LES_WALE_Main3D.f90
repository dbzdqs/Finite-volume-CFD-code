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
!// Developed by: E. Akrami, Mechanical Eng., Amirkabir University of Technology           //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine LES_WALE_Main(Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,DW,Mut)
 Implicit None
!********************************************************************************************* 
 Intent(In   )::Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,DW
 Intent(InOut)::Mut

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE
 Real(8)::MR,Rc,Uc,Vc,Wc,R1,U1,V1,W1,& 
          Sxxhat,Sxyhat,Sxzhat,Syxhat,Syyhat,Syzhat,Szxhat,Szyhat,Szzhat,& 
          Wxxhat,Wxyhat,Wxzhat,Wyxhat,Wyyhat,Wyzhat,Wzxhat,Wzyhat,Wzzhat,& 
          Sdxxhat,Sdxyhat,Sdxzhat,Sdyxhat,Sdyyhat,Sdyzhat,Sdzxhat,Sdzyhat,Sdzzhat,& 
          SijSijhat,WijWijhat,SdijSdij,Ls
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DW,Vol,Mut
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::DUX,DUY,DUZ ,DVX,DVY,DVZ ,DWX,DWY,DWZ
 Real(8),Dimension(1:Dim)::DUXhat,DUYhat,DUZhat,DVXhat,DVYhat,DVZhat,DWXhat,DWYhat,DWZhat
 Real(8),Dimension(1:Dim)::U,V,W,UhatF,VhatF,WhatF,Uhat,Vhat,What
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    Rc = WNP1(1,I)
    Uc = WNP1(2,I)/Rc
    Vc = WNP1(3,I)/Rc
    Wc = WNP1(4,I)/Rc

    U(I)   = Uc
    V(I)   = Vc
    W(I)   = Wc

 End Do

!Part 2:
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,U,Vol,Uhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,V,Vol,Vhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,W,Vol,What)

!Part 3:
 Do I=NF1+1,NF2

    ME = IDS(1,I)
    NE = IDS(2,I)

    UhatF(I) = 0.5*(Uhat(ME)+Uhat(NE))
    VhatF(I) = 0.5*(Vhat(ME)+Vhat(NE))
    WhatF(I) = 0.5*(What(ME)+What(NE))

 End Do

 Do I=NF2+1,NF

    R1 = WB(1,I)
    U1 = WB(2,I)/R1
    V1 = WB(3,I)/R1
    W1 = WB(4,I)/R1

    UhatF(I) = U1
    VhatF(I) = V1
    WhatF(I) = W1

 End Do

!Part 4:
 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,UhatF,DUXhat,DUYhat,DUZhat)
 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,VhatF,DVXhat,DVYhat,DVZhat)
 Call LES_Cell_Gradiant3D(Dim,NC,NF1,NF2,NF,Vol,NX,NY,NZ,IDS,VhatF,DWXhat,DWYhat,DWZhat)
 
 Do I=1,NC
 
   !Part 5: 
    Sxxhat = 0.5 * (DUXhat(I)+DUXhat(I)) 
    Sxyhat = 0.5 * (DUYhat(I)+DVXhat(I))
    Sxzhat = 0.5 * (DUZhat(I)+DWXhat(I))
    
    Syxhat = 0.5 * (DVXhat(I)+DUYhat(I))
    Syyhat = 0.5 * (DVYhat(I)+DVYhat(I))
    Syzhat = 0.5 * (DVZhat(I)+DWYhat(I))
    
    Szxhat = 0.5 * (DWXhat(I)+DUZhat(I))
    Szyhat = 0.5 * (DWYhat(I)+DVZhat(I))
    Szzhat = 0.5 * (DWZhat(I)+DWZhat(I))

   !Part 6:
    Wxxhat = 0.5 * (DUXhat(I)-DUXhat(I))
    Wxyhat = 0.5 * (DUYhat(I)-DVXhat(I))
    Wxzhat = 0.5 * (DUZhat(I)-DWXhat(I))
    
    Wyxhat = 0.5 * (DVXhat(I)-DUYhat(I))
    Wyyhat = 0.5 * (DVYhat(I)-DVYhat(I))
    Wyzhat = 0.5 * (DVZhat(I)-DWYhat(I))
    
    Wzxhat = 0.5 * (DWXhat(I)-DUZhat(I))
    Wzyhat = 0.5 * (DWYhat(I)-DVZhat(I))
    Wzzhat = 0.5 * (DWZhat(I)-DWZhat(I))

   !Part 7:
    SijSijhat = Sxxhat*Sxxhat + Sxyhat*Sxyhat + Sxzhat*Sxzhat + &
                Syxhat*Syxhat + Syyhat*Syyhat + Syzhat*Syzhat + &
                Szxhat*Szxhat + Szyhat*Szyhat + Szzhat*Szzhat

    WijWijhat = Wxxhat*Wxxhat + Wxyhat*Wxyhat + Wxzhat*Wxzhat + &
                Wyxhat*Wyxhat + Wyyhat*Wyyhat + Wyzhat*Wyzhat + &
                Wzxhat*Wzxhat + Wzyhat*Wzyhat + Wzzhat*Wzzhat

   !Part 8:
    Sdxxhat = Sxxhat*Sxxhat+Sxyhat*Syxhat+Sxzhat*Szxhat + Wxxhat*Wxxhat+Wxyhat*Wyxhat+Wxzhat*Wzxhat - (SijSijhat-WijWijhat)/3.0
    Sdxyhat = Sxxhat*Sxyhat+Sxyhat*Syyhat+Sxzhat*Szyhat + Wxxhat*Wxyhat+Wxyhat*Wyyhat+Wxzhat*Wzyhat
    Sdxzhat = Sxxhat*Sxzhat+Sxyhat*Syzhat+Sxzhat*Szzhat + Wxxhat*Wxzhat+Wxyhat*Wyzhat+Wxzhat*Wzzhat
    
    Sdyxhat = Syxhat*Sxxhat+Syyhat*Syxhat+Syzhat*Szxhat + Wyxhat*Wxxhat+Wyyhat*Wyxhat+Wyzhat*Wzxhat
    Sdyyhat = Syxhat*Sxyhat+Syyhat*Syyhat+Syzhat*Szyhat + Wyxhat*Wxyhat+Wyyhat*Wyyhat+Wyzhat*Wzyhat - (SijSijhat-WijWijhat)/3.0
    Sdyzhat = Syxhat*Sxzhat+Syyhat*Syzhat+Syzhat*Szzhat + Wyxhat*Wxzhat+Wyyhat*Wyzhat+Wyzhat*Wzzhat
    
    Sdzxhat = Szxhat*Sxxhat+Szyhat*Syxhat+Szzhat*Szxhat + Wzxhat*Wxxhat+Wzyhat*Wyxhat+Wzzhat*Wzxhat
    Sdzyhat = Szxhat*Sxyhat+Szyhat*Syyhat+Szzhat*Szyhat + Wzxhat*Wxyhat+Wzyhat*Wyyhat+Wzzhat*Wzyhat
    Sdzzhat = Szxhat*Sxzhat+Szyhat*Syzhat+Szzhat*Szzhat + Wzxhat*Wxzhat+Wzyhat*Wyzhat+Wzzhat*Wzzhat - (SijSijhat-WijWijhat)/3.0

   !Part 9:
    SdijSdij = Sdxxhat*Sdxxhat + Sdxyhat*Sdxyhat + Sdxzhat*Sdxzhat + &
               Sdyxhat*Sdyxhat + Sdyyhat*Sdyyhat + Sdyzhat*Sdyzhat + &
               Sdzxhat*Sdzxhat + Sdzyhat*Sdzyhat + Sdzzhat*Sdzzhat

   !Part 10:
    Ls = min(0.4*DW(I) , Vol(I)**(1./3.)*0.325) 

   !Part 11:
    Mut(I)= 1.0/MR*WNP1(1,I)*Ls**2*(SdijSdij**(3.0/2.0)) / (  SijSijhat**(5.0/2.0)  +  SdijSdij**(5.0/4.0)  )

 End Do

!********************************************************************************************* 
    End 
!########################################################################################### 

