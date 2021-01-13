!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Calculate the Convection Terms of 2D Mean Flow Equations by Centeral Scheme!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F009F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine LUSGS(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace,Vol,DT,WNP1&
                 &,WB,P,Mu,Mut,PrL,PrT,X,Y,Z,xc,yc,zc,MR,DW)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace,Vol,DT,WNP1,&
                &WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,X,Y,Z,MR
 Intent (Out  )::DW

 Integer::Dim,I,NC,NF1,NF2,NF,K,KK,KKK,J,Neib,Face,NeqI
 Real(8)::DD,GM,PrL,PrT,MR,IsCN
 Real(8),Dimension(1:5,1:Dim)::Res,WNP1,DW_star,DW
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,eigen,DT,Vol,P,Mu,Mut,xc,yc,zc,X,Y,Z
 Integer,Dimension(1:10,1:Dim)::InonzeroCell,InonzeroFace
 Integer,Dimension(1:Dim)::NnonzeroCell
 Real(8),Dimension(1:5)::Jacobi
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:5,1:5)::Jaco
!*******************************************************************************************
!Part 1: !Calculate eigen value of inviscid and viscous part using Tamador scheme (if "Mu" is zero it works for inviscid)
 Call Calculate_eigValMeanFlow(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,X,Y,Z,xc,yc,zc,MR,eigen)

!Part 2: !If IsCN=2.0 => Second-order Crank-Nicelson scheme. If IsCN=1.0 => First-order Euler scheme.
 IsCN=2.0

!Part 3: !If NeqI=5   => Capable of 3D.                      If NeqI=4   => Capable of 2D.
 NeqI=5

!Part 4: !forward sweep
 Do J=1,NC

   !Part 5:
    Do k=1,NeqI
       DW_star(k,j)=0.0
       DW     (k,j)=0.0
    End do

   !Part 6:
    DO k=1,NnonzeroCell(j)

      !Part 7:
       Neib = InonzeroCell(k,j)
       Face = InonzeroFace(k,j)

      !Part 8:
       IF( Neib<J .and. Neib/=0 )Then

!!!!    !=====================================================================================================
       !Part 9.1: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using Increment scheme
        Call Increment(GM,NX(Face),NY(Face),NZ(Face) ,WNP1(1,Neib),WNP1(2,Neib),&
                       WNP1(3,Neib),WNP1(4,Neib),WNP1(5,Neib),DW_star(1,Neib),&
                       DW_star(2,Neib),DW_star(3,Neib),DW_star(4,Neib),DW_star(5,Neib),Jacobi)

!!!!    !-----------------------------------------------------------------------------------------------------
!!       !Part 9.2: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using analytical Jacobian
!        Call Calculate_Jacobi_neighb(Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Face,Neib,WNP1,WB,Jaco)
!
!        Do kk=1,NeqI
!           Jacobi(kk)=0.0
!           Do kkk=1,NeqI
!              Jacobi(kk) = Jacobi(kk) + Jaco(kk,kkk)*DW_star(kkk,Neib)*DA(Face)
!           End do
!        End do
!!!!    !=====================================================================================================
       !Part 10:
        Do kk=1,NeqI
           DW_star(kk,j) = DW_star(kk,j) + 0.5*(Jacobi(kk) - DA(Face)*eigen(Face)*DW_star(kk,Neib))
        End Do
       End if
    End Do

   !Part 11:
    DD = IsCN*Vol(j)/DT(j)
    DO k=1,NnonzeroCell(j)
       Face = InonzeroFace(k,j)
       DD = DD + 0.5*eigen(Face)*DA(Face)
    End do

   !Part 12:
    Do kk=1,NeqI
       DW_star(kk,j) = ( DW_star(kk,j)-IsCN*Res(kk,J) )/DD
    End do
 End do ! End Sweep Forward


!Part 13: backward sweep
 Do J=NC,1,-1

   !Part 14:
    Do k=1,NeqI
       DW(k,j)=0.0
    End do

   !Part 15:
    DO k=1,NnonzeroCell(j)

      !Part 16:
       Neib = InonzeroCell(k,j)
       Face = InonzeroFace(k,j)

      !Part 17:
       IF( Neib<J .and.  Neib/=0 )then

!!!!   !=====================================================================================================
       !Part 18.1: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using Increment scheme
       Call Increment(GM,NX(Face),NY(Face),NZ(Face) ,WNP1(1,Neib),WNP1(2,Neib),&
                       WNP1(3,Neib), WNP1(4,Neib), WNP1(5,Neib),DW(1,Neib),DW(2,Neib),&
                       DW(3,Neib),DW(4,Neib),DW(5,Neib),Jacobi)
!!!!   !-----------------------------------------------------------------------------------------------------
!       !Part 18.2: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using analytical Jacobian
!        Call Calculate_Jacobi_neighb(Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Face,Neib,WNP1,WB,Jaco)
!
!        Do kk=1,NeqI
!           Jacobi(kk)=0.0
!           Do kkk=1,NeqI
!              Jacobi(kk) = Jacobi(kk) + Jaco(kk,kkk)*DW(kkk,Neib)*DA(Face)
!           End do
!        End do
!!!!   !=====================================================================================================
       !Part 19:
        Do kk=1,NeqI
           DW(kk,j) = DW(kk,j) + 0.5*( Jacobi(kk) - DA(Face)*eigen(Face)*DW(kk,Neib) )
        End do
       Endif
    End do

   !Part 20:
    DD=IsCN*Vol(j)/DT(j)
    DO k=1,NnonzeroCell(j)
       Face = InonzeroFace(k,j)
       DD = DD + 0.5*eigen(Face)*DA(Face)
    End do

   !Part 21:
    Do kk=1,NeqI
       DW(kk,j) = DW_star(kk,j) + DW(kk,j)/DD
    End do
 End do !End Backward sweep

!*******************************************************************************************
 End
!###########################################################################################
