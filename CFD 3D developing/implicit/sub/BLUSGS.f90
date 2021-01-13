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
 Subroutine BLUSGS(Dim,IDS,NC,NF1,NF2,NF,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace&
                            &,Vol,P,Mu,Mut,PrL,PrT,X,Y,Z,xc,yc,zc,DT,MR,WNP1,WB,DW)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,IDS,NC,NF1,NF2,NF,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace,Vol&
                &,P,Mu,Mut,PrL,PrT,X,Y,Z,xc,yc,zc,DT,MR,WNP1,WB
 Intent (Out  )::DW

 Integer::Dim,I,NC,NF1,NF2,NF,K,KK,kkk,J,Neib,Face,NeqI,cell,iter,MaxIter
 Real(8)::DD,GM,PrL,PrT,Er,Cri,MR,IsCN
 Real(8),Dimension(1:5,1:Dim)::Res,WNP1,DW_star,DW,DQ
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,eigen,DT,Vol,Mu,Mut,P,xc,yc,zc,X,Y,Z
 Integer,Dimension(1:10,1:Dim)::InonzeroCell,InonzeroFace
 Integer,Dimension(1:Dim)::NnonzeroCell
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:5)::Jacob_Neib,Var,Var_star,Varj
 Real(8),Dimension(1:5,1:5)::Jacobv,Inv,eigenMatrix,Jaco
 LOGICAL::OK_FLAG
!*******************************************************************************************

!Part 1: !If NeqI=5   => Capable of 3D.                      If NeqI=4   => Capable of 2D.
 NeqI=5

!Part 2: !If IsCN=2.0 => Second-order Crank-Nicelson scheme. If IsCN=1.0 => First-order Euler scheme.
 IsCN=2.0

!Part 3.1:
 Er=10.0
 iter=0

!Part 3.2:
 Cri=0.1    !Error criteria for convergence of forward and backward sweeps
 MaxIter=1  !Max. allowable Iterations for convergence of forward and backward sweeps

!Part 4:
 DQ(:,:)=DW(:,:)
 DW_star(:,:)=0.0
 DW     (:,:)=0.0

!Part 5:
 Do While(Er>Cri)

   !Part 6:
    iter=iter+1
    if(iter>MaxIter)then
        exit
    end if
    Er=0.0

    !Part 7: forward sweep
     Do J=1,NC

       !Part 8:
        Do k=1,NeqI
           Var_star(k)=0.0
           Var     (k)=0.0
        End do

       !Part 9:
        DO k=1,NnonzeroCell(j)

          !Part 10:
           Neib = InonzeroCell(k,j)
           Face = InonzeroFace(k,j)

          !Part 11: !lower triangular
           IF( Neib<J .and. Neib/=0 )Then

           !Part 12:
            Call Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                    &,MR,k,Face,Jacobv)

           !Part 13:
            Call Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face,eigenMatrix)

!!!!       !=====================================================================================================
!          !Part 14.1: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using Increment scheme
!           Call Increment(GM,NX(Face),NY(Face),NZ(Face) ,WNP1(1,Neib),WNP1(2,Neib),&
!                               WNP1(3,Neib),WNP1(4,Neib),WNP1(5,Neib),DW_star(1,Neib),&
!                               DW_star(2,Neib),DW_star(3,Neib),DW_star(4,Neib),DW_star(5,Neib),Jacob_Neib)
!!!!       !-----------------------------------------------------------------------------------------------------
!          !Part 14.2: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using analytical Jacobian
            Call Calculate_Jacobi_neighb(Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Face,Neib,WNP1,WB,Jaco)

            Do kk=1,NeqI
               Jacob_Neib(kk)=0.0
               Do kkk=1,NeqI
                   Jacob_Neib(kk) = Jacob_Neib(kk) + Jaco(kk,kkk)*DW_star(kkk,Neib)*DA(Face)
                End do
            End do
!!!!       !=====================================================================================================

           !Part 15:
            Do kk=1,NeqI
               Varj(kk)=0.0
               Do kkk=1,NeqI
                  Varj(kk) = Varj(kk) + (0.5*(Jacobv(kk,kkk)+eigenMatrix(kk,kkk))*DA(Face))*DW_star(kkk,Neib)
               End Do
            End Do

           !Part 16:
            Do kk=1,NeqI
               Var_star(kk) = Var_star(kk) + 0.5*Jacob_Neib(kk) -Varj(kk)
            End Do
           End if !End lower triangular

          !Part 17: !upper triangular
           IF( Neib>J .and. Neib/=0 )Then

           !Part 18:
            Call Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                    &,MR,k,Face,Jacobv)
           !Part 19:
            Call Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face,eigenMatrix)

!!!!       !=====================================================================================================
!          !Part 20.1: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using Increment scheme
!           Call Increment(GM,NX(Face),NY(Face),NZ(Face) ,WNP1(1,Neib),WNP1(2,Neib),&
!                               WNP1(3,Neib),WNP1(4,Neib),WNP1(5,Neib),DW(1,Neib),&
!                               DW(2,Neib),DW(3,Neib),DW(4,Neib),DW(5,Neib),Jacob_Neib)
!!!!       !-----------------------------------------------------------------------------------------------------
!          !Part 20.2: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using analytical Jacobian
            Call Calculate_Jacobi_neighb(Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Face,Neib,WNP1,WB,Jaco)

            Do kk=1,NeqI
               Jacob_Neib(kk)=0.0
               Do kkk=1,NeqI
                  Jacob_Neib(kk) = Jacob_Neib(kk) + Jaco(kk,kkk)*DW(kkk,Neib)*DA(Face)
               End do
            End do
!!!!       !=====================================================================================================
           !Part 21:
            Do kk=1,NeqI
               Varj(kk)=0.0
               Do kkk=1,NeqI
                  Varj(kk) = Varj(kk) + (0.5*(Jacobv(kk,kkk)+eigenMatrix(kk,kkk))*DA(Face))*DW(kkk,Neib)
               End Do
            End Do

           !Part 22:
            Do kk=1,NeqI
               Var_star(kk) = Var_star(kk) + 0.5*Jacob_Neib(kk)-Varj(kk)
            End Do
           End if !End upper triangular
        End Do

       !Part 23:
        Jaco(:,:)=0.0
        DD = (IsCN*Vol(j))/DT(j)

       !Part 24: !main diagonal
        DO k=1,NnonzeroCell(j)

          !Part 25:
           Face = InonzeroFace(k,j)
           Neib = InonzeroCell(k,j)

          !Part 26:
           Call Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                    &,MR,k,Face,Jacobv)

          !Part 27:
           Call Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face,eigenMatrix)

          !Part 28:
           Jaco(:,:)=Jaco(:,:)+0.5*DA(Face)*(eigenMatrix(:,:)+Jacobv(:,:))
        End do !End main diagonal

       !Part 29:
        Do k=1,NeqI
            Jaco(k,k)=Jaco(k,k)+DD
        End Do

       !Part 30:
        Call M55INV2 (Jaco, INV)

       !Part 31:
        Do kk=1,NeqI
            Var(1)=0.0
            Do k=1,NeqI
                Var(1) =Var(1) + Inv(kk,k)*( -Var_star(k)-IsCN*Res(k,J) )
           End Do
           DW_star(kk,j) = Var(1)
        End do
     End do !End Forward Sweep


    !Part 32: backward sweep
     Do J=NC,1,-1

       !Part 33:
        Do k=1,NeqI
           Var_star(k)=0.0
           Var     (k)=0.0
        End do

       !Part 34:
        DO k=1,NnonzeroCell(j)

          !Part 35:
           Neib = InonzeroCell(k,j)
           Face = InonzeroFace(k,j)

          !Part 36: !lower triangular
           IF( Neib<J .and.  Neib/=0 )then

           !Part 37:
            Call Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                    &,MR,k,Face,Jacobv)

           !Part 38:
            Call Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face,eigenMatrix)

!!!!       !=====================================================================================================
!          !Part 39.1: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using Increment scheme
!           Call Increment(GM,NX(Face),NY(Face),NZ(Face) ,WNP1(1,Neib),WNP1(2,Neib),&
!                               WNP1(3,Neib), WNP1(4,Neib), WNP1(5,Neib),DW_star(1,Neib),DW_star(2,Neib),&
!                               DW_star(3,Neib),DW_star(4,Neib),DW_star(5,Neib),Jacob_Neib)
!!!!       !-----------------------------------------------------------------------------------------------------
!          !Part 39.2: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using analytical Jacobian
            Call Calculate_Jacobi_neighb(Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Face,Neib,WNP1,WB,Jaco)

            Do kk=1,NeqI
               Jacob_Neib(kk)=0.0
               Do kkk=1,NeqI
                  Jacob_Neib(kk) = Jacob_Neib(kk) + Jaco(kk,kkk)*DW_star(kkk,Neib)*DA(Face)
               End do
            End do
!!!!       !=====================================================================================================
           !Part 40:
            Do kk=1,NeqI
               Varj(kk)=0.0
               Do kkk=1,NeqI
                  Varj(kk) = Varj(kk) + (0.5*(Jacobv(kk,kkk)+eigenMatrix(kk,kkk))*DA(Face))*DW_star(kkk,Neib)
               End Do
            End Do

           !Part 41:
            Do kk=1,NeqI
               Var(kk) = Var(kk) + 0.5*Jacob_Neib(kk) -Varj(kk)
            End do
           Endif !End lower triangular

          !Part 42: !upper triangular
           IF( Neib>J .and.  Neib/=0 )then

           !Part 43:
            Call Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                    &,MR,k,Face,Jacobv)

           !Part 44:
            Call Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face,eigenMatrix)

!!!!       !=====================================================================================================
!          !Part 45.1: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using Increment scheme
!           Call Increment(GM,NX(Face),NY(Face),NZ(Face) ,WNP1(1,Neib),WNP1(2,Neib),&
!                               WNP1(3,Neib), WNP1(4,Neib), WNP1(5,Neib),DW(1,Neib),DW(2,Neib),&
!                               DW(3,Neib),DW(4,Neib),DW(5,Neib),Jacob_Neib)
!!!!       !-----------------------------------------------------------------------------------------------------
!          !Part 45.2: !Calculate (Partial Inv. Flux)/(Partial Conserv. Var. Neib. cell) using analytical Jacobian
            Call Calculate_Jacobi_neighb(Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Face,Neib,WNP1,WB,Jaco)

            Do kk=1,NeqI
               Jacob_Neib(kk)=0.0
               Do kkk=1,NeqI
                  Jacob_Neib(kk) = Jacob_Neib(kk) + Jaco(kk,kkk)*DW(kkk,Neib)*DA(Face)
               End do
            End do
!!!!       !=====================================================================================================

           !Part 46:
            Do kk=1,NeqI
               Varj(kk)=0.0
               Do kkk=1,NeqI
                  Varj(kk) = Varj(kk) + (0.5*(Jacobv(kk,kkk)+eigenMatrix(kk,kkk))*DA(Face))*DW(kkk,Neib)
               End Do
            End Do

           !Part 47:
            Do kk=1,NeqI
               Var(kk) = Var(kk) + 0.5*Jacob_Neib(kk) - Varj(kk)
            End do
           Endif !End upper triangular
        End do

       !Part 48:
        Jaco(:,:)=0.0
        DD = (IsCN*Vol(j))/DT(j)

       !Part 49: !main diagonal
        DO k=1,NnonzeroCell(j)
          !Part 50:
           Face = InonzeroFace(k,j)
           Neib = InonzeroCell(k,j)

          !Part 51:
           Call Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                            &,MR,k,Face,Jacobv)

          !Part 52:
           Call Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face,eigenMatrix)

          !Part 53:
           Jaco(:,:)=Jaco(:,:)+0.5*DA(Face)*(eigenMatrix(:,:)+Jacobv(:,:))
        End do !End main diagonal

       !Part 54:
        Do k=1,NeqI
            Jaco(k,k)=Jaco(k,k)+DD
        End Do

       !Part 55:
        Call M55INV2 (Jaco, Inv)

       !Part 56:
        Do kk=1,NeqI
           !Part 57:
            Var_star(1)=0.0
            Do k=1,NeqI
                Var_star(1) =Var_star(1)+  Inv(kk,k)* (Var(k)-IsCN*Res(k,J))
            End Do

           !Part 58:
            Er=max(Er,abs(DW(kk,j)-(Var_star(1)))/abs(DQ(kk,j)))

           !Part 59:
            DW(kk,j)=Var_star(1)
        End do
     End do !End backward sweep
 End Do !End Error loop
!*******************************************************************************************
 End
!###########################################################################################
