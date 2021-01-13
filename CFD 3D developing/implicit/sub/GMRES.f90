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
 Subroutine GMRES(Dim,IDS,NC,NF2,NF,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace,&
                            Vol,P,Mu,Mut,PrL,PrT,xc,yc,zc,DT,MR,WNP1,WB,DW)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,IDS,NC,NF2,NF,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace,&
                            Vol,P,Mu,Mut,PrL,PrT,xc,yc,zc,DT,MR,WNP1,WB
 Intent (out  )::DW

 Integer::Dim,I,NC,NF1,NF2,NF,ColNum,n
 Real(8),Dimension(1:5,1:Dim)::Res
 Real(8),Dimension(1:Dim)::DT,Vol,XC,YC,ZC,NX,NY,NZ,DA,Mu,Mut,P
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:5,1:Dim)::DW
!*******************************************************************************************
 Real(8),Dimension(1:Dim)::SOURCE,x_estimate
 Real(8),Dimension(1:5,1:5)::JACO_CELL,Jaco_Neib,eigenMatrix,Jacobv
 Real(8),Dimension(1:5,1:5,1:10)::ImpMatrix
 Integer,Dimension(1:10000000)::IA,JA
 Real(8),Dimension(1:10000000)::COMP_NZ
 Integer::NUM_NZ,J,K,II,KK,JJ,Face,Neib,NeqI
 Real(8)::tol_abs,tol_rel,GM,MR,PrL,PrT,IsCN,DD
 Integer,Dimension(1:10,1:Dim)::InonzeroCell
 Integer,Dimension(1:10,1:Dim)::InonzeroFace
 Integer,Dimension(1:Dim)::NnonzeroCell,IColNum
!*******************************************************************************************
!Part 1: !If NeqI=5   => Capable of 3D.                      If NeqI=4   => Capable of 2D.
 NeqI=5

!Part 2: !If IsCN=2.0 => Second-order Crank-Nicelson scheme. If IsCN=1.0 => First-order Euler scheme.
 IsCN=2.0   !IsCN=2.0=> CN method or IsCN=1.0=> Euler

!Part 3:
 num_nZ=1
 tol_abs=0.01
 tol_rel=0.01
 x_estimate(:)=0.0001

!Part 4:
 DO I=1,NC
   !Part 5:
    ColNum=0
    IColNum(:)=0

   !Part 6:
    JACO_CELL(:,:)=0.0
    DD = (IsCN*Vol(I))/DT(I)

   !Part 7: !main diagonal
    DO k=1,NnonzeroCell(I)

      !Part 8:
       Face = InonzeroFace(k,I)
       Neib = InonzeroCell(k,I)

      !Part 9:
       Call Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                            &,MR,k,Face,Jacobv)

      !Part 10:
       Call Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face,eigenMatrix)

      !Part 11:
       JACO_CELL(:,:)=JACO_CELL(:,:)+0.5*DA(Face)*(Jacobv(:,:)+eigenMatrix(:,:))
    End do !End main diagonal

   !Part 12:
    Do k=1,NeqI
       JACO_CELL(k,k)=JACO_CELL(k,k)+DD
    End Do

   !Part 13: !Calculating Matrices of one row for one Cell
    DO k=1,NnonzeroCell(I)

       !Part 14:
        Face = InonzeroFace(k,I)
        Neib = InonzeroCell(k,I)

       !Part 15:
        If(Neib/=0)then

        !Part 16:
         ColNum=ColNum+1
         IColNum(ColNum)=Neib

        !Part 17:
         ImpMatrix(:,:,ColNum)=0.0

        !Part 18:
         If(Neib==I)then
          ImpMatrix(:,:,ColNum)=JACO_CELL(:,:)

        !Part 19:
         Else

         !Part 20:
          Call Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,&
                    &P,Face,eigenMatrix)
         !Part 21:
          Call Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                            &,MR,k,Face,Jacobv)

         !Part 22:
          Call Calculate_Jacobi_neighb(Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,face,Neib,WNP1,WB,Jaco_Neib)

         !Part 23:
          ImpMatrix(:,:,ColNum) = 0.5*Jaco_Neib(:,:)*DA(Face) - (0.5*(Jacobv(:,:)+eigenMatrix(:,:))*DA(Face))
         End if
        End if !End Non-Boundary face
    End do !End Calc. Matrices of one row

!=================================================
   !Part 24: !CR Matrix
    Do k=1,NeqI

       !Part 25:
        ia(5*(I-1)+1+(K-1))=num_nZ

       !Part 26:
        Do j=1,ColNum
            Do kk=1,NeqI
                ja(num_nZ)=5*(IColNum(j)-1)+1+(kk-1)
                COMP_NZ(num_nZ)=ImpMatrix(k,kk,j)
                num_nZ=num_nZ+1
            End Do
        End Do

       !Part 27:
        SOURCE((5*(I-1)+1)+(K-1))=-IsCN*(Res(K,I))
    END DO
!!=================================================
 END DO !Cell

!Part 28:
 IA(5*NC+1)=num_nZ
 n=5*NC

 print*,'--------------'
!Part 29:
 Call pmgmres_ilu_cr ( n, num_nZ, ia, ja, COMP_NZ, x_estimate, SOURCE, 5,5, tol_abs, tol_rel )

 print*,'----------2----'
!Part 30:
 Do J=1,NC
    DW(1,J) = x_estimate(5*(j-1)+1)
    DW(2,J) = x_estimate(5*(j-1)+2)
    DW(3,J) = x_estimate(5*(j-1)+3)
    DW(4,J) = x_estimate(5*(j-1)+4)
    DW(5,J) = x_estimate(5*(j-1)+5)
 End Do
 
End
!###########################################################################################
