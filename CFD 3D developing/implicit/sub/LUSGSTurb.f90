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
 Subroutine LUSGSTurb(Dim,NC,NX,NY,NZ,DA,GM,Rest,NnonzeroCell,InonzeroCell,InonzeroFace,eigen,JacobiL,JacobiR,Vol,DT,WTNP1,DW)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NC,NX,NY,NZ,DA,GM,Rest,NnonzeroCell,InonzeroCell,InonzeroFace,eigen,JacobiL,JacobiR,Vol,DT,WTNP1
 Intent (Out  )::DW

 Integer::Dim,I,II,J,NS,ME
 Integer::NC
 Integer::KK
 Real(8)::MR
 Real(8)::GM,sor,DD,IsCN
 Real(8),Dimension(1:2,1:Dim)::WTNP1
 Real(8),Dimension(1:Dim)::Vol
 Real(8),Dimension(1:Dim)::DT
 Real(8),Dimension(1:2,1:Dim)::Rest
 Real(8),Dimension(1:Dim)::NX,NY,NZ
 Real(8),Dimension(1:Dim)::DA
 Integer,Dimension(1:10,1:Dim)::InonzeroCell,InonzeroFace
 Integer,Dimension(1:Dim)::NnonzeroCell
 Real(8),Dimension(1:Dim)::eigen,JacobiL,JacobiR

 Integer::Neib,Face,NeqI,k
 Real(8),Dimension(1:2,1:Dim)::DW_star,DW
!*******************************************************************************************
!Part 1: !If IsCN=2.0 => Second-order Crank-Nicelson scheme. If IsCN=1.0 => First-order Euler scheme.
 IsCN=2.0

!Part 2:
 NeqI=2

!Part 3: forward sweep
 Do J=1,NC

   !Part 4:
    Do k=1,NeqI
       DW_star(k,j)=0.0
       DW     (k,j)=0.0
    End do

   !Part 5:
    DO k=1,NnonzeroCell(j)

      !Part 6:
       Neib = InonzeroCell(k,j)
       Face = InonzeroFace(k,j)

      !Part 7:
       IF( Neib<J .and. Neib/=0 )Then

       !Part 8:
        Do kk=1,NeqI
           DW_star(kk,j) = DW_star(kk,j) + 0.5*( DA(Face)*(JacobiR(Face)-eigen(Face))*DW_star(kk,Neib) )
        End Do

       End if
    End Do

   !Part 9:
    DD = 0.0
    DO k=1,NnonzeroCell(j)
       Face = InonzeroFace(k,j)
       DD = DD + 0.5*(JacobiL(Face)+eigen(Face))*DA(Face)
    End do
	DD=DD+IsCN*Vol(j)/DT(j)

   !Part 10:
    Do kk=1,NeqI
       DW_star(kk,j) = ( DW_star(kk,j)- IsCN*Rest(kk,j) )/DD
    End do

 End do


!Part 11: backward sweep
 Do J=NC,1,-1

   !Part 12:
    Do k=1,NeqI
       DW(k,j)=0.0
    End do

   !Part 13:
    DO k=1,NnonzeroCell(j)

      !Part 14:
       Neib = InonzeroCell(k,j)
       Face = InonzeroFace(k,j)

      !Part 15:
       IF( Neib<J .and.  Neib/=0 )then

       !Part 16:
        Do kk=1,NeqI
            DW(kk,j) = DW(kk,j) + 0.5*( DA(Face)*(JacobiR(Face)-eigen(Face))*DW(kk,Neib) )
        End do

       Endif
    End do

   !Part 17:
    DD=0.0
    DO k=1,NnonzeroCell(j)
       Face = InonzeroFace(k,j)
       DD = DD + 0.5*(JacobiL(Face)+eigen(Face))*DA(Face)
    End do
    DD= DD+IsCN*Vol(j)/DT(j)

   !Part 18:
    Do kk=1,NeqI
       DW(kk,j) = DW_star(kk,j) + DW(kk,j)/DD
    End do
 End do
!*******************************************************************************************
 end
!###########################################################################################
