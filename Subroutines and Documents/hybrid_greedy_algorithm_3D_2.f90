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
!// Developed by: M. A. Zoljanahi, Mechanical Eng., Amirkabir University of Technology     //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine hybrid_greedy_algorithm_3D_2(NBP,Dim,X,Y,Z,IBP,Delx,Dely,Delz,NSBP,ACTV,coxx,coyy,cozz)
 Implicit None
!*********************************************************************************************
 Intent(In)::NBP,Dim,X,Y,Z,IBP,Delx,Dely,Delz
 Intent(Out)::NSBP,ACTV,coyy,coxx,cozz
 Integer::dim,NBP,NSBP,j,i,imax,Itere,imax2,l
 Integer,Parameter::m=9
 Real,Parameter::epsilon=0.0002
 Real(8)::Dx,Dy,Dz,DL,Temp,Sum_DelX,Sum_DelY,Sum_DelZ,Xi,Yi,Zi,fmax,fxmax,fymax,fzmax,betax,betay,betaz,Emax,ffmax,Eave
 Real(8),Dimension(Dim)::X,Y,Z,DelX,DelY,DelZ,Delevx,Delevy,Delevz,E,Ex,Ey,Ez,fevx,fevy,fevz,EE,Exx,Eyy,Ezz,b_y,b_x,b_z
 Integer,Dimension(1:Dim)::IBP
 Integer,Dimension(1:NBP)::ACTV
 Real(8),Dimension(NBP)::Cox,Coy,Coz
 Real(8),Dimension(NBP)::Coxx,Coyy,Cozz
 Real(8),Dimension(1:NBP,1:NBP)::B
!************************************************************************************************************************************************************
!part 1:
 imax=1
 Itere=0
 NSBP=0
 fxmax=Delx(IBP(imax))
 fymax=Dely(IBP(imax))
 fzmax=Delz(IBP(imax))
 Emax=100
 
 !part 2:
 Do While (Emax>(5.e-4))
     
     Itere=Itere+1
     
    !part 3:
     Dl=0
     Call RBF_Function3D(DL,Temp)
     betax=fxmax/Temp
     betay=fymax/Temp
     betaz=fzmax/Temp
     
     NSBP=NSBP+1
     
    !part 4:
     If (NSBP>NBP) Then
         NSBP=NSBP-1
         Goto 9
     End If
     
     ACTV(NSBP)=IBP(imax)
     
    !part 5:
     Do I=1,NSBP
         b_x(I)=DelX(ACTV(I))
         b_y(I)=DelY(ACTV(I))
         b_z(I)=DelZ(ACTV(I))
         
     End Do
     
    !Part 6:
     Call RBF_Coefficient_Matrix_3D(Dim,X,Y,Z,ACTV,NSBP,B)
     
    !Part 7:
     Call solve_lu(NSBP,B,b_x,cox )
     Call solve_lu(NSBP,B,b_y,coy )
     Call solve_lu(NSBP,B,b_z,coz )
     
    !part 8:
     Do i=1,NSBP
         coxx(i)=cox(i)+betax
         coyy(i)=coy(i)+betay
         cozz(i)=coz(i)+betaz
     End Do
     
    !part 9:
     fmax=0
     l=0
     Eave=0.0
     
    !part 10:
     Do i=1,NBP
        !part 11:
         Xi = X(IBP(i))
	     Yi = Y(IBP(i))
         Zi = Z(IBP(i))
         
        !Part 12:
         Sum_DelX = 0.0
	     Sum_DelY = 0.0
         Sum_DelZ = 0.0
         
        !Part 13:            
         Do J=1,NSBP
             
            !part 14:
             If(IBP(i)==ACTV(j)) Goto 7
             
            !Part 15:
             Dx = X( ACTV(j) )-Xi
             Dy = Y( ACTV(j) )-Yi
             Dz = Z( ACTV(j) )-Zi
             DL = Dsqrt( Dx*Dx + Dy*Dy + Dz*Dz )
             
            !Part 16:       
             Call RBF_Function3D(DL,Temp)
             
            !Part 17:
             Sum_DelX = Sum_DelX + Cox(j)*Temp   
             Sum_DelY = Sum_DelY + Coy(j)*Temp
             Sum_DelZ = Sum_DelZ + Coz(j)*Temp
         End Do
         
        !part 18:
         Dx = X( IBP(imax) )-Xi
         Dy = Y( IBP(imax) )-Yi
         Dz = Z( IBP(imax) )-Zi
         DL = Dsqrt( Dx*Dx + Dy*Dy + Dz*Dz )        
         Call RBF_Function3D(DL,Temp)
         
        !part 19:  
         delevx(IBP(i))=Sum_DelX+(betax*Temp)
         delevy(IBP(i))=Sum_DelY+(betay*Temp)
         delevz(IBP(i))=Sum_DelZ+(betaz*Temp)
         
        !part 20:
         Ex(i)=(DelX(IBP(i))-delevx(IBP(i)))
         Ey(i)=(DelY(IBP(i))-delevy(IBP(i)))
         Ez(i)=(DelZ(IBP(i))-delevz(IBP(i)))
         E(i)=Dsqrt( Ex(i)**2 + Ey(i)**2 + Ez(i)**2)
         
        !part 21:
         If ( E(i)>fmax ) Then
             fmax=E(i)
             imax=i
             fevx(IBP(imax))=Sum_DelX
             fevy(IBP(imax))=Sum_DelY
             fevz(IBP(imax))=Sum_DelZ
         End If
         
        !part 22:
         l=l+1
         Eave=(Eave+E(i))/l
         
7    End Do
     
    !part 23:
     IF (Mod(itere,m)==0) Then
         
        !part 24:    
         Do i=1,NBP
             Do j=1,NSBP
                 If(IBP(i)==ACTV(j)) Goto 6
             End Do
             
            !part 25: 
             If (E(i)>(fmax-epsilon)) Then
                 imax2=i
                 NSBP=NSBP+1
                 
                !part 26:
                 If (NSBP>NBP) Then
                     NSBP=NSBP-1
                     Goto 9
                 End If
                 
                 ACTV(NSBP)=IBP(imax2)
             End If
             
6        End Do
         
        !part 27:
         Do I=1,NSBP
             
             b_x(I)=DelX(ACTV(I))
             b_y(I)=DelY(ACTV(I))
             b_z(I)=DelZ(ACTV(I))
             
         End Do
         
        !Part 28:
         Call RBF_Coefficient_Matrix_3D(Dim,X,Y,Z,ACTV,NSBP,B)
         
        !Part 29:
         Call solve_lu(NSBP,B,b_x,cox )
         Call solve_lu(NSBP,B,b_y,coy )
         Call solve_lu(NSBP,B,b_z,coz )
         
        !part 30:
         fmax=0
         Do i=1,NBP
             
            !part 31:
             Xi = X(IBP(i))
	         Yi = Y(IBP(i))
             Zi = Z(IBP(i))
             
            !Part 32:
             Sum_DelX = 0.0
	         Sum_DelY = 0.0
             Sum_DelZ = 0.0
             
            !Part 33:            
             Do J=1,NSBP
                 
                !part 34:
                 If(IBP(i)==ACTV(j)) Goto 88
                 
                !Part 35:
                 Dx = X( ACTV(j) )-Xi
                 Dy = Y( ACTV(j) )-Yi
                 Dz = Z( ACTV(j) )-Zi
                 DL = Dsqrt( Dx*Dx + Dy*Dy + Dz*Dz )
                 
                !Part 36:        
                 Call RBF_Function3D(DL,Temp)
                 
                !Part 37:    
                 Sum_DelX = Sum_DelX + Cox(j)*Temp   
                 Sum_DelY = Sum_DelY + Coy(j)*Temp
                 Sum_DelZ = Sum_DelZ + Coz(j)*Temp
             End Do
             
            !part 38:
             delevx(IBP(i))=Sum_DelX
             delevy(IBP(i))=Sum_DelY
             delevz(IBP(i))=Sum_DelZ
             
            !part 39: 
             Ex(i)=(DelX(IBP(i))-delevx(IBP(i)))
             Ey(i)=(DelY(IBP(i))-delevy(IBP(i)))
             Ez(i)=(DelZ(IBP(i))-delevz(IBP(i)))
             E(i)=Dsqrt( Ex(i)**2 + Ey(i)**2 + Ez(i)**2)
             
            !part 40:
             If ( E(i)>fmax ) Then
                 fmax=E(i)
                 imax=i
                 fevx(IBP(imax))=Sum_DelX
                 fevy(IBP(imax))=Sum_DelY
                 fevz(IBP(imax))=Sum_DelZ
             End If
             
88       End Do
         
     End If

    !part 42:
     Emax=fmax
     
    !part 43:
     fxmax=(Delx(IBP(imax))-fevx(IBP(imax)))
     fymax=(Dely(IBP(imax))-fevy(IBP(imax)))
     fzmax=(Delz(IBP(imax))-fevz(IBP(imax)))
     
 End Do
!*********************************************************************************************
9 End     
!###########################################################################################     
 
 
 