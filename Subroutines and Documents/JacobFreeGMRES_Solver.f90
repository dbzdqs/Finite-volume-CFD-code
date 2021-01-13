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
!// Date: Feb., 20, 2016                                                                   //!
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine JacobFreeGMRES_Solver(Dim,Neq,N_outer,N_inner,N_used,N_limit,S_Pre,P_Pre,&
	                              NC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NX,NY,A,DA,&
                                  DT,Minf,GM,R0,P0,C0,U0,V0,&
					              WNP1,B,X) 
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Neq,N_outer,N_inner,N_limit,&
	            NC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NX,NY,A,DA,&
                DT,Minf,GM,R0,P0,C0,U0,V0,WNP1,B
 Intent(InOut)::n_used,S_Pre,P_Pre,X

 Integer::Dim,Neq,NC,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NF,NIR,N_outer,N_inner,n_used,n_limit,ST
 Integer::Flag,l_iter,I,J,J_iter,I_iter,k_iter
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8)::GM,eps,Betta,Landa,Convergence_Check,Temp,Inner_Product,norm
 Real(8)::Minf,R0,T0,P0,C0,U0,V0,E0
 Real(8),Dimension(1:Dim)::NX,NY,A,DA,DT
 Real(8),Dimension(1:4,1:Dim)::B,X,W,r,bb,XX,Ax
 Real(8),Dimension(1:4,1:Dim)::WNP1 
 Real(8),Dimension(1:n_limit,1:4,1:NC)::P_Pre,S_Pre
 Real(8),Dimension(1:(N_inner+1))::P,y,c,s
 Real(8),Dimension(1:(N_inner+1),1:(N_inner))::H
 Real(8),Dimension(1:(N_inner+1),1:4,1:Dim)::z,v    
!*********************************************************************************************	    
    
 eps = 10e-7
 X = 0.0
 
!Part 1: GMRES outer iterations
 Do L_iter=1,N_outer

    H = 0.0    
   !Part 2: Matrix_Vector_Product  A.x0 = w  input:A=dF(W)/dW , x0=X; output: W
    Call JacobianVector_Product(Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,A,&
                                GM,U0,V0,P0,R0,C0,DT,&
								WNP1,X,Ax)      

   !Part 3: compute initial residual norm   r0=B-A.x0
    Betta = 0.0
    Do I=1,NC
       Do J=1,4
         !compute initial residual: r = b - A.X = Con - A.X= B - W
          r(J,I) = B(J,I) - Ax(J,I)

         !Betta = ||r||2 
          Betta = Betta + r(J,I)*r(J,I)    
       End Do  
    End Do
    Betta = SQRT(Betta) 


   !Part 4: define first Krylov vector :  v1 = r / Betta
    Do I=1,NC
       Do J=1,4
          v(1,J,I) = r(J,I) / Betta   
       End Do  
    End Do

   !Part 5: initialize right hand side : P = Betta . e1
    P = 0.0
    P(1) = Betta
     
   !Part 6: GMRES inner iteration
    J_iter = 0
    Do while( J_iter < N_inner  )
    J_iter = J_iter + 1
          
        
      !Part 7: preconditioning step :  z = ( M^-1 ) . v
       Do I=1,NC
          Do J=1,4  
            !DelX = z : Output of preconditioneng system is Delx and goes to Z            
             bb(J,I) =  V(J_iter,J,I)
          End Do  
       End Do

       Call GMRES_Inner(Dim,Neq,N_outer,5,N_used,N_limit,S_Pre,P_Pre,&
	                    NC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NX,NY,A,DA,&
                        DT,Minf,GM,R0,P0,C0,U0,V0,&
					    WNP1,bb,XX)


       Do I=1,NC
          Do J=1,4   
            !DelX = z : Output of preconditioneng system is Delx and goes to Z            
             z(J_iter,J,I) = XX(J,I)
          End Do  
       End Do    
      
       !Part 8: matrix-vector product : w = A . z , w = A . Temp_arry2
       Call JacobianVector_Product(Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,A,&
                                   GM,U0,V0,P0,R0,C0,DT,&
								   WNP1,XX,W)                   

      !Part 9: Gramm-Schmidt orthogonalization : 
       Do I_iter=1, J_iter
         
         !hij = W . vi 
          Inner_Product = 0.0        
          Do I=1,NC
             Do J=1,4
                 Inner_Product = Inner_Product + W(J,I)*v(I_iter,J,I)                   
             End Do 
          End Do
           
         !hij = W . vi 
          H(I_iter,J_iter) = Inner_Product
           
         !W = W - hij . V
          Do I=1,NC
             Do J=1,4
                W(J,I) = W(J,I) - H(I_iter,J_iter)*v(I_iter,J,I)                
             End Do  
          End Do         
                 
       End Do
        
      !Part 10: hj+1,j = ||w||2
       norm = 0.0
       Do I=1,NC
          Do J=1,4
             norm = norm + W(J,I)*W(J,I) 
          End Do  
       End Do
       H(J_iter+1,J_iter) = Sqrt(norm)
        
      !Part 11: define next Krylov vector : Vj+1 = W/hj+1,j
       Do I=1,NC
          Do J=1,4
             v(J_iter+1,J,I) = W(J,I)/H(J_iter+1,J_iter)
          End Do 
       End Do    
        
       
      !Part : previous Givens rotations on H : 
       Do I_iter=1, J_iter-1  
          Temp= c(I_iter)*H(I_iter,J_iter) + s(I_iter)*H(I_iter+1,J_iter)
          H(I_iter+1,J_iter) = -s(I_iter)*H(I_iter,J_iter) + c(I_iter)*H(I_iter+1,J_iter)
          H(I_iter,J_iter) = Temp
       End Do
        
      !Part 12: compute next rotation : 
       Landa = SQRT ( ( H(J_iter,J_iter)**2.0 ) + ( H(J_iter+1,J_iter)**2.0 ) )
       C(J_iter) = H(J_iter,J_iter) / Landa
       S(J_iter) = H(J_iter+1,J_iter) / Landa

        
      !Part 13: Givens rotation on H:
       H(J_iter,J_iter) = Landa
       H(J_iter+1,J_iter) = 0.0
        
      !Part 14: Givens rotation on P:
       P(J_iter+1) = -S(J_iter)*P(J_iter)
       P(J_iter) = C(J_iter)*P(J_iter)    
        
      !Part 15:inner loop convergence check    
       Convergence_Check = P(J_iter+1)
        
    End Do
     
     
   !Part 16: Back substitution :
    y(J_iter) = P(J_iter) / H(J_iter,J_iter)
    Do k_iter =J_iter-1,1,-1 
       Temp = 0.0 
       Do I_iter = k_iter+1 , J_iter
          Temp = Temp + H(k_iter,I_iter) * y(I_iter)
       End Do
       y(k_iter) = ( P(k_iter) - Temp  ) / H(k_iter,k_iter)
    End Do
     
     
   !Part 17: form approximate solution :
    Do I=1,NC   
       Do J=1,4

          Temp = 0.0
          Do I_iter=1,J_iter            
             Temp = Temp + y(I_iter)*Z(I_iter,J,I)             
          End Do
          X(J,I) = X(J,I) + Temp 

       End Do
    End Do   

     
   !Part 18:inner loop convergence check
   !IF (ABS( Convergence_Check ) <= eps ) Flag = 0    
   
 End Do
!*********************************************************************************************
 End
!###########################################################################################


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:                                                                         //!
!//                                                                                      //!
!// Version:                                                                             //!
!// Date:                                                                                //!
!// Developed by:                                                                        //!
!// Doc ID:                                                                              //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine GMRES_Inner(Dim,Neq,N_outer,N_inner,N_used,N_limit,S_Pre,P_Pre,&
	                    NC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NX,NY,A,DA,&
                        DT,Minf,GM,R0,P0,C0,U0,V0,&
					    WNP1,B,X)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Neq,N_outer,N_inner,N_limit,&
	            NC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NX,NY,A,DA,&
                DT,Minf,GM,R0,P0,C0,U0,V0,B,WNP1
 Intent(InOut)::n_used,S_Pre,P_Pre,X

 Integer::Dim,Neq,NC,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NF,NIR,N_outer,N_inner,n_used,n_limit,FF
 Integer::Flag,l_iter,I,J,J_iter,I_iter,k_iter
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8)::GM,eps,Betta,Landa,Convergence_Check,Temp,Inner_Product,norm
 Real(8)::Minf,R0,T0,P0,C0,U0,V0,E0
 Real(8),Dimension(1:Dim)::NX,NY,A,DA,DT
 Real(8),Dimension(1:4,1:Dim)::B,X,W,r,TempZ
 Real(8),Dimension(1:4,1:Dim)::WNP1 
  Real(8),Dimension(1:n_limit,1:4,1:NC)::P_Pre,S_Pre
 Real(8),Dimension(1:(N_inner+1))::P,y,c,s
 Real(8),Dimension(1:(N_inner+1),1:(N_inner+1))::H
 Real(8),Dimension(1:(N_inner+1),1:4,1:Dim)::z,v    
!*********************************************************************************************	    
 
 eps = 10e-7
 X = 0.0

!Part 1: GMRES outer iterations
 Do L_iter=1,N_outer
 
     H = 0.0    
    !Part 2: Matrix_Vector_Product  A.x0 = w  input:A=dF(W)/dW , x0=X; output: W
     Call JacobianVector_Product(Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,A,&
                                 GM,U0,V0,P0,R0,C0,DT,&
	 						     WNP1,X,W)       
                      


   !Part 3: compute initial residual norm   r0=B-A.x0
    Betta = 0.0
    Do I=1,NC
       Do J=1,4
         !compute initial residual: r = b - A.X = Con - A.X= B - W
          r(J,I) = B(J,I) - W(J,I)
         
		 !Betta = ||r||2 
          Betta = Betta + ( r(J,I)**2.0 )    
       End Do  
    End Do
    Betta = SQRT(Betta) 
  
   !Part 4: define first Krylov vector :  v1 = r / Betta
    Do I=1,NC
       Do J=1,4
          v(1,J,I) = r(J,I) / Betta   
       End Do 
    End Do

   !Part 5: initialize right hand side : P = Betta . e1
    P = 0.0
    P(1) = Betta
     
   !Part 6: GMRES inner iteration
     J_iter = 0
     Do while ( J_iter < N_inner  )
     J_iter = J_iter + 1
           
      !Part 7: preconditioning step :  z = ( M^-1 ) . v 
       Call Preconditioning(Dim,NC,N_inner,J_iter,n_used,n_limit,S_Pre,P_Pre,v,z)
   
       Do I=1,NC
          Do J=1,4  
            !DelX = z : Output of preconditioneng system is Delx and goes to Z            
             TempZ(J,I) =  z(J_iter,J,I)
          End Do  
       End Do    

      !Part 8: matrix-vector product : w = A . z , w = A . Temp_arry2
       Call JacobianVector_Product(Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,A,&
                                   GM,U0,V0,P0,R0,C0,DT,&
								   WNP1,TempZ,W)                   
 

      !Part 9: Gramm-Schmidt orthogonalization : 
       Do I_iter=1,J_iter
         
         !hij = W . vi 
          Inner_Product = 0.0        
          Do I=1,NC
             Do J=1,4
                Inner_Product = Inner_Product + W(J,I)*v(I_iter,J,I)                   
             End Do  
          End Do
           
           
         !hij = W . vi 
          H(I_iter,J_iter) = Inner_Product
           
         !W = W - hij . V
          Do I=1,NC
             Do J=1,4
                W(J,I) = W(J,I) - H(I_iter,J_iter)*v(I_iter,J,I)                
             End Do  
          End Do         
               
       End Do
        
      !Part 10: hj+1,j = ||w||2
       norm = 0.0
       Do I=1,NC
          Do J=1,4
             norm = norm + ( W(J,I)**2.0 ) 
          End Do 
       End Do
       H(J_iter+1,J_iter) = Sqrt(norm)
       
 

      !Part 11: define next Krylov vector : Vj+1 = W/hj+1,j
       Do I=1,NC
          Do J=1,4
             v(J_iter+1,J,I) = W(J,I)/H(J_iter+1,J_iter)
          End Do 
       End Do    
        
       
      !Part : previous Givens rotations on H : 
       Do I_iter=1, J_iter-1     
          Temp= c(I_iter)*H(I_iter,J_iter) + s(I_iter)*H(I_iter+1,J_iter)
          H(I_iter+1,J_iter) = -s(I_iter)*H(I_iter,J_iter) + c(I_iter)*H(I_iter+1,J_iter)
          H(I_iter,J_iter) = Temp

       End Do
             

      !Part 12: compute next rotation : 
       Landa = SQRT ( ( H(J_iter,J_iter)**2.0 ) + ( H(J_iter+1,J_iter)**2.0 ) )
       C(J_iter) = H(J_iter,J_iter) / Landa
       S(J_iter) = H(J_iter+1,J_iter) / Landa
   
      !Part 13: Givens rotation on H:
       H(J_iter,J_iter) = Landa
       H(J_iter+1,J_iter) = 0.0
        
      !Part 14: Givens rotation on P:
       P(J_iter+1) = -S(J_iter)*P(J_iter)
       P(J_iter) = C(J_iter)*P(J_iter)    
       
      !Part 15:inner loop convergence check    
       Convergence_Check = P(J_iter+1)
        
    End Do
   
         

   !Part 16: Back substitution :
    y(J_iter) = P(J_iter) / H(J_iter,J_iter)
     
    Do k_iter=J_iter-1,1,-1     
     
       Temp = 0.0 
       Do I_iter = k_iter+1 , J_iter    
          Temp = Temp + H(k_iter,I_iter) * y(I_iter)
       End Do
       y(k_iter) = ( P(k_iter) - Temp  ) / H(k_iter,k_iter)
        
    End Do
 
   !Part 17: form approximate solution :
    Do I=1,NC 
       Do J = 1, 4

          Temp = 0.0
          Do I_iter = 1 , J_iter            
             Temp = Temp + y(I_iter)*Z(I_iter,J,I)            
          End Do
           
          X(J,I) = X(J,I) + Temp 
           
       End Do
    End Do   

   !Part 18:inner loop convergence check
   !IF (ABS( Convergence_Check ) <= eps ) Flag = 0    

 End Do
!*********************************************************************************************
 End
!###########################################################################################

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:                                                                         //!
!//                                                                                      //!
!// Version:                                                                             //!
!// Date:                                                                                //!
!// Developed by:                                                                        //!
!// Doc ID:                                                                              //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Preconditioning(Dim,NC,m_in,J_iter,n_used,n_limit,S_Pre,P_Pre,v,z)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,m_in,J_iter,n_limit,S_Pre,P_Pre,v
 Intent(InOut)::n_used,z

 Integer::Dim,NC,m_in,J_iter,n_used,n_limit,Column_N,I,J,JJ,Temp  
 Real(8) ::PZ
 Real(8),Dimension(1:(m_in+1),1:4,1:Dim)::v,z       
 Real(8),Dimension(1:n_limit,1:4,1:NC)::P_Pre,S_Pre
!*********************************************************************************************
!Part 1:      
 Do I =1,NC      
    Do J=1,4    
       z(J_iter,J,I) = v(J_iter,J,I)        
    End Do                
 End Do

 IF( n_used > 0 )Then
 
 !Part 2:
  Column_N = mod(n_used,n_limit)
  IF (Column_N == 0) Column_N = n_limit
  
 !Part 3:
  Do JJ=1,n_limit
    
    !Part 4:
     Temp = mod(Column_N+JJ,n_limit)
     IF (Temp == 0) Temp = n_limit
    
    !Part 5:
     PZ = 0.0
     Do I=1,NC 
        Do J= 1 , 4
           PZ = PZ + P_Pre(Temp,J,I) * z(J_iter,J,I)
        End Do
     End Do
    
    !Part 6: 
     Do I = 1,NC 
         Do J= 1 , 4
            z(J_iter,J,I) = z(J_iter,J,I) + (PZ*S_Pre(Temp,J,I))
        End Do
     End Do    
    
  End Do
     
 End IF   
!********************************************************************************************* 
 End
!###########################################################################################

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:                                                                         //!
!//                                                                                      //!
!// Version:                                                                             //!
!// Date:                                                                                //!
!// Developed by:                                                                        //!
!// Doc ID:                                                                              //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Preconditioner_Update(Dim,NC,n_used,n_limit,DelX,DelF,S_Pre,P_Pre)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,n_limit,DelX,DelF
 Intent(InOut)::n_used,S_Pre,P_Pre

 Integer::Dim,NC,n_used,n_limit,Column_N,I,J,JJ,Temp  
 Real(8) :: Kesi,PT
 
 Real(8),Dimension(1:4,1:Dim)::T_Pre,DelX,DelF        
 Real(8),Dimension(1:n_limit,1:4,1:NC)::P_Pre,S_Pre
!*********************************************************************************************	 
!Part 1:
  n_used = n_used + 1 
  Column_N = mod(n_used,n_limit)
  IF (Column_N == 0) Column_N = n_limit
  
  
  !Part 2:
  Do I = 1,NC
    Do J= 1 , 4
        T_Pre(J,I)= DelF(J,I)
    End Do 
  End Do    
  
  !Part 3:
  Do JJ = 1 , n_limit
    
    !Part 4:
    Temp = mod(Column_N+JJ,n_limit)
    IF (Temp == 0) Temp = n_limit
    
    !Part 5:
    PT = 0.0
    Do I = 1,NC 
       Do J= 1 , 4
         PT = PT + P_Pre(Temp,J,I)*T_Pre(J,I)
       End Do
    End Do
    
    !Part 6: 
    Do I = 1,NC 
       Do J= 1 , 4
         T_Pre(J,I) = T_Pre(J,I) + (PT*S_Pre(Temp,J,I))
       End Do
    End Do    
     
  End Do
  
  !Part 7:
  Do I = 1,NC
    Do J= 1 , 4
        S_Pre(Column_N,J,I) = DelX(J,I) - T_Pre(J,I) 
    End Do 
  End Do     
  
  
  !Part 8:
  Kesi = 0.0
  Do I = 1,NC
    Do J= 1 , 4
        Kesi = Kesi + DelX(J,I)*T_Pre(J,I)
    End Do 
  End Do
  
  !Part 9:
  Do I = 1,NC
    Do J= 1 , 4
        P_Pre(Column_N,J,I) = DelX(J,I)/Kesi
    End Do 
  End Do  
!********************************************************************************************* 
 End
!###########################################################################################