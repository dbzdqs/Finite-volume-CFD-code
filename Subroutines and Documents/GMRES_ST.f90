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
!// Developed by: S. Sheikhi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine GMRES_ST(Dim,n,nz_num,ia,ja,a,x,rhs,itr_max,itr_in,tol_abs,tol_rel)
 Implicit None     
!*********************************************************************************************
 Intent(In   )::Dim,N,NZ_NUM,IA,JA,A,RHS,ITR_MAX,itr_in,TOL_ABS,TOL_REL
 Intent(InOut)::X

 Integer::n,nz_num,itr_in,i,itr,itr_max,itr_used,j,k,k_copy,Dim
 Real(8)::av(1:4),c(1:4,1:itr_in),g(1:4,1:itr_in+1),gg(1:4,1:itr_in+1),h(1:4,1:itr_in+1,1:itr_in),htmp,mu(1:4),s(1:4,1:itr_in)
 Real(8)::tol_abs,tol_rel,v(1:4,1:n,1:itr_in+1),y(1:4,1:itr_in+1),yy(1:4,1:itr_in+1),rho(1:4),rho_tol(1:4)
 Real(8), parameter :: delta = 1.0D-03
 logical, parameter :: verbose = .false.
 Real(8),Dimension(1:Dim)::a
 Real(8),Dimension(1:4,1:Dim)::x,rhs,r,r22,r23
 Integer,Dimension(1:Dim)::ia,ja
!*********************************************************************************************
!Part 1:
 itr_used = 0

!Part2:
 Do itr = 1, itr_max
print*,itr
   !Part3:
	Call AX_st(Dim,n,nz_num,ia,ja,a,x,r)

   !Part4:
	r(1:4,1:n)= rhs(1:4,1:n)-r(1:4,1:n)
	
   !Part5:
	rho(1) = sqrt ( dot_product ( r(1,1:n), r(1,1:n) ) )
	rho(2) = sqrt ( dot_product ( r(2,1:n), r(2,1:n) ) )
	rho(3) = sqrt ( dot_product ( r(3,1:n), r(3,1:n) ) )
	rho(4) = sqrt ( dot_product ( r(4,1:n), r(4,1:n) ) )
	

   !Part6:
    if ( itr == 1 ) rho_tol(1:4) = rho(1:4) * tol_rel

   !Part7:
    v(1,1:n,1) = r(1,1:n) / rho(1)
	v(2,1:n,1) = r(2,1:n) / rho(2)
	v(3,1:n,1) = r(3,1:n) / rho(3)
	v(4,1:n,1) = r(4,1:n) / rho(4)

   !Part8:
    g(1:4,1) = rho(1:4)
    g(1:4,2:itr_in+1) = 0.0D+00

   !Part9:
    h(1:4,1:itr_in+1,1:itr_in) = 0.0D+00

    Do k = 1, itr_in

	 !Part10:
      k_copy = k
	  
	  r22(1:4,1:n)=v(1:4,1:n,k)
	  call AX_st (Dim, n, nz_num, ia, ja, a, r22,r23 )
	  v(1:4,1:n,k+1)=r23(1:4,1:n)
	  
      av(1) = sqrt ( dot_product ( v(1,1:n,k+1), v(1,1:n,k+1) ) )
	  av(2) = sqrt ( dot_product ( v(2,1:n,k+1), v(2,1:n,k+1) ) )
	  av(3) = sqrt ( dot_product ( v(3,1:n,k+1), v(3,1:n,k+1) ) )
	  av(4) = sqrt ( dot_product ( v(4,1:n,k+1), v(4,1:n,k+1) ) )

	 !Part11: 
      Do j = 1, k

	   !Part12:
        h(1,j,k) = dot_product ( v(1,1:n,k+1), v(1,1:n,j) )
		h(2,j,k) = dot_product ( v(2,1:n,k+1), v(2,1:n,j) )
		h(3,j,k) = dot_product ( v(3,1:n,k+1), v(3,1:n,j) )
		h(4,j,k) = dot_product ( v(4,1:n,k+1), v(4,1:n,j) )

	   !Part13:
        v(1,1:n,k+1) = v(1,1:n,k+1) - h(1,j,k) * v(1,1:n,j)
		v(2,1:n,k+1) = v(2,1:n,k+1) - h(2,j,k) * v(2,1:n,j)
		v(3,1:n,k+1) = v(3,1:n,k+1) - h(3,j,k) * v(3,1:n,j)
		v(4,1:n,k+1) = v(4,1:n,k+1) - h(4,j,k) * v(4,1:n,j)
      end Do

	 !Part14:
	  h(1,k+1,k) = sqrt ( dot_product ( v(1,1:n,k+1), v(1,1:n,k+1) ) )
	  h(2,k+1,k) = sqrt ( dot_product ( v(2,1:n,k+1), v(2,1:n,k+1) ) )
	  h(3,k+1,k) = sqrt ( dot_product ( v(3,1:n,k+1), v(3,1:n,k+1) ) )
	  h(4,k+1,k) = sqrt ( dot_product ( v(4,1:n,k+1), v(4,1:n,k+1) ) )
	  
	  
     Do I=1,4
	    if ( av(I) + delta * h(I,k+1,k) == av(I) ) then

           Do j = 1, k
              htmp = dot_product ( v(I,1:n,k+1), v(I,1:n,j) )
              h(j,k,I) = h(j,k,I) + htmp
              v(I,1:n,k+1) = v(I,1:n,k+1) - htmp * v(I,1:n,j)
           end Do

           h(I,k+1,k) = sqrt ( dot_product ( v(I,1:n,k+1), v(I,1:n,k+1) ) )

         end if
	  	 
         if ( h(I,k+1,k) /= 0.0D+00 ) then
            v(I,1:n,k+1) = v(I,1:n,k+1) / h(I,k+1,k)
         end if
	  end Do
	  
	 !Part15:	 	  
      if ( 1 < k ) then
				
        y(1:4,1:k+1) = h(1:4,1:k+1,k)
	    Do j = 1, k - 1
          call mult_givens ( c(1:4,j), s(1:4,j), j,y(1:4,1:k+1) )
        end Do

        h(1:4,1:k+1,k) = y(1:4,1:k+1)
		
      end if

	 !Part16:
      mu(1) = sqrt ( h(1,k,k)**2 + h(1,k+1,k)**2 )
	  mu(2) = sqrt ( h(2,k,k)**2 + h(2,k+1,k)**2 )
	  mu(3) = sqrt ( h(3,k,k)**2 + h(3,k+1,k)**2 )
	  mu(4) = sqrt ( h(4,k,k)**2 + h(4,k+1,k)**2 )
	  c(1,k) = h(1,k,k) / mu(1)
	  c(2,k) = h(2,k,k) / mu(2)
      c(3,k) = h(3,k,k) / mu(3)
      c(4,k) = h(4,k,k) / mu(4)

      s(1,k) = -h(1,k+1,k) / mu(1)
	  s(2,k) = -h(2,k+1,k) / mu(2)
	  s(3,k) = -h(3,k+1,k) / mu(3)
	  s(4,k) = -h(4,k+1,k) / mu(4)

	 !Part17:
      h(1,k,k) = c(1,k) * h(1,k,k) - s(1,k) * h(1,k+1,k)
	  h(2,k,k) = c(2,k) * h(2,k,k) - s(2,k) * h(2,k+1,k)
	  h(3,k,k) = c(3,k) * h(3,k,k) - s(3,k) * h(3,k+1,k)
	  h(4,k,k) = c(4,k) * h(4,k,k) - s(4,k) * h(4,k+1,k)
      h(1:4,k+1,k) = 0.0D+00
	  	  
	 !Part18:
      call mult_givens ( c(1:4,k), s(1:4,k), k, g(1:4,1:k+1) )

	 !Part19:
      rho(1:4) = abs ( g(1:4,k+1)) 
	  
      itr_used = itr_used + 1
	  
	 !Part20:
      if (( rho(1) <= rho_tol(1) .and. rho(1) <= tol_abs ).and.( rho(2) <= rho_tol(2) .and. rho(2) <= tol_abs ).and.( rho(3) <= rho_tol(3) .and. rho(3) <= tol_abs ).and.( rho(4) <= rho_tol(4) .and. rho(4) <= tol_abs )) then
       exit
      end if

    end Do
   	
    k = k_copy - 1
   
   !Part21:
    y(1:4,k+1) = g(1:4,k+1) / h(1:4,k+1,k+1)
	Do j=1,4
       Do i = k, 1, -1
         y(j,i) = ( g(j,i) - dot_product ( h(j,i,i+1:k+1), y(j,i+1:k+1) ) ) / h(j,i,i)
       end Do

	  !Part22:
       Do i = 1, n
         x(j,i) = x(j,i) + dot_product ( v(j,i,1:k+1), y(j,1:k+1) )
	   end Do
	end Do

   !Part23:
    if (( rho(1) <= rho_tol(1) .and. rho(1) <= tol_abs ).and.( rho(2) <= rho_tol(2) .and. rho(2) <= tol_abs ).and.( rho(3) <= rho_tol(3) .and. rho(3) <= tol_abs ).and.( rho(4) <= rho_tol(4) .and. rho(4) <= tol_abs )) then
      exit
    end if

 end Do
!*********************************************************************************************   
 End
!###########################################################################################