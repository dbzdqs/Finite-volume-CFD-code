!DDDDDDDDDDDDDDDDDDDDDDDDDDDAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNSSSSSSSSSSSSSSSSSSSSSSSSSS
!//             /////////////       ////////////////    ////////     //////    ////////////////        //!
!//             /////////////      ////////////////    //////////   //////    ////////////////         //!
!//            /////    /////     //////    //////    //////////// //////    /////                     //!
!//           /////    //////    ////////////////    ///////////////////    ////////////////           //!
!//          /////    //////    ////////////////    ////// ////////////               /////            //!
!//         ///////////////    //////    //////    //////   //////////    ////////////////             //!
!//       ///////////////     //////    //////    //////     ////////    ////////////////              //!
!//          Developer            Assistant    in      Numerical             Sciences                  //!
!//----------------------------------------------------------------------------------------------------//!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace Department, Amirkabir University of Technology           //!
!// Chief Developer: N. msnkre, Aerospace eng., Amirkabir University of Technology                     //!
!// Date: October, 14, 2013                                                                            //!
!//                                                                                                    //!
!// The Program is Available Through the Website: www.DANS.ir                                          //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                               //!
!//----------------------------------------------------------------------------------------------------//!
!// Description:                                                                                       //!
!// This code are designed to implement the PSO algorithm. However this code is used to solve an 
!// academic test case but it could be used for any optimization problem by modifying objective function/!                                                                                                  //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
Program DANS_OptAlgorithmPSO
 Implicit None
!===============================
 Integer,Parameter::Dim=10000  
 Integer,Parameter::SwarmSize=5 
!===============================
 
 Integer::I,J,S,Iter,Check
 Integer::ArrayLength
 Integer::Maxiter
 Real(8)::ConvergERROR 
 Real(8),Allocatable,Dimension(:)::InitialArray
 Real(8),Allocatable,Dimension(:)::gBest 
 Real(8),Allocatable,Dimension(:)::Pmax
 Real(8),Allocatable,Dimension(:)::Pmin
 Real(8),Allocatable,Dimension(:)::ObjPbest
 Real(8),Allocatable,Dimension(:)::ObjFuncS                                     
 Real(8),Allocatable,Dimension(:,:)::Velocity                                    
 Real(8),Allocatable,Dimension(:,:)::Pos          
 Real(8),Allocatable,Dimension(:,:)::pBest           
 Real(8),Allocatable,Dimension(:,:)::Delta             
 Real(8),Allocatable,Dimension(:,:)::TempMatrix                   
 Real(8),Allocatable,Dimension(:,:)::CheckMatrix
 Real(8)::ObjGbest
 Real(8)::ObjectiveFun
 Logical::Converged=.False.
!*********************************************************************************************
!Part 1:
 ArrayLength = 2  !x,y
 Maxiter = 10000
 ConvergERROR= 1e-6

!Part 2:
 Allocate( InitialArray(1:ArrayLength)             )
 Allocate( Pos         (1:SwarmSize,1:ArrayLength) )
 Allocate( Velocity    (1:SwarmSize,1:ArrayLength) )
 Allocate( pBest       (1:SwarmSize,1:ArrayLength) )
 Allocate( gBest       (1:ArrayLength)             )
 Allocate( Pmax        (1:ArrayLength)             )
 Allocate( Pmin        (1:ArrayLength)             )
 Allocate( ObjpBest    (1:SwarmSize)               )
 Allocate( CheckMatrix (1:SwarmSize,1:ArrayLength) )
 Allocate( ObjFuncS    (1:SwarmSize)               )

!Part 3:
 Do I=1,ArrayLength
    InitialArray(I) = 1.0
 End Do 

!Part 4:
 Call PSO_Initial(SwarmSize,ArrayLength,InitialArray,Pmax,Pmin,Pos,Velocity,pBest)
      
!Part 5:
 Open(15,File='Check.dat')
 Open(55,File='Converge Check.dat')
 Iter = 0
 Do While ( (Iter <= Maxiter).And.(.Not.Converged) )

	Iter = Iter+1
    Print*,'Iteration:',Iter 

	Do S=1,SwarmSize
        
      !Part 6:
       ObjectiveFun = 	(Pos(S,1) + 2* Pos(S,2) - 7)**2 + (2*Pos(S,1) + Pos(S,2) - 5)**2

       write(15,'(a,I4,I4,F10.6)') 'Iter - S - ObjectiveFun' , Iter , S , ObjectiveFun
       write(* ,'(a,I4,I4,F10.6)') 'Iter - S - ObjectiveFun' , Iter , S , ObjectiveFun

       ObjFuncS(S) =  ObjectiveFun

      !Part 7:
       If(Iter==1)  ObjPbest(S) = ObjectiveFun

      !For Maximizing >
	  !for minimizing <
       If( ObjectiveFun < ObjPbest(S) )Then
        Do I=1,ArrayLength
           pBest(S,I) = Pos(S,I)
        End Do
        ObjPbest(S) = ObjectiveFun
       End If
	
       write(15,'(a,I4,F10.6)') 'S - ObjPbest' , S , ObjPbest(S)
       write(* ,'(a,I4,F10.6)') 'S - ObjPbest' , S , ObjPbest(S)

      !Part 8: 
       If((Iter==1).And.(S==1)) Then
        ObjGbest =  ObjPbest(S)
        Do I=1,ArrayLength
	       gBest(I) = pBest(S,I)
        End Do
       End If

      !For Maximizing >
      !for minimizing <
       If( ObjPbest(S) < ObjGbest )Then
        Do I = 1,ArrayLength
           gBest(I) = pBest(S,I)
	    End Do
        ObjGbest = ObjPbest(S)
       End If

       write(15,'(a,I4,F10.6)') ' Iter - ObjGbest' , Iter , ObjGbest
       write(* ,'(a,I4,F10.6)') ' Iter - ObjGbest' , Iter , ObjGbest

	  !Part 9:
	   Call PSO_Update(Iter,Maxiter,SwarmSize,ArrayLength,Pmax,Pmin,pBest,gBest,ObjFuncS,ObjGbest,Pos,Velocity)
   
    End Do
   
   !Part 10:
    Check = 0
    Do S=1,SwarmSize
       Do I = 1,ArrayLength
          If( ABS(Pos(S,I) - gBest(I)) <= ConvergERROR ) Check = Check + 1
       End Do
    End Do
  
    If(Check == (SwarmSize* ArrayLength) ) converged =.True.
			
   !Part 11:
    write(55,'(I4,F10.6,F10.6)') Iter , ObjFuncS(1) , ObjFuncS(2)
 
 End Do

!Part 12:
 write(15,*) '===================================='
 write(15,'(a)') 'gBest'
 Do I = 1,ArrayLength
    write(15,'(F10.6)')  gBest(I)
 End Do
 	
 Write(15,'(a)') 'Final ObjGbest' 
 write(15,'(F10.6)') ObjGbest
!*********************************************************************************************
 End Program DANS_OptAlgorithmPSO
!###########################################################################################
 