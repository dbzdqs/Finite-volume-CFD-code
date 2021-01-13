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
!// Developed by: M. H. Saadat, Aerospace Eng., Amirkabir University of Technology         //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeYang_Source(Dim,IDS,NC,A,MR,Wnp1,Wntp1,Mu,Mut,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,&
                          DDEPSX,DDEPSY,Ce1,Ce2,f1,f2,Lk,Le,St)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,IDS,NC,A,MR,Wnp1,Wntp1,Mu,Mut,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDEPSX,DDEPSY,&
                Ce1,Ce2,f1,f2,Lk,Le
 Intent(Out  )::St
 
 Integer::Dim,I,NC
 Real(8)::K,Epsilon,Rho,MR,Ce1,Ce2
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::A,Mu,Mut,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDEPSX,DDEPSY,toxx,toxy,toyy,prok,Proe,Lk,Le,f1,f2
 Real(8),Dimension(1:2,1:Dim)::Wntp1,St
!***************************************************************************************************  
 !Part 1:
 Do I=1,NC
     St(1,I)=0.0
     St(2,I)=0.0
 End Do
 
!part 2:
 Do I=1,NC
     
    !part 3:
    Rho     = Wnp1(1,I)
    k       = Wntp1(1,I)/Rho
    Epsilon = Wntp1(2,I)/Rho
    
    !part 4:
    toxx(I) = (MR)*Mut(I)*( (4.0/3.0)*DDUX(I)-(2.0/3.0)*DDVY(I)  ) -(2.0/3.0)*Rho*K
    toxy(I) = (MR)*Mut(I)*( DDUY(I)+DDVX(I)  )
    toyy(I) = (MR)*Mut(I)*( (4.0/3.0)*DDVY(I)-(2.0/3.0)*DDUX(I)  ) -(2.0/3.0)*Rho*K
    
    !part 5:
    Proe(I) =  Dabs (toxx(I)*DDUX(I) + toxy(I)*(DDUY(I)+DDVX(I)) + toyy(I)*DDVY(I) )
    
    !part 6:
    Prok(I) = min (Proe(I),10.0*Rho*Epsilon)
    
    !part 7:
    St(1,I) = ( -Prok(I) + Rho*Epsilon -Rho*Lk(I) )  
    St(2,I) = ( -Ce1*f1(I)*Epsilon*Proe(I)/k + Ce2*f2(I)*Rho*Epsilon*Epsilon/k - Rho*Le(I) )   
    
    St(1,I) = St(1,I)*A(I)
    St(2,I) = St(2,I)*A(I)
 END DO
!*********************************************************************************************
 End
!###########################################################################################

