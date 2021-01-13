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
!// Developed by: M. Valadkhani, Mechanical Eng., Amirkabir University of Technology       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeChien_Main_DualTim(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                         NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,DT_Real,A,MR,NRKS,RKJ,Wb,WNP1,Mu,DUY,WTNM1,WTN,WTNP1,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,&
                NFF2,NP,IDS,INW,XC,YC,DW,DT,DT_Real,A,MR,NRKS,RKJ,Wb,WNP1,Mu,DUY,WTNM1,WTN
 Intent(Out  )::Mut
 Intent(InOut)::WTNP1

 Integer::Dim,I,II,J,NC,NS,NP,NRKS,ME,NF,NFW1,NFW2,NF1,NF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2
 Real(8)::RKco,K,Epsilon,Rho,MR,Co,tauwall,ustar,yplus,Yn,fmu,Sigk,Sige,Ce1,Ce2,Cmu,DT_Real,AA,Temp  ,Rm_K,Rm_E
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNM1,WTN,WTNP1,WTC,WTB,Cont,Dift,St,Rest
 Real(8),Dimension(1:Dim)::X,Y,NX,NY,XC,YC,DW,A,Mu,Mut,DT,DUY,DUX_C,DUY_C,DVX_C,DVY_C,&
                           DKX_F,DKY_F,DEPSX_F,DEPSY_F
 Real(8),Dimension(1:5)::RKJ
!*********************************************************************************************
!Part 1:
 Sigk=1.0d0
 Sige=1.3
 Ce1=1.35d0
 Ce2=1.8d0
 Cmu=0.09d0

!Part 2:
 Do I=1,NC
    WTC(1,I) =WTNP1(1,I)
    WTC(2,I) =WTNP1(2,I)
 End Do
 
!Part 3:
 Do NS=1,NRKS
      
   !Part 4:    
	RKco=RKJ(NS)
       
   !Part 5:
    Call KeChien_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS,WTB)

   !Part 6:  
    Call Velocity_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C)
            
   !Part 7:  
    Call KFi_FaceGrad(Dim,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,Wnp1,WTNP1,Wb,WTB,DKX_F,DKY_F,DEPSX_F,DEPSY_F)
    Do I=NFS1+1,NFS2
       DKY_F(I)=0.0
	   DEPSY_F(I)=0.0
    End do

   !Part 8:
	Call KFi_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,Wnp1,WTNP1,Wb,WTB,Cont)
    
   !Part 9:
	Call KFi_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DEPSX_F,DEPSY_F,Sigk,Sige,MR,Mu,Mut,Dift)

   !Part 10:
    Call KeChien_Source(Dim,NC,IDS,DW,A,INW,MR,Ce1,Ce2,WNP1,WTNP1,Mu,Mut,DUY,DUX_C,DUY_C,DVX_C,DVY_C,St)
      
   !Part 11:
    Do J=1,NC
        
        AA = A(J)
        Rest(1,J) = -(Cont(1,J)+Dift(1,J)+St(1,J))/AA
        Rest(2,J) = -(Cont(2,J)+Dift(2,J)+St(2,J))/AA
        
        Temp = 2*DT_Real+3*RKco*DT(J)
        
        WTNP1(1,J) = ( 2*WTC(1,J)*DT_Real + RKco*DT(J)*(4*WTN(1,J) - WTNM1(1,J) + 2*Rest(1,J)*DT_Real) ) / Temp
        WTNP1(2,J) = ( 2*WTC(2,J)*DT_Real + RKco*DT(J)*(4*WTN(2,J) - WTNM1(2,J) + 2*Rest(2,J)*DT_Real) ) / Temp
        
      !Part 12: 
       if(WTNP1(1,J)<0.0 )   WTNP1(1,J)=WTC(1,J)    !        WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 )   WTNP1(2,J)=WTC(2,J)    !        WTNP1(2,J)=WTN(2,J)
       
      !Part 13:
       Rho     = WNP1(1,J)
       K       = WTNP1(1,J)/Rho
       Epsilon = WTNP1(2,J)/Rho

      !Part 14:
       II = INW(J)     ! II: Wall Face
       Yn = DW(J) 
       ME = IDS(1,II)

       Tauwall = Mu(ME)*DUY(II)
       Ustar   = Dsqrt(abs(MR*Tauwall/WNP1(1,ME)))
       Yplus   = (1.0/MR)*Rho*Ustar*Yn/Mu(J)

	   Fmu =  1.0 - exp( -0.0115*Yplus )
		       
      !Part 15:  
       Mut(J) = (1.0/MR)*Cmu*Fmu*Rho*k*k / Epsilon 
              
      End Do !J
      
 End Do !NS
 
!*********************************************************************************************
 End
!###########################################################################################
