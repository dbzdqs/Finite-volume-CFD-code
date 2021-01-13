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
!// Chief Developer: M. Namvar, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. A. Jahangirian, Aerospace eng. Amirkabir University of Technology      //!
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: M. Namvar, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: S. Kavoosi, Mechanical Eng., Amirkabir University of Technology          //!
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Solver_DANSTurb_WithSource(&
 Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,&
 NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,Xc,Yc,NX,NY,DA,A,Dw,INW,&
 ERmx,CFLx,NRKS,NWrite,RKJ,&
 Minf,Rinf,MR,ALF,GM,R0,P0,T0,B0,C0,U0,V0,Tt,PrL,PrT,Mu0,Mut0, &
 WB ,WNP1,WTNP1,P,Mu,Mut,DUY,Source_X,Source_Y) 
 Implicit None
!**********************************************************************************************************
 Intent(In   )::Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,&
 NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,Xc,Yc,NX,NY,DA,A,Dw,INW,&
 ERmx,CFLx,NRKS,NWrite,&
 Minf,Rinf,MR,ALF,GM,R0,P0,T0,B0,C0,U0,V0,Tt,PrL,PrT,Mu0,Mut0,RKJ
 
 Intent(InOut)::WNP1,WTNP1,P,Mu,Mut,DUY

 Intent(Out  ):: WB
!**********************************************************************************************************
 Integer::Dim,I,J,NC,NP,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
          NR,NRKS,NWrite,Init,Ncyc,NS  
 Real(8)::GM,Co,ALF,R0,P0,T0,C0,U0,V0,E0,B0,ERmx,CFLx,Minf,Rinf,Tt,MR,PrL,PrT,Rm,U,V,RKco,Mu0,&
          Mut0,Temp ,K2,K4,time
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1,WC,Con,Dif
 Real(8),Dimension(1:Dim)::X,Y,XC,YC,A,NX,NY,DA,DW,DT,P,Mu,Mut,DUX,DUY,DVX,DVY,DTX,DTY,Taukk
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:2,1:Dim)::WTNP1
 Real(8),Dimension(1:5)::RKJ
 
 Real(8),Dimension(1:Dim)::Source_X,Source_Y
!***************************************** Main ********************************************
!Part 1:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

!Part 2:
 Ncyc = 0
 Rm   = 1000.0

!Part 3:
 Do While(Rm > ERmx)

   !Part 4:
    Ncyc=Ncyc+1
     
   !Part 5:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)   
    End Do

   !Part 6:
	Call TimSTP_Turb(Dim,NC,NF,NF1,NF2,IDS,NX,NY,DA,A,CFLx,GM,P,WNP1,WB,Mu,Mut,PrL,PrT,MR,DT)

   !Part 7:
    Do NS=1,NRKS
   
      !Part 8:
	   RKco=RKJ(NS)

      !Part 9:
	   Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)

      !Part 10:
       Call VelTemp_GradFace(Dim,NC,NF1,NF2,NFW1,NFW2,NF,NP,IDS,X,Y,Xc,Yc,WNP1,WB,GM,P,DUX,DUY,DVX,DVY,DTX,DTY)

      !Part 11:
	   Call DifMeanFlow_TurbNoWallFu(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,GM,PrL,PrT,NX,NY,MR,Mu,Mut,&
                                WNP1,WTNP1(1,:),WB,DUX,DUY,DVX,DVY,DTX,DTY,Dif)

      !Part 12:
       Do J=1,NC

		  Co = RKco*DT(J)/A(J)

          WNP1(1,J) = WC(1,J) - Co*( Con(1,J)                        )
		  WNP1(2,J) = WC(2,J) - Co*( Con(2,J)+Dif(2,J) - Source_X(J) ) 
          WNP1(3,J) = WC(3,J) - Co*( Con(3,J)+Dif(3,J) - Source_Y(J) ) 
          WNP1(4,J) = WC(4,J) - Co*( Con(4,J)+Dif(4,J)               )

         !Part 13:
	      U    = WNP1(2,J)/WNP1(1,J)
          V    = WNP1(3,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V) -0.5*WTNP1(1,J) )

         !Part 14:
          Temp = GM*P(J)/WNP1(1,J) 
          Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)         

       End Do

      !Part 15: 
       Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

     End Do !Ns	

    !Part 16:
     Call KWSST_Main(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
                     NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut)

    !Part 17:
 	 Call ResMass(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
    !Print*,Ncyc,Rm

    !Part 18:
     IF( Mod(Ncyc,NWrite)==0 )Then
      !Print*,'Writing Results... ',Ncyc,Rm
	 End If

 End Do !Do While
!*********************************************************************************************
 End 
!###########################################################################################