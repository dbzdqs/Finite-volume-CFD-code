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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KFi_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,WNP1,WTNP1,WB,WTB,Cont)
 Implicit None
!**********************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,WNP1,WTNP1,WB,WTB
 Intent(Out  )::Cont

 Integer::Dim,I,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE,P1,P2
 Real(8)::U,V,DX,DY,F1,F2,Q,R,RU,RV,RK,RFi
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB,Cont
 Integer,Dimension(1:4,1:Dim)::IDS
!**********************************************************************************************	
!Part 1:
 DO I=1,NC
    Cont(1,I) = 0.0
    Cont(2,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF
 
   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I)+V*NY(I)

   !Part 6:
    F1 = Q*WTB(1,I)
    F2 = Q*WTB(2,I)

   !Part 7:
    Cont(1,ME) = Cont(1,ME) + F1
    Cont(2,ME) = Cont(2,ME) + F2

 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 10:  
    U = 0.5*( WNP1(2,ME)/WNP1(1,ME) + WNP1(2,NE)/WNP1(1,NE) )
    V = 0.5*( WNP1(3,ME)/WNP1(1,ME) + WNP1(3,NE)/WNP1(1,NE) )

   !Part 11:
    Q  = U*NX(I)+V*NY(I)
    
    !Part 12:
    if( Q>=0.0 )then
     RK  = WTNP1(1,ME)
     RFi = WTNP1(2,ME)
    else
     RK  = WTNP1(1,NE)
     RFi = WTNP1(2,NE)
    end if
    
    F1 = Q * RK
    F2 = Q * RFi

   !Part 13:
    Cont(1,ME) = Cont(1,ME) + F1
    Cont(2,ME) = Cont(2,ME) + F2

    Cont(1,NE) = Cont(1,NE) - F1
    Cont(2,NE) = Cont(2,NE) - F2

 End Do
!*********************************************************************************************
 End
!###########################################################################################
