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
!// Developed by: M. Valadkhani, Mechanical Eng., Amirkabir University of Technology       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine PressLiftDragCo_Viscous(Dim,NFW1,NFW2,GM,Rinf,Minf,ALF,IDS,WNP1,WB,NX,NY,DA,Mu,DUY,CL,CD)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,DA,IDS,Rinf,WNP1,Mu,DUY
 Intent(Out  )::CL,CD

 Integer::Dim,I,NFW1,NFW2,ME
 Real(8)::CPP,CFC,GM,Rinf,Minf,ALF,Cnx,Cny,DAA,NXX,NYY,CP,CL,CD,ALFa,PI,U,V,VT,Direct,CF
 Real(8),Dimension(1:Dim)::NX,NY,DA ,Mu,DUY
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
!*********************************************************************************************
!Part 1:
 PI=4.*Atan(1.0)
 
!Part 2:
 ALFa = ALF*(PI/180.)
 
!Part 3:
 CPP = 0.5*GM*Minf*Minf
 
!Part 4:
 CFC = 2/(Rinf*Minf)
 
!Part 5:
 Cnx = 0.0
 Cny = 0.0
 
!Part 6:
 DO I=NFW1+1,NFW2
     
    !Part 7:
     DAA = 1.0 !DA(I)
     NXX = NX(I)/DAA
     NYY = NY(I)/DAA
     
    !Part 8:
     ME= IDS(1,I)
     U = WNP1(2,ME)/WNP1(1,ME)
     V = WNP1(3,ME)/WNP1(1,ME)
     
    !Part 9:
     VT = -NYY*U + NXX*V
     Direct = VT / Dabs(VT)
     
    !Part 10:
     CP = (GM*WB(5,I)-1.0) / CPP
     
    !Part 11:
     CF = CFC * Mu(ME) * DUY(I)
     
    !Part 12:
     Cnx = Cnx + NXX*CP - NYY*CF*Direct
     Cny = Cny + NYY*CP + NXX*CF*Direct
     
 END DO
 
!Part 13:
 CD = Cnx*DCOS(ALFa) + Cny*DSIN(ALFa)
 CL = Cny*DCOS(ALFa) - Cnx*DSIN(ALFa)
 
!For 2M5,2M6 test Cases
!CD = Cnx
!CL = Cny
 
!*********************************************************************************************
 End
!###########################################################################################
