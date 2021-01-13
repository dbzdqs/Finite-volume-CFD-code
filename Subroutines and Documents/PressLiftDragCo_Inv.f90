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
!// Developed by: *//*-+/                       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine PressLiftDragCo(Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,DA,CL,CD)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,DA
 Intent(Out  )::CL,CD

 Integer::Dim,I,NFW1,NFW2
 Real(8)::CPP,GM,Minf,ALF,Cnx,Cny,DAA,NXX,NYY,CP,CL,CD,ALFa,PI
 Real(8),Dimension(1:Dim)::NX,NY,DA
 Real(8),Dimension(1:5,1:Dim)::WB
!*********************************************************************************************
!Part 1:
 PI=4.*Atan(1.0)
 
!Part 2:
 ALFa = ALF*(PI/180.)
 
!Part 3:
 CPP = 0.5*GM*Minf*Minf
 
!Part 4:
 Cnx = 0.0
 Cny = 0.0
 
!Part 5:
 DO I=NFW1+1,NFW2
     
    !Part 6:
     DAA = 1.0 !DA(I)
     NXX = NX(I)/DAA
     NYY = NY(I)/DAA
     
    !Part 7:
     CP = (GM*WB(5,I)-1.0)/CPP
     
    !Part 8:
     Cnx = Cnx + NXX*CP
     Cny = Cny + NYY*CP
    
 END DO
 
!Part 9:
 CD=Cnx*DCOS(ALFa)+Cny*DSIN(ALFa)
 CL=Cny*DCOS(ALFa)-Cnx*DSIN(ALFa)
 
!For 2M5,2M6 test Cases
!CD = Cnx
!CL = Cny
 
!*********************************************************************************************
 End
!###########################################################################################
