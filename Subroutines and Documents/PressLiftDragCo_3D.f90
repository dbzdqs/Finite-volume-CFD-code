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
 Subroutine PressLiftDragCo_3D(Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,NZ,DA,Unsteady_Moving,PlanForm_Area,CL,CD)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,NZ,DA,Unsteady_Moving
 Intent(Out  )::CL,CD
 Intent(InOut)::PlanForm_Area

 Integer::Dim,I,NFW1,NFW2,Unsteady_Moving
 Real(8)::CPP,GM,Minf,ALF,Cnx,Cny,Cnz,DAA,NXX,NYY,NZZ,CP,CL,CD,ALFa,PI  ,PlanForm_Area
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA
 Real(8),Dimension(1:6,1:Dim)::WB
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
!Cnz = 0.0
  
!Part 5:
 DO I=NFW1+1,NFW2
     
    !Part 6:
     DAA = 1.0 !DA(I)
     NXX = NX(I)/DAA
     NYY = NY(I)/DAA
 !   NZZ = NZ(I)/DAA
     
    !Part 7:
     CP = (GM*WB(6,I)-1.0)/CPP
     
    !Part 8:
     !Cnx = Cnx + NXX*CP
     !Cny = Cny + NYY*CP
    !Cnz = Cnz + NZZ*CP
     
     Cnx = Cnx + NXX !*WB(6,I)
     Cny = Cny + NYY !*WB(6,I)
     
 END DO
 
!Part 9:
 If(Unsteady_Moving == 0) Then
     
    !Part 10:
     PlanForm_Area=0.0
     
    !Part 11:
     DO I=NFW1+1,NFW2
         
        !Airfoil in X direction
         PlanForm_Area = PlanForm_Area + DABS(NY(I))
         
        !Airfoil in Y direction
        !PlanForm_Area = PlanForm_Area + DABS(NX(I))
        
     End Do
     
    !Part 12:
     PlanForm_Area = PlanForm_Area/2.
     
 End If
 
!Part 13:
 Cnx = Cnx !/ PlanForm_Area
 Cny = Cny !/ PlanForm_Area
 
!Part 14:
 CL=Cny*DCOS(ALFa)-Cnx*DSIN(ALFa)
 CD=Cnx*DCOS(ALFa)+Cny*DSIN(ALFa)
 
!*********************************************************************************************
 End
!###########################################################################################
