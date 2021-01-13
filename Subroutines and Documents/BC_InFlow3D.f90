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
!// Date: Mar., 10, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine	BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
 Implicit None
!**********************************************************************************************
 Intent (In   )::Dim,NFI1,NFI2,GM,U0,V0,W0,P0,R0,IDS,Wnp1,NX,NY,NZ,DA,P,ALF,Minf
 Intent (Out  )::Wb

 Integer::Dim,I,NFI1,NFI2,ME
 Real(8)::GM,GM1,MLoc,NXX,NYY,NZZ,PE,RE,CE,VE,HTE,QNE,QTE,RIE,SE,UE,WE,UB,VB,WBB,RB,PB,REB,U0,V0,W0,P0,R0,&
     ITT,ITP,ALF,RIB,HTB,EP1,EP2,EP3,C1,C2,MB,TB,Minf,UBB,C
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,P
 Real(8),Dimension(1:6,1:Dim)::WB
!**********************************************************************************************	
!Part 1:
 ITT=(1.+0.2*Minf*Minf)
 ITP= (ITT**3.5)   * P0
 GM1=GM-1
!Part 2:
 DO I=NFI1+1,NFI2

   !Part 3:
    ME  = IDS(1,I)

    NXX = NX(I)/DA(I)
    NYY = NY(I)/DA(I)
    NZZ = NZ(I)/DA(I)

   !Part 4:
    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    WE = WNP1(4,ME)/RE
    PE = P(ME)
    CE = SQRT(ABS(GM*PE/RE))
    
   !Part 5:
    MLoc = SQRT((UE*UE+VE*VE+WE*WE)/CE)
    
   !Part 6:
	IF(MLoc<1.)then
              
    !Part 7: 
     QNE = UE*NXX+VE*NYY+WE*NZZ
            
    !Part 8:
     RIE =-QNE-2.*CE/GM1
     HTE = (PE/RE)*(GM/GM1)+0.5*(UE*UE+VE*VE+WE*WE)
           
    !Part 9:
	 RIB=RIE
	 HTB=HTE
       
    !Part 10:
	 EP1=1+2/GM1
	 EP2=2*RIB
	 EP3=(GM1/2)*(RIB*RIB-2*HTB)

	 C1=(-EP2+SQRT(EP2*EP2-4*EP1*EP3))/(2*EP1)
	 C2=(-EP2-SQRT(EP2*EP2-4*EP1*EP3))/(2*EP1)

	 C=MAX(C1,C2)
        
    !Part 11:
	 UBB= (2*C/GM1)+RIB

     MB=UBB/C

	 PB=ITP*(1+(GM1/2)*MB**2)**(-GM/GM1)
!!!if(I==NFI2)then
!!!    print*,MB
!!!    pause
!!!endif
	 TB=ITT*(PB/ITP)**(GM1/GM)

	 RB=GM*PB/TB

     UB=Cos(ALF)*UBB
	 VB=Sin(ALF)*UBB
	 WBB=0.0
        
   !Part 12:
    ElseIF(MLoc>=1.)then
        
     RB  = R0
     UB  = U0
     VB  = V0
     WBB = W0
     PB  = P0

	END IF
    
   !Part 13:
	REB= PB/GM1 + 0.5*RB*(UB*UB + VB*VB + WBB*WBB)
    
   !Part 14:
	WB(1,I) = RB
    WB(2,I) = RB*UB
    WB(3,I) = RB*VB
    WB(4,I) = RB*WBB
    WB(5,I) = REB
    WB(6,I) = PB

 END do
!**********************************************************************************************
 End
!############################################################################################## 

