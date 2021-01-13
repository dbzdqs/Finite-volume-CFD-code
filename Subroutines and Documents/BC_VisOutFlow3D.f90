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
 Subroutine BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)
 Implicit None
!**********************************************************************************************
 Intent(In   )::Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P
 Intent(Out  )::WB

 Integer::Dim,I,ME,NFO1,NFO2
 Real(8)::GM,GM1,NXX,NYY,NZZ,U,V,W,CC,MLoc,PE,RE,UE,VE,WE,TE,RB,UB,VB,WBB,PB,TB,REB,P0
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,P
 Real(8),Dimension(1:6,1:Dim)::WB
!**********************************************************************************************	
 GM1= GM-1

!Part 1:
 DO I=NFO1+1,NFO2

   !Part 2:
    ME  = IDS(1,I)

    U = WNP1(2,ME)/WNP1(1,ME)
    V = WNP1(3,ME)/WNP1(1,ME)
    W = WNP1(4,ME)/WNP1(1,ME)

    CC = GM*P(ME)/WNP1(1,ME)
	MLoc = SQRT((U*U+V*V+W*W)/CC)
	
   !Part3 3:	
    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    WE = WNP1(4,ME)/RE
    PE = P(ME)
	TE = PE*GM/RE
	
   !Part 4:
	IF(MLoc<1.)then
	PB=P0
	UB=UE
	VB=VE
	WBB=WE
	TB=TE
	RB=PB*GM/TB
	
   !Part 5:
	Elseif(MLoc>1.)then
	PB=PE
	UB=UE
	VB=VE
	WBB=WE
	TB=TE
	RB=PB*GM/TB

	END IF
	
   !Part 6:
	REB= PB/GM1 + 0.5*RB*(UB*UB + VB*VB + WBB*WBB)

   !Part 7:    
	Wb(1,I) = RB
    Wb(2,I) = RB*UB
    Wb(3,I) = RB*VB
    Wb(4,I) = RB*WBB
    Wb(5,I) = REB
    Wb(6,I) = PB

 END do
!**********************************************************************************************
 End
!##############################################################################################