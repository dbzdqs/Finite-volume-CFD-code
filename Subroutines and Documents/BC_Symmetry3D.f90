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
!// Date: Mar., 05, 2013                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine	 BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Implicit None
!**********************************************************************************************
 Intent (In   )::Dim,NFS1,NFS2,GM,U0,V0,P0,R0,C0,IDS,Wnp1,NX,NY,DA,P
 Intent (Out  )::Wb

 Integer::Dim,I,NFS1,NFS2,ME,P1,P2
 Real(8)::GM,GM1,U,V,CC,MLoc,C0,U0B,V0B,P0B,R0B,C0B,S0,DX,DY,DH,STH,CTH,QN0,QT0,RI0,PE,RE,CE,&
          VE,QNE,QTE,QNB,QTB,RIE,SE,UE,QNN,C,QTT,SB,UB,VB,RB,PB,REB,U0,V0,P0,R0,ITT,ITP,AII,HTE,RIB,HTB,EP1,&
		  EP2,EP3,C1,C2,MB,TB,TE,WE,REE
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:6,1:Dim)::Wb
!**********************************************************************************************	
 GM1= GM-1
 DO I=NFS1+1,NFS2

    ME  = IDS(1,I)


    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    WE = WNP1(4,ME)/RE
    PE = P(ME)

	REE= PE/GM1 + 0.5*RE*(UE*UE + VE*VE + WE*WE)

	Wb(1,I) = RE
    Wb(2,I) = RE*UE
    Wb(3,I) = RE*VE
    Wb(4,I) = RE*WE
    Wb(5,I) = REE
    Wb(6,I) = PE

	END do

!**********************************************************************************************
 End
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS