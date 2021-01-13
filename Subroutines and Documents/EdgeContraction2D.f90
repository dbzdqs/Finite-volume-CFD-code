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
!// Date: June, 10, 2017                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: A. Hemati zadeh, Mechanical Eng., Amirkabir University of Technology     //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine EdgeContraction2D(Dim,EI,Dead,Heir,NConectEdge,IConectEdge,NEdgeOfCell,IEdgeOfCell,NF,NC,IDS,NR,NFR,BeginOfReg)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,EI,NEdgeOfCell,IEdgeOfCell
 Intent(InOut)::Dead,Heir,NConectEdge,IConectEdge,NF,NC,IDS

 Integer::Dim,J,NF,NC,EI,jj,NR,Region,RegionLastFace
 Integer::ME,NE,N1,N2,N3,N4,EP1,EP2,Heir,Dead,E1,E2,E3,E4,I,P1,P2,Edge,Tmp
 Integer,Dimension(1:Dim)::NEdgeOfCell
 Integer,Dimension(1:4,1:DIM)::IEdgeOfCell
 Integer,Dimension(1:Dim)::NConectEdge
 Integer,Dimension(1:Dim,1:100)::IConectEdge
 Integer,Dimension(1:3)::Edg
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
 LOGICAL::MEIsTriangle=.FALSE.,NEIsTriangle=.FALSE.
 Integer,Dimension(1:100)::NFR,BeginOfReg
!********************************************************************************************* 
!Part 1:
 MEIsTriangle = .FALSE.
 NEIsTriangle = .FALSE.
 
!Part 3:
 ME = IDS(1, EI)
 NE = IDS(2, EI)
 P1 = IDS(3, EI)
 P2 = IDS(4, EI)

!Part 4:
 CALL Find_EdgePointNeib(Dim,Dead,IDS,IEdgeOfCell,EI,EP1,EP2,E1,E2,E3,E4,N1,N2,N3,N4)
 
!Part 5:
 If(IEdgeOfCell(4,ME)==0) Then
     If(IDS(N1,E1)==0 .AND. N2==1)Then
         Tmp=IDS(1,E2)
         IDS(1,E2)=IDS(2,E2)
         IDS(2,E2)=Tmp
         N2=2
     Endif
     
      MEIsTriangle  = .TRUE.
      IDS(N2,E2)    = IDS(N1,E1)
 EndIf
 
!Part 6:
 If(NE/=0)Then     
  If(IEdgeOfCell(4,NE)==0) Then
      If(IDS(N3,E3)==0 .AND. N4==1)Then
         Tmp=IDS(1,E4)
         IDS(1,E4)=IDS(2,E4)
         IDS(2,E4)=Tmp
         N4=2
      Endif
      
      NEIsTriangle = .TRUE.
      IDS(N4,E4)   = IDS(N3,E3)
  Endif
 Endif
 
!Part 7:
 Do I=1,NConectEdge(Dead)
    Edge = IConectEdge(Dead,I)
    IF (IDS(3,Edge)==Dead) IDS(3,Edge)  = Heir
    IF (IDS(4,Edge)==Dead) IDS(4,Edge)  = Heir
 END Do

!Part 8:
 If(ME>NE)Then
     
 !Part 9:
  IF(MEIsTriangle) Call Remove_CELL(Dim,ME,NF,NC,IDS)
  IF(NE/=0 .AND. NEIsTriangle) Call Remove_CELL(Dim,NE,NF,NC,IDS)
  
 Else
     
 !Part 10:
  IF(NE/=0 .AND. NEIsTriangle) Call Remove_CELL(Dim,NE,NF,NC,IDS)
  IF(MEIsTriangle) Call Remove_CELL(Dim,ME,NF,NC,IDS)
  
 EndIf

!Part 11:
 If(.NOT. (MEIsTriangle)) E1 = 0
 If(.NOT. (NEIsTriangle)) E3 = 0

!Part 12:
 Edg(1) = EI
 Edg(2) = E1
 Edg(3) = E3

!Part 13:
 Do I=1,3
    Do J=I+1,3
       IF( Edg(I)<Edg(J) )Then
        Tmp     = Edg(I)
        Edg(I)  = Edg(J)
        Edg(J)  = Tmp
       EndIF
    End Do
 End Do
 
!Part 14:
 Do I=1,3
    IF(Edg(I)/=0)Then
    
     Region = 0
     Do J=1,NR
        If(Edg(I)>=BeginOfReg(J) .AND. (Edg(I)<= (BeginOfReg(J)+NFR(J)-1)))Then
         Region = J
         Exit
        Endif
     EndDo    
     
     RegionLastFace = (BeginOfReg(Region) + NFR(Region) - 1)  
     IDS(:,Edg(I))  = IDS(:,RegionLastFace)
     NFR(Region)    = NFR(Region) - 1
    Endif
 EndDo

!*********************************************************************************************
 End
!###########################################################################################