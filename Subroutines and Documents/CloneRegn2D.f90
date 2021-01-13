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
 Subroutine CloneRegn2D(Dim,IndxReg,NP,NF,NR,NFR,IDS,X,Y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,IndxReg
 Intent(InOut)::NP,NF,NR,NFR,IDS,X,Y

 Integer::Dim,I,J,Pt,ME,NE,P1,P2,Edge,cnt
 Real(8)::CrosProduct
 
 Integer::NP
 Integer::NF
 Integer::NR
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y

 Integer,Dimension(1:Dim)::NConectEdge
 Integer,Dimension(1:Dim,1:100)::IConectEdge

 Integer::FirstFace,LastFace
 Integer::IndxReg
 
 Integer::NCloneEdg
 Integer::NClonePt
 Integer,Dimension(1:Dim)::ClonePt
!********************************************************************************************* 
!Part 1: 
 FirstFace = 0
 Do I=1,IndxReg-1
    FirstFace = FirstFace + NFR(I)
 End Do
 LastFace = FirstFace + NFR(IndxReg)

!Part 2:
 NCloneEdg = NFR(IndxReg)

!Part 3: 
 NF = NF + NCloneEdg
   
!Part 4:  
 Do i=NF,LastFace+1,-1   
     IDS(:,i) = IDS(:,i-NCloneEdg)
 End Do  

!Part 5:
 Do i=NR,IndxReg,-1  
    NFR(i+1) = NFR(i)
 End Do 

!Part 6:
 NR = NR + 1 

!Part 7: 
 NClonePt = 1
 Pt = IDS(3,FirstFace+1)
 ClonePt(Pt) = NP + NClonePt
     
 X( ClonePt(Pt) ) = X(Pt)
 Y( ClonePt(Pt) ) = Y(Pt)

!Part 8:
 Do i=FirstFace+1,LastFace 

    NClonePt =  NClonePt + 1
    Pt = IDS(4,i)
    ClonePt(Pt) = NP + NClonePt

    X( ClonePt(Pt) ) = X(Pt)
    Y( ClonePt(Pt) ) = Y(Pt)

 End Do  

!Part 9: 
 NP = NP + NClonePt

!Part 10: 
 cnt=LastFace
 Do i=FirstFace+1,LastFace 

    ME = IDS(1,i)
    NE = IDS(2,i)
    P1 = IDS(3,i)
    P2 = IDS(4,i)

    cnt=cnt+1
    IDS(1,cnt) = NE
    IDS(2,cnt) = 0
    IDS(3,cnt) = ClonePt(P2)
    IDS(4,cnt) = ClonePt(P1)

    IDS(2,i) = 0
   
 End Do  
   
!Part 11: 
 Call ConectedEdgeOfPoint(Dim,NF,NP,IDS,NConectEdge,IConectEdge)
 
!Part 12: 
 P1 = IDS(3,LastFace)  
 P2 = IDS(4,LastFace)
  
!Part 13: 
 Do I=1,NConectEdge(P2) 
     
   !Part 14: 
    Edge = IConectEdge(P2,I)     
 
   !Part 15: 
    IF(IDS(3,Edge)==P2) Pt = IDS(4,Edge) 
    IF(IDS(4,Edge)==P2) Pt = IDS(3,Edge) 
           
   !Part 16: 
    CrosProduct = (X(P2)-X(P1))*(Y(Pt)-Y(P1)) - (X(Pt)-X(P1))*(Y(P2)-Y(P1)) 
    if(CrosProduct<0.0)Then
     if(IDS(3,Edge)==P2)   IDS(3,Edge) = ClonePt(P2)
     if(IDS(4,Edge)==P2)   IDS(4,Edge) = ClonePt(P2)
    endif
    
 End Do
    
!Part 17: 
 Do J=FirstFace+1,LastFace 
   
    P1 = IDS(3,J)  
    P2 = IDS(4,J)

    Do I=1,NConectEdge(P1)
        Edge = IConectEdge(P1,I)     
       
        if( Edge>FirstFace+1 .and. Edge<LastFace)cycle
        
        IF(IDS(3,Edge)==P1) Pt = IDS(4,Edge) 
        IF(IDS(4,Edge)==P1) Pt = IDS(3,Edge) 

        CrosProduct = (X(P2)-X(P1))*(Y(Pt)-Y(P1)) - (X(Pt)-X(P1))*(Y(P2)-Y(P1)) 
        if(CrosProduct<0.)Then
         if(IDS(3,Edge)==P1) IDS(3,Edge) = ClonePt(P1)
         if(IDS(4,Edge)==P1) IDS(4,Edge) = ClonePt(P1)
        endif
          
     End Do

 End Do

!*********************************************************************************************
 End
!###########################################################################################