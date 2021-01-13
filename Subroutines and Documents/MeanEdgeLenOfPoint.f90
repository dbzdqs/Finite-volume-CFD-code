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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine MeanEdgeLenOfPoint(Dim,NF,NP,IDS,X,Y,EdgeLen)
 Implicit None
!********************************************************************************************* 
 Intent(In    ):: Dim,NF,NP,IDS,X,Y
 Intent(Out )::EdgeLen

 Integer::Dim,I,NF,NP,NBP,P1,P2,NE,N
 Real(8)::X1,Y1,X2,Y2,EgeLenght 
 Real(8),Dimension(1:Dim)::X,Y,EdgeLen
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::NConectEdge
!********************************************************************************************* 
!Part 1:
 NConectEdge(:) = 0
 EdgeLen(:) = 0.0
 
!Part 2:
 Do I=1,NF
 
      !Part 3: 
       P1 = IDS(3,I)
       P2 = IDS(4,I)

      !Part 4:       
       X1 = X(P1) ; Y1 = Y(P1)
       X2 = X(P2) ; Y2 = Y(P2)

      !Part 5:       
       EgeLenght = sqrt( (X2-X1)**2 + (Y2-Y1)**2 )
       
      !Part 6:	   
       EdgeLen(P1) = EdgeLen(P1) + EgeLenght  
       EdgeLen(P2) = EdgeLen(P2) + EgeLenght  

      !Part 7:	   
       NConectEdge(P1) = NConectEdge(P1) + 1
       NConectEdge(P2) = NConectEdge(P2) + 1
 End Do

!Part 8:
 Do I=1,NP
     EdgeLen(I) = EdgeLen(I)/NConectEdge(I)
 End Do
!*********************************************************************************************
    End
!###########################################################################################
    
    
    
