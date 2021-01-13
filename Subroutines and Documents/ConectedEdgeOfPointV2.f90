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
!// Date: May., 15, 2016                                                                   //!
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
 Subroutine ConectedEdgeOfPointV2(Dim,NF,NP,BeginOfReg,NR,NFR,IDS,NConectEdge,IConectEdge)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF,NP,IDS
 Intent(Out  )::NConectEdge,IConectEdge

 Integer::Dim,I,J,NF,NP,P1,P2,NR
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::NConectEdge
 Integer,Dimension(1:Dim,1:100)::IConectEdge
 Integer,Dimension(1:100)::NFR,BeginOfReg
!*********************************************************************************************
!Part 1:
 Do I=1,NP
    NConectEdge(I) = 0
 End Do 

!Part 2:
 Do J=1,NR
     Do I=BeginOfReg(J),BeginOfReg(J)+NFR(J)-1
        P1  = IDS(3,I)
        P2  = IDS(4,I)
 
       !Part 4:
        NConectEdge(P1) = NConectEdge(P1) + 1
        IConectEdge(P1, NConectEdge(P1) ) = I
 
       !Part 5:
        NConectEdge(P2) = NConectEdge(P2) + 1
        IConectEdge(P2, NConectEdge(P2) ) = I
     EndDo
 EndDo
 
!*********************************************************************************************
 End
!###########################################################################################
