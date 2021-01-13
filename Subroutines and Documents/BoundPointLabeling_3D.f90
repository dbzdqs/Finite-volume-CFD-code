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
 Subroutine BoundPointLabeling_3D(Dim,IDS,NR,NFR,BC,FaceType,NBP,IBP)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,IDS,NR,NFR,BC,FaceType
 Intent(Out  )::NBP,IBP

 Integer::DIM,NBP,NR,NBF,I,J,K,M,P1,P2,Face
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType,IBP
 Integer,Dimension(1:100)::NFR,BC
!*********************************************************************************************
!part 1:
 NBF=0
 NBP=0
 IBP=0
 Face=1
 
!part 2: 
 Do I=1,NR
    !part 3:
    !If( BC(I)==1 ) Goto 10
     If( BC(I)==1 .OR. BC(I)==6 ) Goto 10
     
    !part 4: 
     Do J=NBF+1,NBF+NFR(I)
        !part 5:
         Do K=3,FaceType(J)+2
             P1 = IDS(K,J)
             
            !part 6:
             If(Face==1) Then
                 NBP = NBP + 1
                 IBP(NBP) = P1
                 
            !part 7:
             Else If(Face>1) Then
                 Do M=1,NBP
                     P2=IBP(M)
                     If(P1==P2)  GoTo 15
                 End Do
                 NBP = NBP + 1
                 IBP(NBP) = P1
             End If
15       End Do
         
        !part 8:
         Face = Face + 1 
     End Do
     
    !part 9:
10   NBF = NBF + NFR(I)
     
 End Do
!*********************************************************************************************
 End
!###########################################################################################

