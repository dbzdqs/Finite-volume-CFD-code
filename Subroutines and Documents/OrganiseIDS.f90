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
 Subroutine OrganiseIDS(Dim,NBoundCrv,NFacCrv,BFacPt,NF,IDS)
 Implicit None
!*********************************************************************************************
 Intent(In   )                  ::  Dim,NBoundCrv,NFacCrv,BFacPt,NF
 Intent(InOut)                  ::  IDS

 Integer                        ::  Dim,NBoundCrv,NF,i,SF,j,P1,P2,k,FCounter
 Integer,Dimension(1:100)       ::  NFacCrv
 Integer,Dimension(1:Dim,1:2)   ::  BFacPt
 Integer,Dimension(1:4,1:Dim)   ::  IDS,TIDS
 Integer,Dimension(1:Dim)       ::  FCon
!*********************************************************************************************
 !Part 1: Define Temporary IDS
 TIDS(1:4,1:NF) = IDS(1:4,1:NF)
 IDS(1:4,1:NF)  = 0
 FCon(1:NF)     = 0
 
 !Part 2: Cycle all Boundary Edge Curves
 SF=0
 FCounter=0
 Do i=1,NBoundCrv
     Do j=(SF+1),(SF+NFacCrv(i))
         P1=BFacPt(j,1)
         P2=BFacPt(j,2)
         Do k=1,NF
             If ((TIDS(3,k)==P1 .AND. TIDS(4,k)==P2) .OR. &
                (TIDS(3,k)==P2 .AND. TIDS(4,k)==P1)) Then
                 FCounter = FCounter + 1
                 IDS(1:4,FCounter) = TIDS(1:4,k)
                 FCon(k)=1
                 GOTO 100
             End If
         End Do
100 End Do
     SF = SF + NFacCrv(i)
 End Do
 
 !Part 3: Put the non-boundary faces (NE/=0)
 Do i=1,NF
     If (FCon(i)==0) Then
         FCounter = FCounter + 1
         IDS(1:4,FCounter) = TIDS(1:4,i)
     End If
 End Do
 
!*********************************************************************************************
 End Subroutine OrganiseIDS
!###########################################################################################