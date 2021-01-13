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
!// Developed by: A. Moslemi Pak, Mechanical Eng., Amirkabir University of Technology      //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine DualMeshGenerator2D(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,XCen,YCen,NewNP,NewNC,NewNF,NewNR,NewBC,NewNFR,NewIDS,NewX,NewY)
Implicit None
!*********************************************************************************************
Intent(In   )                       ::  Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,XCen,YCen
Intent(Out  )                       ::  NewNP,NewNC,NewNF,NewNR,NewBC,NewNFR,NewIDS,NewX,NewY

Integer                             ::  Dim,i,j,NC,NR,NF,NP,NewNP,NewNC,NewNF,NewNR,SF,MBPNo
Integer                             ::  P1,P2,ME,NE,BP1,BP2,BF1,BF2,IntF,IntBCNO,NewPoint
Integer,Dimension(1:4,1:Dim)        ::  IDS,NewIDS
Real(8),Dimension(1:Dim)            ::  X,XCen,Y,YCen,NewX,NewY
Integer,Dimension(1:100)            ::  NFR,BC,NewBC,NewNFR
Integer,Dimension(1:Dim)            ::  MBP,CheckNewP
!*********************************************************************************************
!Part 1: Predefining the Voronoi Diagram variables
NewNP       = NC
NewNC       = NP
NewNF       = 0
NewIDS      = 0
NewBC(1:NR) = 0
NewNR       = 1     !"NewNR" counts the number of BC for V.Diagram
CheckNewP   = 0     !"CheckNewP" prevents to produce new nodes on the boundaries
    
!Part 2: Create only the new interior faces and considering their region and their ME and NE
SF  = 0
MBPNo = 0           !"MBPNo" is defined for considering the number of middle boundary point creation
NewBC(1) = 1        !The V.Diagram 1st BC must be interior
               
Do i=1,NR
    If (BC(i)==1) Then
        !Part 3: This condition creates the interior faces from the interior faces of the input mesh
        Do j=SF+1,SF+NFR(i)
            !One new face adds
            NewNF       = NewNF     + 1
            NewNFR(1)   = NewNFR(1) + 1
            !First define the local variables of the current new face using the 2D input Mesh IDS
            P1 = IDS(2,j)
            P2 = IDS(1,j)
            ME = IDS(3,j)
            NE = IDS(4,j)
            !Create the New IDS for the new added face Using the above local variables
            !Here "NewIDS" defines a temporary variable for NewIDS's interior faces
            NewIDS(1,NewNF) = ME
            NewIDS(2,NewNF) = NE
            NewIDS(3,NewNF) = P1
            NewIDS(4,NewNF) = P2
            !Define the Coordinates of "P1" and "P2" using their correspond Cell centroids
            NewX(P1) = XCen(P1)
            NewY(P1) = YCen(P1)
            NewX(P2) = XCen(P2)
            NewY(P2) = YCen(P2)
        End Do
        !Add the total faces left to its previous cumulation
        SF=SF+NFR(i)
    Else
        !Part 4: This condition creates the interior faces using the boundary faces of the 2D input Mesh
        Do j=SF+1,SF+NFR(i)
            !One face creates by connecting the centroid of the boundary cell to center of its boundary edge
            NewNF           = NewNF + 1                 !Total face number has increased one number
            NewNFR(1)       = NewNFR(1) + 1             !NFR for interior faces adds one more
                
            !One new points add to the rest of the Voronoi Diagram points. (The middle of boundary edge)
            NewNP = NewNP + 1
            MBPNo = MBPNo + 1
            MBP(MBPNo) = NewNP                          !"MBP" is the middle of the boundary points BP1 and BP2
                
            !First define the local variables of the current new face using the 2D input Mesh IDS
            P1 = NewNP
            P2 = IDS(1,j)
            ME = IDS(3,j)
            NE = IDS(4,j)
                
            !Define the interior face (IntF)
            NewIDS(1,NewNF) = ME
            NewIDS(2,NewNF) = NE
            NewIDS(3,NewNF) = P1
            NewIDS(4,NewNF) = P2
                
            !Define the Coordinates
            NewX(NewNP) = 0.50000*(X(NE)+X(ME))
            NewY(NewNP) = 0.50000*(Y(NE)+Y(ME))
                
        End Do
        !Add the total faces left to its previous cumulation
        SF=SF+NFR(i)
    End If
End Do
    
!Part 5: Build the not interior faces from the not interior faces of the 2D input Mesh
SF  = 0
MBPNo = 0
Do i=1,NR
    If (BC(i)==1) Then
        !Skip for only not interior faces of the 2D input Mesh
        SF=SF+NFR(i)
    Else
        !Creating the not interior faces of the V.Diagram
        NewNR           = NewNR + 1
        NewBC(NewNR)    = BC(i)
        Do j=SF+1,SF+NFR(i)
            !Two faces in the current region of the V.Diagram creates
            NewNF           = NewNF + 2                 !Total face number has increased two number
            NewNFR(NewNR)   = NewNFR(NewNR) + 2       !NFR for boundary faces adds two more (current region)
            BF1             = NewNF - 1                 !"BF1" is the number of 1st Boundary Faces
            BF2             = NewNF                     !"BF2" is the number of 2nd Boundary Faces
                
            !According to the "CheckNewP" may add one or two or none new nodes
            !"BP1" and "BP2" are the 1st and 2nd boundary points,
            !which they are getting from the input mesh boundary points.
                
            !Check whether the "BP1" and "BP2" has been considered before or not
            If (CheckNewP(IDS(3,j))==0) Then
                !The "BP1" has not been considered yet. So the total points will increase one.
                NewNP   = NewNP + 1
                BP1     = NewNP
            Else
                !The "BP1" has been considered. None new point is created.
                BP1     = CheckNewP(IDS(3,j))
            End If
            If (CheckNewP(IDS(4,j))==0) Then
                !The "BP2" has not been considered yet. So the total points will increase one.
                NewNP   = NewNP + 1
                BP2     = NewNP
            Else
                !The "BP2" has been considered. None new point is created.
                BP2     = CheckNewP(IDS(4,j))
            End If
            
            MBPNo   = MBPNo + 1
            NewPoint = MBP(MBPNo)                       !Use the middle boundary point creates in the previous loop
            !Define the boundary faces (BF1 and BF2)
            NewIDS(1,BF1)       = IDS(3,j)
            NewIDS(2,BF1)       = 0
            NewIDS(3,BF1)       = BP1
            NewIDS(4,BF1)       = NewPoint
            CheckNewP(IDS(3,j)) = BP1
                
            NewIDS(1,BF2)       = IDS(4,j)
            NewIDS(2,BF2)       = 0
            NewIDS(3,BF2)       = NewPoint
            NewIDS(4,BF2)       = BP2
            CheckNewP(IDS(4,j)) = BP2
              
            !Define the Coordinates
            NewX(BP1) = X(IDS(3,j))
            NewY(BP1) = Y(IDS(3,j))
            NewX(BP2) = X(IDS(4,j))
            NewY(BP2) = Y(IDS(4,j))
                
        End Do
        !Add the total faces left to its previous cumulation
        SF=SF+NFR(i)
    End If
        
End Do
    
!*********************************************************************************************
 End Subroutine DualMeshGenerator2D
!###########################################################################################