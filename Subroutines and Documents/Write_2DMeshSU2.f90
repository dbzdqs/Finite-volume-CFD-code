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
Subroutine Write_2DMeshSU2(Dim,MeshDim,NC,Corn,NP,X,Y,NBoundCrv,BCName,NFacCrv,BFacPt,Neib,CellType)
    Implicit None
!*********************************************************************************************
    Intent(In   )                           ::  Dim,MeshDim,NC,Corn,NP,X,Y,NBoundCrv,BCName,NFacCrv,BFacPt,Neib,CellType

    Integer                                 ::  Dim,IO,MeshDim,NC,i,j,k,j1,NP,NBoundCrv,NBE,Dumy
    Character*100                           ::  Temp
    Character*100,Dimension(1:10)           ::  BCName
    Integer,Dimension(1:4,1:Dim)            ::  Neib,Corn
    Integer,Dimension(1:Dim,1:2)            ::  BFacPt
    Integer,Dimension(1:100)                ::  NFacCrv
    Real(8),Dimension(1:Dim)                ::  X,Y
    Integer,Dimension(1:Dim)                ::  CellType
!*********************************************************************************************
    !Part 1: Opening the *.txt file for writting the mesh information
    !        Creating the *.plt file for using tecplot
    Open(3,File='2DMeshSU2.Txt')
    Open(2,File='2DMeshSU2.Plt')
    
    !Part 2: Write the number of Nodes in *.txt file
    Write(3,*) NP  ,' Number of Points'
    
    !Part 3: Write the number of Elements in *.txt file
    Write(3,*) NC  ,' Number of Cells'

    !Part 4:Write the number of Boundary Curves in *.txt file
    Write(3,*) NBoundCrv ,' Number of Boundary Curves' 
    
    !Part 5: Write the number of edges belong to each curve in *.txt file
    Do i=1,NBoundCrv
        write(3,'(10x,I4,1x,A,3x,A)') NFacCrv(i) , ' Number of Edges Belong to Each Curves ', BCName(i)
    End Do
    
    !Part 6: Writting the number of nodes building edges in *.txt file
    NBE=0
    Do j=1,NBoundCrv
        Do j1=NBE+1,NBE+NFacCrv(j)
            Write(3,*) BFacPt(j1,1) , BFacPt(j1,2)
        End Do
	        NBE = NBE + NFacCrv(j)
    End Do
    
    !Part 7: Writting the elements and their neighbours properties in *.txt file 
     Do i=1,NC
        If (CellType(i)==3) Then
            !This is a triangular element and it contains only three nodes
            Write(3,'(3x,I5,2x,I5,2x,I5,2x,I5, 5x,I5)') Corn(1,i),Corn(2,i),Corn(3,i),0 , i
        Else If (CellType(i)==4) Then
            !This is a quadrilateral element and it contains only four nodes
            Write(3,'(3x,I5,2x,I5,2x,I5,2x,I5, 5x,I5)') Corn(1,i),Corn(2,i),Corn(3,i),Corn(4,i) , i
        End If
	    Write(3,'(3x,I5,2x,I5,2x,I5,2x,I5, 5x,I5)') Neib(1,i),Neib(2,i),Neib(3,i),Neib(4,i) , i
     End do
    
     !Part 8: Writting the coordinates of the nodes in *.txt file
     Do i=1,NP
	    Write(3,*) X(i),Y(i)
     End do
     
    !Part 9: Writting the cells properties in *.plt file
    If(NC/=0) Then
        Write(2,*) 'Variables="X","Y"'
        Write(2,*) 'Zone T="Grid"'
        Write(2,*) ' N=  ', NP, ',E= ' , NC, ',F=FEPOINT ET=QUADRILATERAL'
        Do i=1,NP
            Write(2,*) X(i),Y(i)
        End Do
        Do i=1,NC
            If (CellType(i)==3) Then
                Dumy = Corn(3,i)
                If (Dumy==0) Dumy=Corn(2,i)
                Write(2,*) Corn(1,i),Corn(2,i),Corn(3,i),Dumy
            Else If (CellType(i)==4) Then
                Dumy = Corn(4,i)
                If (Dumy==0) Dumy=Corn(3,i)
                Write(2,*) Corn(1,i),Corn(2,i),Corn(3,i),Dumy
            End If
        End Do
    End If
    
    !Part 10: Writting the boundary properties in *.plt file
    NBE=0
    Do i=1,NBoundCrv
    
        Write(2,*) 'Variables="X","Y"'
        Write(2,*) ' ZONE T="', Trim(BCName(I)) , '"'
        Write(2,*) ' N=',NP,' E=',NFacCrv(I),',Datapacking=Point, Zonetype=Fetriangle'
        Do j=1,NP
	        Write(2,*) X(j),Y(j)
        End Do
        Do j = (NBE+1) , (NBE+NFacCrv(i))
  	        Write(2,*) BFacPt(j,1),BFacPt(j,2),BFacPt(j,2)
        End Do

        NBE = NBE + NFacCrv(i)
    End Do
    
    Close(1)
    Close(2)
!*********************************************************************************************    
End Subroutine Write_2DMeshSU2
