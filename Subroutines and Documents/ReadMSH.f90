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
!// Developed by: M. Vakili, Computer Science, Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ReadMSH(Dim,DimIDS,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BCType,BCTitle,CellType,FaceType,MeshDim)
 Implicit None
!*********************************************************************************************
 Intent(In )::Dim,DimIDS 
 Intent(Out)::NP,NC,NF,NR,NFR,IDS,X,Y,Z,BCType,BCTitle,CellType,FaceType,MeshDim

 Integer::Dim,NP,NC,NF,NR,MeshDim,IO,IndexCode,RowNum,SpaceCount,RegionNum,FaceCount,&
          CellCount,BCTitleCount,MixFace,I,K,J,L,Count,M,N,CellFlag,DimIDS
 Real(8),Dimension(1:Dim)::X,Y,Z
 Integer,Dimension(1:Dim)::FaceType,CellType 
 Integer,Dimension(1:100)::NFR,BCType
 Integer,Dimension(1:Dim,1:10)::TempIDS
 Integer,Dimension(1:DimIDS,1:Dim)::IDS
 Integer,Dimension(1:5)::NodeIndex,FaceIndex,CellIndex
 
 Character*200 TempStr,CommentStr
 Character*100 IndexChars,TempChars
 Character*100,Dimension(1:100)::BCTitle
 Character Char
 Character,PARAMETER::Space=' '
!*********************************************************************************************
!Part1
 Open(1,File='fluent.msh')

!Part2
 FaceIndex='' ;BCTitle=''  ;

 RowNum=100 ; SpaceCount=0 ; RegionNum=0 ; FaceCount=1 ; CellCount=1 ; BCTitleCount=1
 CellFlag=0 ; CellType=0   ; L=0         ; M=1         ; N=1         ; Z=0

!Part3
 Do While(RowNum>0)

    I=2
    K=1
    J=1
    Count=0
    MixFace=0
    IndexChars=''
    TempChars=''
    CommentStr=''   
    Char=''

   !Part4
    Read(1,'(A200)',IOSTAT = IO) TempStr
    TempStr=Adjustl(TempStr)
 	
   !Part5   
	If(TempStr(1:1)=='(') Then

    !Part6
     Do While(TempStr(I:I)/=Space .And. TempStr(I:I)/='(')
        IndexChars=IndexChars(:J)//TempStr(I:I)
        I=I+1
        J=J+1
     End Do 

    Else

     RowNum=RowNum-1 	
	 Cycle
		
    EndIF
	
   !Part7
    Read(IndexChars,*) IndexCode
    
   !Part8
    Select Case(IndexCode)

   !Part9
    Case(0)
     J=1

    !Part10
     If(Index(TempStr,"Faces of zone")>0 .Or. Index(TempStr,"Interior faces of zone")>0) Then
      print*,'Faces...'
     !Part11
      Do I=3,Len_Trim(TempStr)                	
         If(TempStr(I:I)/=')' .And. TempStr(I:I)/='"') Then        		   
		  CommentStr=CommentStr(:J)//TempStr(I:I)                       
         End If
         
	     J=J+1
 	  End Do

     !Part12
      BCTitle(BCTitleCount)=CommentStr
      BCTitleCount=BCTitleCount+1

     End If

   !Part13
    Case(2)

    !Prt14
     Do While(TempStr(I:I)/=')') 
	          	
   	    If(TempStr(I:I)/=Space) Then
		 TempChars=TempStr(I:I)
         Read(TempChars,'(I1)') MeshDim            		               
        End If 

        I=I+1

     End Do

   !Part15
    Case(10)
        
    !Part16
     Read(1,'(A200)',IOSTAT = IO) TempStr

    !Part17
	 Do While(TempStr(I:I)/=')')

       !Part18
        If(TempStr(I:I)/='(' .And. TempStr(I:I)/=Space) Then
		 TempChars=TempChars(:L)//TempStr(I:I)
         L=L+1

        ElseIf(TempStr(I:I)==Space .And. Len_Trim(TempChars)>0) Then
         Read(TempChars,'(Z20)') NodeIndex(K)
         K=K+1
         L=0
        End If

        I=I+1

     End Do
            
    !Part19
     If( TempStr(I+1:I+1)/='(' )Read(1,*)

    !Part20
     If(MeshDim==2) Then
      Do Count=NodeIndex(2),NodeIndex(3)
	     Read(1,*) X(Count),Y(Count)
      End Do
     Else If(MeshDim==3) Then
      Do Count=NodeIndex(2),NodeIndex(3)
         Read(1,*) X(Count),Y(Count),Z(Count)
      End Do
     End If

    !Part21
	 NP=NodeIndex(3)-NodeIndex(2)+1
                    
   !Part22
    Case(12)

    !Part23
  	 If(CellFlag==1) Then
      RowNum=RowNum-1      	
      Cycle
     End If
            
    !Part24
     Read(1,'(A200)',IOSTAT = IO) TempStr

    !Part25
     Do While(K<=5)
        If(TempStr(I:I)/='(' .And. TempStr(I:I)/=Space .And. TempStr(I:I)/=')') Then
		 TempChars=TempChars(:L)//TempStr(I:I)
         L=L+1

        Else If(TempStr(I:I)==Space .And. Len_Trim(TempChars)>0) Then
         Read(TempChars,'(Z20)') CellIndex(K)
         K=K+1
         L=0
              
        End If
                
        I=I+1
     End Do

    !Part26
     NC=CellIndex(3)-CellIndex(2)+1

    !Part27
     CellFlag=1

    !Part28
     If(CellIndex(5)==0) Then

     !Part29
      Read(1,'(A200)',IOSTAT = IO) TempStr
      TempStr=Adjustl(TempStr)
                
     !Part30
      Do While(Trim(TempStr)/=')' .And. Trim(TempStr)/='))')
                  	
        !Part31
         K=Len_Trim(TempStr)                    
         Do J=1,K
            If(TempStr(J:J)/=Space .And. TempStr(J:J)/='(') Then
             Read(TempStr(J:J),'(I1)') CellType(CellCount)
             CellCount=CellCount+1
            End If                            
         End Do
                    
        !Part32
         Read(1,'(A200)',IOSTAT = IO) TempStr
         TempStr=Adjustl(TempStr)
         
      End Do

    !Part33      
     Else
      CellType=CellIndex(5)
      
     End If

   !Part34
    Case(13)

    !Part35
     Do M=I,Len_Trim(TempStr)
                
       !Part36
        If(TempStr(I:I)/='(' .And. TempStr(I:I)/=')' .And. TempStr(I:I)/=Space) Then
         TempChars=TempChars(:L)//TempStr(I:I)
         I=I+1
         L=L+1
    
        Else If(Len_Trim(TempChars)>0) Then
         Read(TempChars,'(Z20)') FaceIndex(K)
         TempChars=''
         K=K+1
         I=I+1
         L=0
     
        Else If(TempStr(I:I)/=')' .And. I<Len_Trim(TempStr)) Then
         I=I+1
                                            
        End If
                    
     End Do 

    !Part37
     If(FaceIndex(1)==0) Then
      NF=FaceIndex(3)-FaceIndex(2)+1

    !Part38
     Else
      RegionNum=RegionNum+1
      NFR(RegionNum)=FaceIndex(3)-FaceIndex(2)+1
      BCType(RegionNum)=FaceIndex(4)

     End If

    !Part39
     NR=RegionNum
   	
    !Part40
     If(TempStr(I:I)/='(' .And. FaceIndex(1)/=0) Read(1,'(A1)') Char

    !Part41
     If((TempStr(I:I)=='(') .Or. Char=='(')Then
      Do Count=FaceIndex(2),FaceIndex(3)

        !Part42                  	
         Read(1,'(A200)',IOSTAT = IO) TempStr
         M=1

        !Part43
         Do While(Len_Trim(TempStr) >= M-1)

           !Part44
            If(TempStr(M:M)/=Space) Then
             TempChars=TempChars(:N)//TempStr(M:M)
             N=N+1
 
            Else If(TempStr(M:M)==Space .And. Len_Trim(TempChars)>0) Then
             SpaceCount=SpaceCount+1
             Read(TempChars,'(Z20)')TempIDS(FaceCount,SpaceCount) 
             TempChars=''                               
             N=1                               
       
            End If 
            M=M+1
                            
         End Do

        !Part45
         If(FaceIndex(5)/=0) Then
          FaceType(FaceCount)=SpaceCount-2
         Else
          MixFace=1
          FaceType(FaceCount)=TempIDS(FaceCount,1)
         End If

        !Part46
         IDS(1,FaceCount)=TempIDS(FaceCount,SpaceCount-1)
         IDS(2,FaceCount)=TempIDS(FaceCount,SpaceCount)
    
         Do i=1,FaceType(FaceCount)
            IDS(i+2,FaceCount)=TempIDS(FaceCount,i+MixFace)
         End do
     
        !Part47
         FaceCount=FaceCount+1
         SpaceCount=0
         M=1
         N=1
    
      End Do

     End If
          
    End Select  

   !Part48
	RowNum=RowNum-1
    
   
 End Do !End of Main While
   
 Close(1)
!*********************************************************************************************
 End 
!###########################################################################################
