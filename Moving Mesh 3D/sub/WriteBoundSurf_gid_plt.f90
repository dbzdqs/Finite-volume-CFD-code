!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:To Write Edge Based Mesh in *.Txt and *.Plt Files                        //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: November, 6, 2015                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Doc ID: MC5F079F1                                                                    //!
!//                                                                                      //!
!// This Program is Available Through the Website: WWW.MarketCode.ir                     //!
!// It May be Copied, Modified, and Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine WriteBoundSurf_gid_plt(Dim,NP,NR,NFR,BC,IDS,X,Y,Z,FaceType)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NP,NR,NFR,BC,IDS,X,Y,Z,FaceType

 Integer::Dim,I,J,JJ,NP,NR,SFace
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:Dim)::X,Y,Z
!*******************************************************************************************
 
 SFace=0
 Do I=1,NR
     
    IF(BC(I)==1) goto 1  
    

     Open(4,File='BoundSurf.Plt')
 
     Write(4,*) ' TITLE = "Title" '  
     Write(4,*) ' VARIABLES  = X , Y , Z'
     Write(4,*) ' ZONE T="Title", N= ', NP , ' , E= ', NFR(I) ,', ET=QUADRILATERAL, F=FEBLOCK'

     Do J=1,NP
	    Write(4,*) X(J)
     End Do
     Do J=1,NP
	    Write(4,*) Y(J)
     End Do
     Do J=1,NP
	    Write(4,*) Z(J)
     End Do
   
    
    Do J=SFace+1,SFace+NFR(I)
	   Do JJ=3,FaceType(J)+2            
          Write(4,'(I15)',Advance='No') IDS(JJ,J)         
       End Do
       IF(FaceType(J)==3) Write(4,'(I15)',Advance='No') IDS(FaceType(J)+2,J)
       Write(4,*)
    End Do
1   SFace=SFace+NFR(I)

 End do
  
 !Close(4)
!*******************************************************************************************
 End
!###########################################################################################

