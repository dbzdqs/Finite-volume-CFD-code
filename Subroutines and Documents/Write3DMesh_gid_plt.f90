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
 Subroutine Write3DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,Z,FaceType)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NP,NF,IDS,X,Y,Z,FaceType

 Integer::Dim,JJ,J,NP,NC,NF
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:Dim)::X,Y,Z
!*******************************************************************************************
!Part 1:
 Open(2,File='3DMesh.Plt')
 
!Part 2: 
 Write(2,*) ' TITLE = "Title" '
 Write(2,*) ' VARIABLES  = X , Y , Z'
 Write(2,*) ' ZONE T="Title", N= ', NP , ' , E= ', NF ,', ET=QUADRILATERAL, F=FEBLOCK'

!Part 3:
 Do J=1,NP
    Write(2,*) X(J)
 End Do
 Do J=1,NP
    Write(2,*) Y(J)
 End Do
 Do J=1,NP
    Write(2,*) Z(J)
 End Do

!Part 4: 
 Do J=1,NF
    Do JJ=3,FaceType(J)+2            
       Write(2,'(I15)',Advance='No') IDS(JJ,J)         
    End Do
    IF(FaceType(J)==3) Write(2,'(I15)',Advance='No') IDS(FaceType(J)+2,J)
    Write(2,*)
 End Do
 
 Close(2)
!*******************************************************************************************
 End
!###########################################################################################

