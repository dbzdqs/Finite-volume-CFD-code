!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:                                                   //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Doc ID: MC2F000F1                                                                    //!
!//                                                                                      //!
!// This Program is Available Through the Website: www.MarketCode.ir                     //!
!// It May be Copied, Modified, and Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine MeshPreprocForImplicit(Dim,NC,NEdgeOfCell,IEdgeOfCell,INeib,NnonzeroCell,InonzeroCell,InonzeroFace)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NEdgeOfCell,IEdgeOfCell,INeib
 Intent(Inout)::NnonzeroCell,InonzeroCell,InonzeroFace

 Integer::Dim,NC,I,J,II,TEMP,MinInx
 Integer,Dimension(1:6,1:Dim)::INeib
 Integer,Dimension(1:Dim)::NEdgeOfCell
 Integer,Dimension(1:6,1:Dim)::IEdgeOfCell
 Integer,Dimension(1:10,1:Dim)::InonzeroCell
 Integer,Dimension(1:10,1:Dim)::InonzeroFace
 Integer,Dimension(1:Dim)::NnonzeroCell
!*******************************************************************************************

 Do I=1,NC

    NnonzeroCell(I)   = 1
    InonzeroCell(1,I) = I
    InonzeroFace(1,I) = 1

    Do J=1,NEdgeOfCell(I)

!       IF( INeib(J,I)/=0 )Then
            NnonzeroCell(I) = NnonzeroCell(I) + 1
            InonzeroCell(NnonzeroCell(I),I) = INeib(J,I)
            InonzeroFace(NnonzeroCell(I),I) = IEdgeOfCell(J,I)
!       EndIF

    End Do

    DO II=1,NnonzeroCell(I)

       MinInx = InonzeroCell(II,I)
       DO J=II,NnonzeroCell(I)
          MinInx=MIN(InonzeroCell(J,I),MinInx)
       END do

       DO J=II,NnonzeroCell(I)

          IF(MinInx==InonzeroCell(J,I))then
           TEMP=InonzeroCell(II,I)
           InonzeroCell(II,I)=InonzeroCell(J,I)
           InonzeroCell(J,I)=TEMP

           TEMP=InonzeroFace(II,I)
           InonzeroFace(II,I)=InonzeroFace(J,I)
           InonzeroFace(J,I)=TEMP

           Exit
          END if
       END do

    END do

 END DO
!*******************************************************************************************
 End
!###########################################################################################





