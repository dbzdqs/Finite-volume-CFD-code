!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:                                                    //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Developed by: E. FarhadKhani, Iran, Tehran, Efarhadkhani@gmail.com                   //!
!// Doc ID: MC5F000F1                                                                    //!
!//                                                                                      //!
!// This Program is Available Through the Website: www.MarketCode.ir                     //!
!// It May be Copied, Modified, and Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine NeibOfCell3D(Dim,NC,IDS,NFace_Cell,IFace_Cell,INeib)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,IDS,NFace_Cell,IFace_Cell
 Intent(Out  )::INeib
 
 Integer::Dim,NF,NC,ME,NE,I,J,face
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:6,1:Dim)::INeib
 Integer,Dimension(1:6,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::NFace_Cell
!*******************************************************************************************
 INeib(:,:) = 0
 
 Do I=1,NC
    Do J=1,NFace_Cell(I)
        
       Face = IFace_Cell(J,I)
    
       ME = IDS(1,Face)
       NE = IDS(2,Face)  

       IF( ME==I ) INeib(J,ME) = NE

    End Do 
 End Do 
!*******************************************************************************************
 End
!###########################################################################################





