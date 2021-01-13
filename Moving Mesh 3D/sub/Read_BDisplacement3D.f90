!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate Displacement of non-Boundary Points by Radfial Basis Function //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: Agust, 03, 2015                                                                //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Developed by: A.R. Rezaei, Iran, Tehran, a.r.rezaei@aut.ac.ir                        //!
!// Doc ID: MC5F033F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Read_BDisplacement3D(Dim,NStp,NBP,IBP,DelX,DelY,DelZ)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim
 Intent(Out  )::NStp,NBP,IBP,DelX,DelY,DelZ

 Integer::Dim,I,NBP,NStp
 Integer,Dimension(1:dim)::IBP
 Real(8),Dimension(1:dim)::DelX,DelY,DelZ
!*******************************************************************************************
!Part 1:
 Open(11,File='BDisplacement.Txt')

!Part 2:
 Read(11,*) NStp

!Part 3:
 Read(11,*)
 Read(11,*) NBP 

!Part 4:	  
 Do I=1,NBP
    Read(11,'(I15)',Advance='No') IBP(I) 
    Read(11,*)   DelX( IBP(I) ),DelY( IBP(I) ),DelZ( IBP(I) )
  
 End Do

!*******************************************************************************************
 End
!###########################################################################################

