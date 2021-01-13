!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Radidal Basis Function Which has Used for Interpolation                  //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: Agust, 03, 2015                                                                //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Developed by: A.R. Rezaei, Iran, Tehran, a.r.rezaei@aut.ac.ir                        //!
!// Doc ID: MC5F030F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine RBF_Function(DL,Temp)
 Implicit None
!******************************************************************************************* 
 Intent(In  )::DL
 Intent(Out )::Temp

 Real(8)::Temp,DL
!******************************************************************************************* 

!Temp=((1-(DL/5.))**4)*((4*(DL/5.))+1)   !toop

!Temp=Dl
!Temp=(1-(DL/2.))**2                      !wing
 
!Temp=((1-DL/2.)**6)*((35*((Dl/2.)**2))+(18*Dl/2.)+3)
!Temp=((1-DL)**8)*((32*(Dl**3))+(25*(Dl**2))+(8*Dl)+1)
!Part 1:

!IF(DL>1.0e-8)Then    
! Temp = (DL*DL)*(dlog(DL))
!Else 
!Temp = 0.0       
!End If    
 
!Part 2:
!Temp =1./( 1. + (DL/2.)**2)**2              !wing bending-NSTP=40-H=50   **
Temp =1./( 1. + (DL)**2)             !wing bending-NSTP=40-H=50   **

!Part 3:
!Temp = 1. / (1. + (DL/1.)**2)**2            ! NACA0012 Plunging 

!Temp = 1. / (1. + (DL/6.)**2)**2           !   beam bending -NSTP=30-H=100

!Temp = 1. / (1. + (DL/4.)**2)**2           !   ball plunging -NSTP=40-H=10  **

!Temp = 1. / (1. + (DL/5.)**2)**2           !   wing bending -NSTP=30-H=50

!Temp = 1. / (1. + (DL/2.)**2)**2           !   wing plunging1 -NSTP=40-H=15  **
!Temp = 1. / (1. + (DL/4.)**2)**2           !   wing plunging2 -NSTP=40-H=15  **

!Part 4:
!Temp = exp( -DL*DL )                           !   wing bending -NSTP=30-H=50
!*******************************************************************************************
 End
!###########################################################################################
 