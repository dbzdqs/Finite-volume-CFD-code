!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Calculate the Convection Terms of 2D Mean Flow Equations by Centeral Scheme!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F009F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Multiply(A,B,C)
 Implicit None
!*******************************************************************************************
 Intent (In   )::A,B
 Intent (Out  )::C

 Integer::n,m,k
 Real(8)::X
 Real(8),Dimension(1:5,1:5)::A,B,C
!*******************************************************************************************
 Do n=1,5
    Do m=1,5
       X=0.0
       Do k=1,5
          X=X+A(n,k)*B(k,m)
       end do
       C(n,m)=X
    End Do
 End Do
!*******************************************************************************************
 End
!###########################################################################################
