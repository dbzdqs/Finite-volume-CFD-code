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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine CG_Matrix_Solver(Dim,NBP,A,b,x)
 Implicit None
!********************************************************************************************* 
 Intent(In   )::Dim,NBP,A,b
 Intent(Out  )::x

 Integer::Dim,NBP,I,J,Iter
 Real(8)::norm0,norm,Sumation,Num,Denom,alpha,beta
 Real(8),Dimension(1:NBP+3)::x,b,p,r
 Real(8),Dimension(1:NBP+3,1:NBP+3)::A
!********************************************************************************************* 
!Part :
 Do I=1,NBP+3
    x(I)=0.0
 End Do

!Part :
 norm0=0

 Do I=1,NBP+3
   
    Sumation=0   
   
    Do J=1,NBP+3
       Sumation=Sumation+(A(I,J)*x(J))
    End Do
   
    r(I)=b(I)-Sumation
    p(I)=r(I)
    norm0=norm0+r(I)*r(I)
   
 End Do

 norm=norm0

!Part :
 Iter=0
 if(norm0<1e-11) iter=30000

!Part :
 Do  while ((Iter<900))

    Iter=Iter+1
    Num=0
    Denom=0
   
   !Part :
    Do I=1,NBP+3
       Num=Num+(r(i)*r(i))
       Sumation=0              
       Do J=1,NBP+3,1
          Sumation=(A(I,J)*p(J))+Sumation
       End Do
       Denom=Denom+(Sumation*p(I))
    End Do
       
    alpha=Num/Denom
   
   !Part :
    Do I=1,NBP+3
       x(I)=x(I)+(alpha*p(I))
    End Do
   
   !Part :
    Do I=1,NBP+3
       Sumation=0              
       Do J=1,(NBP+3),1
          Sumation=(A(I,J)*x(J))+Sumation
       End Do
       r(I)=b(I)-Sumation
    End Do
   
   !Part :
    Denom=Num
    Num=0
    Do I=1,NBP+3
       Num=(r(I)*r(I))+Num                 
    End Do
       
    Beta=Num/Denom
   
   !Part :
    Do I=1,NBP+3
       p(I)=r(I)+(Beta*p(I))    
    End Do
   
 End Do
!*********************************************************************************************
 End
!###########################################################################################   