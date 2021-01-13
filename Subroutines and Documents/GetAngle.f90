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
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function GetAngle(Dim,FrontBegining,FrontEnd,Corner,X,Y)
Implicit None
!===========================================================================================
Intent(In)::Dim,FrontBegining,FrontEnd,Corner,X,Y

Integer::Dim,FrontBegining,FrontEnd,Corner
Real(8)::GetAngle,MFU,MFV,U,V,InnerProduct,Temp,NMF,Norm
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
MFU=X(FrontEnd) - X(FrontBegining)	       ! Considering Main Front as a vector in form of (U + V) to give direction to it
MFV=Y(FrontEnd) - Y(FrontBegining)		   

U=X(Corner) - X(FrontBegining)	           ! Considering other edge as a vector 
V=Y(Corner) - Y(FrontBegining)

!Part 2:
InnerProduct = MFU*U + MFV*V               ! Calculating Inner Product

Temp=MFU*MFU + MFV*MFV
NMF=DSQRT(Temp)			                   ! Norm of Main Front

Temp=U*U + V*V
Norm=DSQRT(Temp)		                   ! Norm of other vector

Temp=InnerProduct/(NMF*Norm)

if(Temp > 1) then
	Temp = 1.0
elseif(Temp < -1) then
	Temp = -1.0
endif

GetAngle = ACOS(Temp)
!===========================================================================================
End Function GetAngle
!*********************************************************************************************
