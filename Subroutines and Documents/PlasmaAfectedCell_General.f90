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
!// Developed by: S. Kabirian, Aerospace Eng., Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine PlasmaAfectedCell_General(Dim,NC,Xc,Yc,Delta)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,Xc,Yc 
 Intent(Out  )::Delta

 Integer :: Dim,I,NC,Intersect,M,J,JP1
 Real(8),Dimension(1:Dim)::Xc,Yc
 Integer,Dimension(1:Dim)::Delta
 Real(8),Dimension(1:5)::XX,YY
 Real(8)::X1,X2,X3,X4,Y1,Y2,Y3,Y4,Xcros,Ycros
!=============================================================================================================
 Delta(:)=0

 XX(1)=0.02
 YY(1)=0.0294
 
 XX(2)=0.0170
 YY(2)=0.0334
 
 XX(3)=0.0213
 YY(3)=0.0362
 
 XX(4)=0.0242
 YY(4)=0.0320
 
 XX(5)=0.0209
 YY(5)=0.0331
 
 DO i=1,NC
    
    !For Flat Plate, Test Case 01
	!IF (Xc(i)>-0.481 .AND. Xc(i)<-0.45121 .AND. Yc(i)>0.034413 .AND. Yc(i)<0.0500 ) Delta(i)=1
  
    !For Flat Plate, Test Case 02
     !IF (Xc(i)>0.58536 .AND. Xc(i)<0.61536 .AND. Yc(i)>0.0 .AND. Yc(i)<0.015 .AND. Yc(i)<-0.5*Xc(i)+0.30768 ) Delta(i)=1  
	 
    !For NACA0012
    
    !For NACA0015
    
     !IF (Xc(i)>0.284392 .AND. Xc(i)<0.337302 .AND. Yc(i)>0.07 .AND. Yc(i)<0.128 ) Delta(i)=1  
	 
     X1=Xc(i)
     Y1=Yc(i)
     X2=XX(5)
     Y2=YY(5)
     
     M=0
     Do J=1,4
        X3=XX(J)
        Y3=YY(J)
        JP1=J+1
        IF(J==4)JP1=1
        X4=XX(JP1)
        Y4=YY(JP1)
        
        Call Cross_Lines(X1,X2,X3,X4,Y1,Y2,Y3,Y4,Xcros,Ycros,Intersect)
        If(Intersect==1) M=M+1

     End Do

     If( (M/2)/=(M/2.0) ) Delta(i)=1 
    
 END DO
!*********************************************************************************************
 End
!###########################################################################################
