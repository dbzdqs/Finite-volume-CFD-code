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
Subroutine CalcGradientDirction(Dim,Corn,X,Y,V,DELTA,Gradient,PreDistortionMetrics,Mu_min,E_min,QElms,QEC,COINCIDENT_TOLERANCE)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,V,DELTA,PreDistortionMetrics,Mu_min,E_min,QElms,QEC,COINCIDENT_TOLERANCE
Intent(Out)::Gradient
Intent(InOut)::X,Y

Integer,Parameter::Gx = 1
Integer,Parameter::Gy = 2

Integer::Dim,V,I,J,State,E_min,QEC,choice
Integer,Dimension(1:1000)::QElms
Integer,Dimension(1:Dim,1:4)::Corn
Real(8)::DELTA,Improve,Mu,Mu_min,MuMinPlus,MuPlus,x_value,y_value,QuadDistortionMetric,COINCIDENT_TOLERANCE,randValue
Real(8),Dimension(1:8)::Mu_min_plus
Real(8),Dimension(1:100)::PreDistortionMetrics,PostDistortionMetrics
Real(8),Dimension(1:Dim)::X,Y
Real(8),Dimension(1:100,1:2)::Gradient
!===========================================================================================
!Part 1:

x_value = X(V)
y_value = Y(V)

!Part 2:

do I=1,QEC !---------------------- Perturbing x-coordinates of elements -------------------
    
    Mu = QuadDistortionMetric(Dim,Corn,X,Y,QElms(I),COINCIDENT_TOLERANCE)
        
    Call RANDOM_NUMBER(randValue)
    randValue = randValue*2
    choice = 1 + FLOOR(randValue)

    Select Case(choice)
    
        Case(1)
        
            X(V) = x_value + DELTA
            
            !--------------------- Calculate post distortion metric ------------------
                
            MuPlus = QuadDistortionMetric(Dim,Corn,X,Y,QElms(I),COINCIDENT_TOLERANCE)     

            Gradient(I,Gx) = (MuPlus-Mu)/DELTA
            
        Case(2)
            
            X(V) = x_value - DELTA
            
            !--------------------- Calculate post distortion metrics ------------------
                
            MuPlus = QuadDistortionMetric(Dim,Corn,X,Y,QElms(I),COINCIDENT_TOLERANCE)     

            Gradient(I,Gx) = (MuPlus-Mu)/DELTA
            
    End Select

end do

X(V) = x_value

!Part 3:

do I=1,QEC !---------------------- Perturbing y-coordinates of elements -------------------
    
    Mu = QuadDistortionMetric(Dim,Corn,X,Y,QElms(I),COINCIDENT_TOLERANCE)
        
    Call RANDOM_NUMBER(randValue)
    randValue = randValue*2
    choice = 1 + FLOOR(randValue)

    Select Case(choice)
    
        Case(1)
        
            Y(V) = y_value + DELTA
            
            !--------------------- Calculate post distortion metric ------------------
                
            MuPlus = QuadDistortionMetric(Dim,Corn,X,Y,QElms(I),COINCIDENT_TOLERANCE)     

            Gradient(I,Gy) = (MuPlus-Mu)/DELTA
            
        Case(2)
            
            Y(V) = y_value - DELTA
            
            !--------------------- Calculate post distortion metrics ------------------
                
            MuPlus = QuadDistortionMetric(Dim,Corn,X,Y,QElms(I),COINCIDENT_TOLERANCE)     

            Gradient(I,Gy) = (MuPlus-Mu)/DELTA
            
    End Select

end do

Y(V) = y_value

!===========================================================================================
End Subroutine CalcGradientDirction 
!*********************************************************************************************
