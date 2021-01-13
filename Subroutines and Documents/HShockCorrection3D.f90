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
!// Date: May., 15, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine HShockCorrection3D(Dim,NF,NP,IDS,FaceType,X,Y,Z,NodeMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NF,NP,I,J,K,Face,P,Q,temp,counter
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:3,1:3,1:Dim)::NodeMetric
 Real(8)::Length,HShock,Beta,Factor,hP,hQ,Norm !hP: Local Size of P in the Q direction
 Real(8),Dimension(1:Dim)::X,Y,Z
 Logical::correction
 Integer,Dimension(1:4)::FacePoints
!*********************************************************************************************
!Part 1:
 Beta = 2.0
 correction = .TRUE.
 counter = 0
 Do while (correction)
   !Part 2:
    Do I=1,NF
       Face = I
       FacePoints(1:FaceType(Face)) = IDS(3:(2+FaceType(Face)),Face)
       Do J=1,FaceType(Face)    
          P = FacePoints(J)
          Q = FacePoints(Mod(J,FaceType(Face)) + 1)
           
          Do K=1,2
             If(K==2)Then
              temp = P
              P    = Q
              Q    = temp  
             Endif
           
            !Part 3:
             Norm = sqrt(( (X(Q)-X(P)) **2)+( (Y(Q)-Y(P)) **2)+( (Z(Q)-Z(P)) **2))
             Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,NodeMetric(:,:,P),P,Q,hP)
             Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,NodeMetric(:,:,Q),Q,P,hQ)
          
            !Part 4:
             If(hQ<hP) Cycle
           
             If(abs(hP-hQ)<=0.00000000001)Then
              Length = Norm/hP
             Else
              Length = Norm*((hQ-hP) / (hP*hQ*( Log(hQ/hP) )))
             Endif
           
            !Part 5:
             HShock = (MAX( (hQ/hP) , (hP/hQ) )**(1.0/Length))
             If((HShock-Beta)>0.00000000001)Then
              Factor   =   1.0/ ( (Beta/HShock)**(Length**2.0) )
              NodeMetric(:,:,Q) = Factor * NodeMetric(:,:,Q)
             Endif
           
          EndDo
       EndDo   
    EndDo
    
   !Part 6:
    correction = .FALSE.
    Do I=1,NF
       Face = I
       FacePoints(1:FaceType(Face)) = IDS(3:(2+FaceType(Face)),Face)
       Do J=1,FaceType(Face)    
          P = FacePoints(J)
          Q = FacePoints(Mod(J,FaceType(Face)) + 1)
          
          Do K=1,2
             If(K==2)Then
              temp = P
              P    = Q
              Q    = temp  
             Endif
          
            !Part 7:
             Norm = sqrt(( (X(Q)-X(P)) **2)+( (Y(Q)-Y(P)) **2)+( (Z(Q)-Z(P)) **2))
             Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,NodeMetric(:,:,P),P,Q,hP)
             Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,NodeMetric(:,:,Q),Q,P,hQ)
            
             If(abs(hP-hQ)<=0.00000000001)Then
              Length = Norm/hP
             Else
              Length = Norm*((hQ-hP) / (hP*hQ*( Log(hQ/hP) )))
             Endif
           
            !Part 8:
             HShock = (MAX( (hQ/hP) , (hP/hQ) )**(1.0/Length))
             If((HShock-Beta)>0.00000000001)Then
              correction = .TRUE.
              exit
             Endif
           
          EndDo
          If(correction==.TRUE.)exit
       EndDo
       If(correction==.TRUE.)exit
    EndDo
    counter = counter + 1
    If(counter>10) exit
 EndDo
!*********************************************************************************************
 End
!###########################################################################################