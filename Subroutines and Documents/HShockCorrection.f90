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
!// Date: Mar., 05, 2013                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine HShockCorrection(Dim,NF,NP,IDS,X,Y,NodeMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NF,NP,I,J,K,Face,P,Q,counter=0
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Real(8)::hP,hQ,Norm
 Real(8)::Length,HShock,Beta,Factor
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1,2)::UnitV
 Real(8),Dimension(2,1)::UnitVt
 Integer,Dimension(1:2)::MDim1,MDim2
 Real(8),Dimension(1)::TempMatrix1,TempMatrix3
 Real(8),Dimension(1:2,1)::TempMatrix2
 Logical::correction
!*********************************************************************************************
!Part 1:
 Beta = 2
 correction = .TRUE.
 Do while (correction)
    counter=counter+1
    If(counter>200)exit
   !Part 2:
    Do I=1,NF
       Face = I
       Do K=1,2
          If(K==1)Then
           P = IDS(3,Face)
           Q = IDS(4,Face)
          Else
           Q = IDS(3,Face)
           P = IDS(4,Face)
          Endif
          
         !Part 3:
          Norm = sqrt(( (X(Q)-X(P)) **2)+( (Y(Q)-Y(P)) **2))
          Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,NodeMetric(:,:,P),P,Q,hP)
          Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,NodeMetric(:,:,Q),Q,P,hQ)
          
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
    
   !Part 6:
    correction = .FALSE.
    Do I=1,NF
       Face = I
       Do K=1,2
          If(K==1)Then
           P = IDS(3,Face)
           Q = IDS(4,Face)
          Else
           Q = IDS(3,Face)
           P = IDS(4,Face)
          Endif
          
         !Part 7:
          Norm = sqrt(( (X(Q)-X(P)) **2)+( (Y(Q)-Y(P)) **2))
          Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,NodeMetric(:,:,P),P,Q,hP)
          Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,NodeMetric(:,:,Q),Q,P,hQ)
            
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
 EndDo
!*********************************************************************************************
 End
!###########################################################################################