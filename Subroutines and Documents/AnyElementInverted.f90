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
!// Date: Feb., 5, 2017                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function AnyElementInverted(Dim,Corn,X,Y,TEC,TElms,QEC,QElms)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,TEC,TElms,QEC,QElms

Integer::Dim,TEC,QEC,Elm,I
Integer,Dimension(1:1000)::TElms,QElms
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::AnyElementInverted,QuadIsInverted,TriIsInverted
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
AnyElementInverted = .False.
!Part 1:
do I=1,QEC
    
    Elm = QElms(I)
    
    if(QuadIsInverted(Dim,Corn,X,Y,Elm)) then
        
        AnyElementInverted = .True.
        exit
        
    endif
    
end do
!Part 2:
if(.Not. AnyElementInverted) then

    do I=1,TEC
    
        Elm = TElms(I)
        
        if(TriIsInverted(Dim,Corn,X,Y,Elm)) then
        
            AnyElementInverted = .True.
            exit
            
        endif
        
    end do
    
endif
!===========================================================================================
End Function AnyElementInverted
!*********************************************************************************************
