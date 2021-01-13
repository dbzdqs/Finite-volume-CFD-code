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
Subroutine AddElements(Dim,NC,NP,NBE,BFP,Corn,Neib,X,Y)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP
Intent(InOut)::NC,NP,Corn,Neib,X,Y

Integer::Dim,NC,NP,NBE,I
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::done,ElementInversionOccured,inversion
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1: 
do 
    do I=1,NC
        if(Corn(I,4)/=0 .And. Corn(I,1)/=Corn(I,2)) then
            Call ElementOpen(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,I,done)

            if(done) then
                print *,'I: ',I,' ElementOpen'
                exit
            endif
            
        endif
        
    end do
    
    if(I > NC) exit
    
end do
!===========================================================================================
End Subroutine AddElements 
!*********************************************************************************************
