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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine SortFacesOfRegion2D(Dim,IndxReg,NR,NFR,IDS)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,IndxReg,NR,NFR
 Intent(InOut)::IDS

 Integer::Dim,I,J,J1,P1_2,P2_2,ME_2,NE_2,P2_1,Temp,Indx
 Integer::FirstFace,LastFace
 Integer::IndxStartingFace
 Integer::NR
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer::IndxReg
!*********************************************************************************************
!Part 1:
 FirstFace = 0
 Do I=1,IndxReg-1
    FirstFace = FirstFace + NFR(I)
 End Do
 LastFace = FirstFace + NFR(IndxReg)
 
 print*,'please check if direction of the first edge is toward mid of curve (find input.plt file)'
 print*,'First Pt:',IDS(3,FirstFace+1),'Last Pt:',IDS(4,FirstFace+1)
 print*,'(Yes=1/No=0)'
 read*,i

!Part 2: 
 if(i==0)then
 Temp               = IDS(1,FirstFace+1)
 IDS(1,FirstFace+1) = IDS(2,FirstFace+1)
 IDS(2,FirstFace+1) = Temp 

 Temp               = IDS(3,FirstFace+1)
 IDS(3,FirstFace+1) = IDS(4,FirstFace+1)
 IDS(4,FirstFace+1) = Temp
 endif
 
!Part 3:
 Do I=FirstFace+1,LastFace-1

   !Part 4:
    P2_1=IDS(4,I)

   !Part 5:
	Do J=I+1,LastFace

      !Part 6:         
       ME_2 = IDS(1,J)
       NE_2 = IDS(2,J)
       P1_2 = IDS(3,J)  
       P2_2 = IDS(4,J)    

      !Part 7:
	   If(P2_1==P1_2) Exit

      !Part 8:
	   If(P2_1==P2_2)then  
        IDS(1,J) = NE_2
        IDS(2,J) = ME_2
        IDS(3,J) = P2_2
        IDS(4,J) = P1_2
		Exit
	   Endif

	End do

 End do
!*********************************************************************************************
 End 
!###########################################################################################

