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
 Subroutine DetectSepetRegnOfMesh2D(Dim,NF,IDS,NR,NFR)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF
 Intent(InOut)::IDS
 Intent(Out  )::NR,NFR

 Integer::Dim,I,J,J1,M,P1,P2,NBoundFac,NF,NBoundCurves,Temp,NR
 Integer,Dimension(1:100)::NFeceOneachCurves
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:100)::NFR
!*********************************************************************************************
!Part 1:
 NBoundFac=0
 Do I=1,NF

    If( IDS(2,I)==0 ) Then
	 NBoundFac=NBoundFac+1
     Do J=1,4
        Temp       = IDS(J,I)
        IDS(J,I)   = IDS(J,NBoundFac)
        IDS(J,NBoundFac) = Temp
     End Do 
    Endif

 End Do


 !!!J=40
 !!!Do J1=1,4
 !!!   M           = IDS(1,J1)
 !!!   IDS(1,J1)   = IDS(J,J1)
 !!!   IDS(J,J1) = M
 !!!End do

!Part 2:
 NBoundCurves=1
 NFeceOneachCurves(1) = 1

!Part 3:
 Do I=1,NBoundFac-1

   !Part 4:
    P2=IDS(4,I)

   !Part 5:
	Do J=I+1,NBoundFac

      !Part 6:
       P1=IDS(3,J)     

      !Part 7:
	   If(P2==P1)then  

	    NFeceOneachCurves(NBoundCurves) = NFeceOneachCurves(NBoundCurves)+1

        Do J1=1,4
           M           = IDS(J1,I+1)
           IDS(J1,I+1) = IDS(J1,J)
           IDS(J1,J)   = M
		End do

		Go to 10
	   Endif

	End do

   !Part 8:
    NBoundCurves=NBoundCurves+1
    NFeceOneachCurves(NBoundCurves) = 1
    
10 End do
   
 NFR(:) = NFeceOneachCurves(:)
 NR = NBoundCurves+1
 NFR(NR) = NF - NBoundFac
   
!*********************************************************************************************
 End 
!###########################################################################################

