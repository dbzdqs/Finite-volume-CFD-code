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
 Subroutine Extrude2DMesh(Dim,Nlayer,NF,NR,NFR,BC,NP,NC,SIDS,Corn,NEdgOfCell,FaceType,IDS)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Nlayer,SIDS,Corn,NEdgOfCell
 Intent(InOut)::NF,NR,NFR,NP,NC,BC
 Intent(Out  )::FaceType,IDS

 Integer::Dim,I,J,cnt,SF,NF,NR,NP,NC,Nlayer,FacInx,ILayer
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:4,1:Dim)::SIDS
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:4,1:Dim)::Corn
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:Dim)::NEdgOfCell
!*********************************************************************************************
!Part 1:
 FacInx = 0
 SF=0
 Do I=1,NR  
    Do ILayer=1,Nlayer

      Do J=SF+1,SF+NFR(I)	
       
	    !Part 2:  
	     FacInx = FacInx + 1
         FaceType(FacInx)=4

		!Part 3:
	     IDS(1,FacInx) = SIDS(1,J) + (ILayer-1)*NC
         IDS(2,FacInx) = SIDS(2,J) + (ILayer-1)*NC
         IF( SIDS(2,J)==0 ) IDS(2,FacInx) = 0

		!Part 4:
         IDS(3,FacInx) = SIDS(3,J) + (ILayer-1)*NP
	     IDS(4,FacInx) = SIDS(4,J) + (ILayer-1)*NP          
	     IDS(5,FacInx) = SIDS(4,J) + (ILayer  )*NP
	     IDS(6,FacInx) = SIDS(3,J) + (ILayer  )*NP
  
      End do
	 
    End do
	
    SF=SF+NFR(I)
 End do

!Part 5:
 Do ILayer=2,NLayer
    Do I=1,NC	
       
      !Part 6:	   
	   FacInx = FacInx + 1
       FaceType(FacInx)=NEdgOfCell(I)

	  !Part 7:
	   IDS(1,FacInx) = I + (ILayer-2)*NC
       IDS(2,FacInx) = I + (ILayer-1)*NC

	  !Part 8:
       IDS(3,FacInx) = Corn(1,I) + (ILayer-1)*NP
	   IDS(4,FacInx) = Corn(2,I) + (ILayer-1)*NP          
	   IDS(5,FacInx) = Corn(3,I) + (ILayer-1)*NP          
	   IDS(6,FacInx) = Corn(4,I) + (ILayer-1)*NP
     
    End do
 End do

!Part 9:
 Do I=1,NC	
    
   !Part 10:	
    FacInx = FacInx + 1
    FaceType(FacInx)=NEdgOfCell(I)

   !Part 11:
	IDS(1,FacInx) = I
    IDS(2,FacInx) = 0

   !Part 12:
    cnt=2
    Do J=NEdgOfCell(I),1,-1
       cnt=cnt+1
       IDS(cnt,FacInx) = Corn(J,I)
    End Do
 
 End do

!Part 13:
 Do I=1,NC	
     
   !Part 14:     
    FacInx = FacInx + 1
    FaceType(FacInx)=NEdgOfCell(I)

   !Part 15:
	IDS(1,FacInx) = I + (NLayer-1)*NC
    IDS(2,FacInx) = 0

   !Part 16:
    IDS(3,FacInx) = Corn(1,I) + NLayer*NP
	IDS(4,FacInx) = Corn(2,I) + NLayer*NP          
	IDS(5,FacInx) = Corn(3,I) + NLayer*NP         
	IDS(6,FacInx) = Corn(4,I) + NLayer*NP

 End do

!Part 17:
 Do I=1,NR  
    NFR(I) = NFR(I)*NLayer  
    BC(I) = BC(I)
 End do

!Part 18:
 IF(NLayer>1)Then
  NR=NR+1
  NFR(NR) = NC*(NLayer-1)
    BC(I) = 1
 EndIF

!Part 19:
 NR=NR+1
 NFR(NR) = NC
 BC(NR) = 6

 NR=NR+1
 NFR(NR) = NC
 BC(NR) = 6

!Part 20:
 NP = NP*(Nlayer+1)
 NF = NF*NLayer + NC*(NLayer+1)
 NC = NC*Nlayer
!*********************************************************************************************
 End
!###########################################################################################