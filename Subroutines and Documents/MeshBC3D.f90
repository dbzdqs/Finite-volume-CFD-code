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
MODULE MeshBC3DHeader  
IMPLICIT NONE
INTERFACE
SUBROUTINE MeshBC3D(IntegerLocalArray,NR,NFR,BC,IDS,FaceType,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

 Integer,Dimension(:,:),pointer::IntegerLocalArray
 Integer::NR,NF,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NF1,NF2,NFIF1,NFIF2,NFN,NFW,NFF,NFI,NFO,NFS,NFIF
 Integer,Dimension(:)::NFR,BC
 Integer,Dimension(:,:)::IDS
 Integer,Dimension(:)::FaceType
 
END SUBROUTINE MeshBC3D
END INTERFACE   
END MODULE MeshBC3DHeader  
!*********************************************************************************************
 Subroutine MeshBC3D(IntegerLocalArray,NR,NFR,BC,IDS,FaceType,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 Implicit None
!*********************************************************************************************
 Intent(In   )::IntegerLocalArray,NR,NF
 Intent(Out  )::NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFIF1,NFIF2
 Intent(InOut)::NFR,BC,IDS,FaceType

 Integer,Dimension(:,:),pointer::IntegerLocalArray
 Integer::J,JJ,J1,I,SF,N,M,NR,NF,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,&
          NF1,NF2,NFIF1,NFIF2,NFN,NFW,NFF,NFI,NFO,NFS,NFIF
 
 Integer,Dimension(:)::NFR,BC
 Integer,Dimension(:,:)::IDS
 Integer,Dimension(:)::FaceType
 
 Integer,Dimension(1:NR)::TNFR,TBC
 Integer,Dimension(:,:),Pointer::TIDS
 Integer,Dimension(:),Pointer::TFaceType
 
 TFaceType=>IntegerLocalArray(1,:)
 TIDS=>IntegerLocalArray(2:7,:)
!*********************************************************************************************
!Part 1:
 Do J=1,NF
	Do J1=1,6
       TIDS(J1,J) = IDS(J1,J)
    End do
    TFaceType(J) = FaceType(J)
 End Do

!Part 2:
 Do J=1,NR
    TNFR(J) = NFR(J)
	TBC(J)  = BC(J) 
 End do

!Part 3:
 N=0
 M=0
 Do JJ=1,100

   !Part 4:
    SF=0
    Do J=1,NR
       IF(TBC(J)==JJ)Then

	    Do I=SF+1,SF+TNFR(J)

           N=N+1
           Do J1=1,6
              IDS(J1,N) = TIDS(J1,I)
		   End do
           FaceType(N) = TFaceType(I)

        Enddo
		
	   !Part 5:		   
	    M=M+1
	    NFR(M) = TNFR(J)
		BC(M)  = TBC(J)  

	   Endif
	   SF=SF+TNFR(J)
    End Do

 End Do

!Part 6:
 NFN  = 0
 NFW  = 0
 NFF  = 0
 NFI  = 0
 NFO  = 0
 NFS  = 0
 NFIF = 0
 Do J=1,NR
    IF( BC(J)==1 ) NFN  = NFN  + NFR(J)
    IF( BC(J)==2 ) NFW  = NFW  + NFR(J)
    IF( BC(J)==3 ) NFF  = NFF  + NFR(J)
    IF( BC(J)==4 ) NFI  = NFI  + NFR(J)
    IF( BC(J)==5 ) NFO  = NFO  + NFR(J)
    IF( BC(J)==6 ) NFS  = NFS  + NFR(J)
    IF( BC(J)==7 ) NFIF = NFIF + NFR(J)
 End Do

!Part 7:
 NF1=0
 NF2=NF1+NFN

 NFW1=NF2
 NFW2=NFW1+NFW

 NFF1=NFW2
 NFF2=NFF1+NFF

 NFI1=NFF2
 NFI2=NFI1+NFI

 NFO1=NFI2
 NFO2=NFO1+NFO

 NFS1=NFO2
 NFS2=NFS1+NFS

 NFIF1=NFS2
 NFIF2=NFIF1+NFIF
!*********************************************************************************************
 End
!###########################################################################################

 