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
!// Developed by: S. Sheikhi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine VectorizeMesh(Dim,NPassCell,IPassCell,NCELL_EDGE,CELL_EDGE,Tmp_IDS,Tmp_NX,Tmp_NY,Tmp_DA,Tmp_NF,FaceNum,&
                     IDS,NX,NY,DA,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFIF1,NFIF2)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NPassCell,IPassCell,NCELL_EDGE,CELL_EDGE,Tmp_IDS,Tmp_NX,Tmp_NY,Tmp_DA,Tmp_NF,FaceNum
 Intent(Out  )::IDS,NX,NY,DA,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFIF1,NFIF2

 Integer::Dim,NC,J,I,K,NPassCell,NFN,NFW,NFF,NFI,NFO,NFS,NFIF,Tmp_NF,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,&
          NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFIF1,NFIF2,E
  
 Real(8),Dimension(1:Dim)::NX,NY,DA,Tmp_NX,Tmp_NY,Tmp_DA
 Integer,Dimension(1:4,1:Dim)::IDS,Tmp_IDS,CELL_EDGE
 Integer,Dimension(1:Dim)::IPassCell,NCELL_EDGE
 Integer,Dimension(1:Dim)::NonFreezFace
 Integer,Dimension(1:10,1:2)::FaceNum
!*********************************************************************************************
 !Part 1:
 Do J=1,Tmp_NF
    NonFreezFace(J) = 0
 End Do
  	
 Do I=1,NPassCell
   !Part2:
    K=IPassCell(I)

	!Part3:
    Do J=1,NCELL_EDGE(K)
	   E=CELL_EDGE(J,K)
	   NonFreezFace( abs(E) ) = 1
	End Do
 End Do

!Part 4:
 NF=0
 Do I=1,Tmp_NF
    
    IF(NonFreezFace(I)==1)Then
	NF=NF+1

    NX(NF)=Tmp_NX(I)
    NY(NF)=Tmp_NY(I)
    DA(NF)=Tmp_DA(I)
    IDS(1,NF)=Tmp_IDS(1,I)
    IDS(2,NF)=Tmp_IDS(2,I)
    IDS(3,NF)=Tmp_IDS(3,I)
    IDS(4,NF)=Tmp_IDS(4,I)
	EndIF
 End Do
     

!Part 5: 
 NF1  = FaceNum(1,1)   ;    NF2  = FaceNum(1,2)
 NFW1 = FaceNum(2,1)   ;    NFW2 = FaceNum(2,2)
 NFF1 = FaceNum(3,1)   ;    NFF2 = FaceNum(3,2)
 NFI1 = FaceNum(4,1)   ;    NFI2 = FaceNum(4,2)
 NFO1 = FaceNum(5,1)   ;    NFO2 = FaceNum(5,2)
 NFS1 = FaceNum(6,1)   ;    NFS2 = FaceNum(6,2)
 NFIF1= FaceNum(7,1)   ;    NFIF2= FaceNum(7,2) 

!Part 6:
 NFN=0
 Do I=NF1+1,NF2
    IF(NonFreezFace(I)==1) NFN=NFN+1
 End Do

 NFW=0
 Do I=NFW1+1,NFW2
    IF(NonFreezFace(I)==1) NFW=NFW+1
 End Do

 NFF=0
 Do I=NFF1+1,NFF2
    IF(NonFreezFace(I)==1) NFF=NFF+1
 End Do

 NFI=0
 Do I=NFI1+1,NFI2
    IF(NonFreezFace(I)==1) NFI=NFI+1
 End Do

 NFO=0
 Do I=NFO1+1,NFO2
    IF(NonFreezFace(I)==1) NFO=NFO+1
 End Do

 NFS=0
 Do I=NFS1+1,NFS2
    IF(NonFreezFace(I)==1) NFS=NFS+1
 End Do

 NFIF=0
 Do I=NFIF1+1,NFIF2
    IF(NonFreezFace(I)==1) NFIF=NFIF+1
 End Do


!Part 7: 
 NF1=0      ; NF2=NF1+NFN
 NFW1=NF2   ;  NFW2=NFW1+NFW
 NFF1=NFW2  ; NFF2=NFF1+NFF
 NFI1=NFF2  ; NFI2=NFI1+NFI
 NFO1=NFI2  ; NFO2=NFO1+NFO
 NFS1=NFO2  ; NFS2=NFS1+NFS
 NFIF1=NFS2 ; NFIF2=NFIF1+NFIF
!*********************************************************************************************
 End
!###########################################################################################