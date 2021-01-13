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
 Subroutine CopyMesh(Dim,NF,IDS,NX,NY,DA,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
                     Tmp_IDS,Tmp_NX,Tmp_NY,Tmp_DA,Tmp_NF,FaceNum)
 Implicit None
!********************************************************************************************* 
 Intent(In   )::Dim,NF,IDS,NX,NY,DA,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2
 Intent(Out  )::Tmp_IDS,Tmp_NX,Tmp_NY,Tmp_DA,Tmp_NF,FaceNum

 Integer::Dim,I,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,Tmp_NF,NF
 Real(8),Dimension(1:Dim)::NX,NY,DA,Tmp_NX,Tmp_NY,Tmp_DA
 Integer,Dimension(1:4,1:Dim)::IDS,Tmp_IDS
 Integer,Dimension(1:10,1:2)::FaceNum
!*********************************************************************************************	
 !Part 1:
 Tmp_NF = NF

!Part 2:
 Do I=1,NF
   Tmp_NX(I)=NX(I)
   Tmp_NY(I)=NY(I)
   Tmp_DA(I)=DA(I)
   Tmp_IDS(1,I)=IDS(1,I)
   Tmp_IDS(2,I)=IDS(2,I)
   Tmp_IDS(3,I)=IDS(3,I)
   Tmp_IDS(4,I)=IDS(4,I)
 End Do	   

!Part 3: 
 FaceNum(1,1)  = NF1    ;   FaceNum(1,2)  = NF2
 FaceNum(2,1)  = NFW1   ;   FaceNum(2,2)  = NFW2
 FaceNum(3,1)  = NFF1   ;   FaceNum(3,2)  = NFF2
 FaceNum(4,1)  = NFI1   ;   FaceNum(4,2)  = NFI2
 FaceNum(5,1)  = NFO1   ;   FaceNum(5,2)  = NFO2
 FaceNum(6,1)  = NFS1   ;   FaceNum(6,2)  = NFS2
 FaceNum(7,1)  = NFIF1  ;   FaceNum(7,2)  = NFIF2
   
!*********************************************************************************************
 End
!###########################################################################################