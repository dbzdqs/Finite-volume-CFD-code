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
!// Developed by: M. Keley, Mechanical Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine NearstCellProfile(Dim,NC,NFW1,NFW2,DW,INW,NNearstCell,INearstCell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,DW,INW
 Intent(Out  )::NNearstCell,INearstCell

 Integer::Dim,J,I,II,P1,P2,ME,NC,NFW1,NFW2,NFS1,NFS2,K,JJ,I1,SS
 Real(8)::Dmin,Dis,Xj,Yj,Xi,Yi,DX,DY,Mu0,Mut0,Y1,Y2,E,Yn,dm
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:Dim)::DW,S
 Integer,Dimension(NFW1:NFW2,1:12000)::INearstCell
 Integer,Dimension(NFW1:NFW2)::NNearstCell
!*********************************************************************************************	
!Part 1:
 NNearstCell(:)=0
 INearstCell(:,:)=0

!Part 2:
 DO K=NFW1+1,NFW2  
    DO J=1,NC

        I  = INW(J)
	   IF(I/=K) CYCLE

       !Part 3:
	    NNearstCell(K)=NNearstCell(K)+1
        INearstCell(K,NNearstCell(K)) = J

     !Part 4:
      DO JJ=NNearstCell(K),2,-1
          
       Y1 = DW( INearstCell(K,JJ  ) )
       Y2 = DW( INearstCell(K,JJ-1) )
	    
	   IF( Y1<Y2 ) THEN
		I1=INearstCell(K,JJ  )
		INearstCell(K,JJ  ) = INearstCell(K,JJ-1)
		INearstCell(K,JJ-1)=I1
	   Endif

      END DO

    END DO

 END DO
!*********************************************************************************************
 End
!###########################################################################################
