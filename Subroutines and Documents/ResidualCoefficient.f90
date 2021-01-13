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
!// Developed by: S. Sheikhi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ResidualCoefficient(Dim,NC,NF1,NF2,IDS,EPS,ia,ja,a,cnt)
 Implicit None
!*********************************************************************************************
 Intent(In)::Dim,NC,NF1,NF2,IDS,EPS
 Intent(Out  )::ia,ja,a,cnt

 Integer::Dim,I,ME,NE,NC,NF1,NF2,cnt
 Real(8)::EPS
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::ia,ja,Nb
 Real(8),Dimension(1:Dim)::a
!*********************************************************************************************
!Part 1:
 cnt=0

!Part 2: 
 Do I=1,NC
    Nb(I)=0    
 end Do

!Part 3:
 Do I=NF1+1,NF2
 
   !Part 4: 
    ME=IDS(1,I)
    NE=IDS(2,I)
   
    Nb(ME)=Nb(ME)+1
    Nb(NE)=Nb(NE)+1
  
   !Part 5: 
    cnt=cnt+1
    a(cnt)=-EPS
    ia(cnt)=ME
    ja(cnt)=NE
   
   !Part 6:
    cnt=cnt+1
    a(cnt)=-EPS
    ia(cnt)=NE
    ja(cnt)=ME		 
   
 end Do

!Part 7:
 Do I=1,NC
    cnt=cnt+1
    a(cnt)=1+Nb(I)*EPS
    ia(cnt)=I
    ja(cnt)=I 
 end Do
!********************************************************************************************* 
 END
!###########################################################################################