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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine LES_DSmag_Eddy3D(Dim,NC,NF1,NF2,IDS,Vol,MR,Rho,Rhohat,Skk,Skkhat,Sabs,Sabshat,Lkk,Mut,&
                              Sxx,Syy,Szz,Sxy,Sxz,Szy,&
                              Sxxhat,Syyhat,Szzhat,Sxyhat,Sxzhat,Szyhat,&
                              Lxx,Lyy,Lzz,Lxy,Lxz,Lzy)
 
 Implicit None
!********************************************************************************************* 
 Intent(In   )::Dim,NC,NF1,NF2,IDS,Vol,MR,Rho,Rhohat,Skk,Skkhat,Sabs,Sabshat,Lkk,&
                              Sxx,Syy,Szz,Sxy,Sxz,Szy,&
                              Sxxhat,Syyhat,Szzhat,Sxyhat,Sxzhat,Szyhat,&
                              Lxx,Lyy,Lzz,Lxy,Lxz,Lzy
 Intent(InOut)::Mut

 Integer::Dim,I,NC,NF1,NF2,Allocatestatus,DeAllocatestatus
 Real(8)::MR,Part1,Part2,Cd,Axx,Ayy,Azz,Axy,Axz,Azy,Mxx,Myy,Mzz,Mxy,Mxz,Mzy,Lxxd,Lyyd,Lzzd,Lxyd,Lxzd,Lzyd
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::Rho,Rhohat,Lkk,Vol,Mut,Sabshat,Sabs,Skk,Skkhat,&
                           Sxx,Syy,Szz,Sxy,Sxz,Szy,&
                           Sxxhat,Syyhat,Szzhat,Sxyhat,Sxzhat,Szyhat,&
                           Lxx,Lyy,Lzz,Lxy,Lxz,Lzy,&
                           Bxx,Byy,Bzz,Bxy,Bxz,Bzy,Bxxhat,Byyhat,Bzzhat,Bxyhat,Bzyhat,Bxzhat,&
                           LijMij,MijMij,LijMijhat,MijMijhat
!********************************************************************************************* 
!Part 1:
 Do I=1,NC

    Part1 =-2.0*Rho(I) * (Vol(I)**0.6666) * Sabs(I)
	Part2 = 0.3333*Skk(I)

    Bxx(I) = Part1 * ( Sxx(I) - Part2 )
    Byy(I) = Part1 * ( Syy(I) - Part2 )
    Bzz(I) = Part1 * ( Szz(I) - Part2 )
    
    Bxy(I) = Part1 *   Sxy(I)
    Bxz(I) = Part1 *   Sxz(I)
    Bzy(I) = Part1 *   Szy(I)

 End Do

!Part 2:
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Bxx,Vol,Bxxhat)
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Byy,Vol,Byyhat)
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Bzz,Vol,Bzzhat)
 
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Bxy,Vol,Bxyhat)
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Bxz,Vol,Bxzhat)
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Bzy,Vol,Bzyhat)

!Part 3:
 Do I=1,NC

   !Part 4:
    Part1 =-2.0*Rhohat(I) * 4*(Vol(I)**0.6666) * Sabshat(I)
	Part2 = 0.3333*Skkhat(I)

    Axx = Part1 * ( Sxxhat(I) - Part2 )
    Ayy = Part1 * ( Syyhat(I) - Part2 )
    Azz = Part1 * ( Szzhat(I) - Part2 )
    
    Axy = Part1 *   Sxyhat(I)
    Axz = Part1 *   Sxzhat(I)
    Azy = Part1 *   Szyhat(I)

   !Part 5:
    Mxx = Axx-Bxxhat(I)
    Myy = Ayy-Byyhat(I)
    Mzz = Azz-Bzzhat(I)
    
    Mxy = Axy-Bxyhat(I)
    Mxz = Axz-Bxzhat(I)
    Mzy = Azy-Bzyhat(I)

   !Part 6:
    Lxxd= Lxx(I) - 0.3333*Lkk(I)
    Lyyd= Lyy(I) - 0.3333*Lkk(I) 
    Lzzd= Lzz(I) - 0.3333*Lkk(I)
    
    Lxyd= Lxy(I) 
    Lxzd= Lxz(I) 
    Lzyd= Lzy(I) 

   !Part 7:
    LijMij(I) = Lxxd*Mxx + Lyyd*Myy + Lzzd*Mzz + 2*Lxyd*Mxy + 2*Lxzd*Mxz + 2*Lzyd*Mzy
    MijMij(I) = Mxx*Mxx  + Myy*Myy  + Mzz*Mzz  + 2*Mxy*Mxy  + 2*Mxz*Mxz  + 2*Mzy*Mzy    

 End Do

!Part 8:
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,LijMij,Vol,LijMijhat)
 Call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,MijMij,Vol,MijMijhat)

 Do I=1,NC
   
   !Part 9: 
    Cd=LijMijhat(I)/MijMijhat(I) 
   
   !Part 10:    
    If (Cd<0.) Then
      Cd=0.
    Else If (Cd>=0.0324) Then
      Cd=0.0324;
    Else If (Cd>=0.0 .And. Cd<=0.0324) Then
      Cd=Cd
    Else
      Cd=0.0324;
    End If     
   
   !Part 11: 
    Mut(I) = (1.0/MR) * Cd * (Vol(I)**0.6666) * Rho(I)*Sabs(I)

 End Do

!*********************************************************************************************
 End 
!########################################################################################### 
