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
!// Date: Mar., 10, 2015                                                                   //!
!// Developed by: S. Kabirian, Aerospace Eng., Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine PlasmaShayy(Dim,NC,freq,PlasmaHeight,Plasmawidth,appliedVoltage,PlasmaGap,&
                        Roh_c,e_c,Eb,DischargeTime,Xc,Yc,Delta, F_DBD_x,F_DBD_y)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,freq,PlasmaHeight,Plasmawidth,appliedVoltage,PlasmaGap,Roh_c,e_c,Eb,DischargeTime,Xc,Yc,Delta   
 Intent(Out  )::F_DBD_x,F_DBD_y

 Integer::Dim
 Integer::NC
 Integer::J
 Real(8)::kx,ky,E,E0,k1,k2 
 Real(8)::freq
 Real(8)::PlasmaHeight      !height of Plasma in Y-Axis (ND)
 Real(8)::Plasmawidth       !width of Plasma in X-Axis (ND)
 Real(8)::appliedVoltage    !applied voltage (kv)
 Real(8)::PlasmaGap         !Distance Between the plates (cm)
 Real(8)::Roh_c		        !Density of Electron (/cm^2)
 Real(8)::e_c		        !charge of Electron (c)
 Real(8)::Eb			    !Breakdown Electric Field Strength (kv/cm)
 Real(8)::DischargeTime     !Discharge Time (sec)
 Real(8)::alpha	
 Real(8),Dimension(1:Dim)::Xc,Yc
 Real(8),Dimension(1:Dim)::F_DBD_x,F_DBD_y
 Integer,Dimension(1:Dim)::delta
!********************************************************************************************* 
!Part 1:
 alpha = 1	
 E0    = appliedVoltage/PlasmaGap
 k1    = (E0-Eb)/Plasmawidth      !Plasma Parameter
 k2    = (E0-Eb)/PlasmaHeight     !Plasma Parameter
 
!Part 2:
 Do j=1,NC
 
   !Part 3:
    E=abs(E0-k1*Xc(j)-k2*Yc(j))

   !Part 4:
    kx=k2/sqrt(k1*k1 + k2*k2)
    ky=k1/sqrt(k1*k1 + k2*k2)

   !Part 5:
    F_DBD_x(j)=kx*freq*alpha*Roh_c*e_c*DischargeTime*E*delta(J)
    F_DBD_y(j)=ky*freq*alpha*Roh_c*e_c*DischargeTime*E*delta(J)

 End do
!*********************************************************************************************
 End
!###########################################################################################
