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
 Subroutine PlasmaShayyParameters(freq,PlasmaHeight,Plasmawidth,appliedVoltage,PlasmaGap,Roh_c,e_c,Eb,DischargeTime)
 Implicit None
!********************************************************************************************* 
 Intent(Out)::freq,PlasmaHeight,Plasmawidth,appliedVoltage,PlasmaGap,Roh_c,e_c,Eb,DischargeTime

 Real(8)::freq
 Real(8)::PlasmaHeight      !height of Plasma in Y-Axis (ND)
 Real(8)::Plasmawidth       !width of Plasma in X-Axis (ND)
 Real(8)::appliedVoltage    !applied voltage (kv)
 Real(8)::PlasmaGap         !Distance Between the plates (cm)
 Real(8)::Roh_c		        !Density of Electron (/cm^2)
 Real(8)::e_c		        !charge of Electron (c)
 Real(8)::Eb			    !Breakdown Electric Field Strength (kv/cm)
 Real(8)::DischargeTime     !Discharge Time (sec)
!*********************************************************************************************
!INPUTS (All parameters must be in SI format)
 PlasmaHeight   = 0.018          
 Plasmawidth    = 0.024           
 AppliedVoltage = 4.0*sqrt(2.0)   
 PlasmaGap      = 0.025           
 Roh_c          = 1.0*10.0**11	
 e_c            = 1.602*10.0**-19 
 Eb             = 30.0			  
 DischargeTime  = 67.0*10.0** -6.0 
 freq           = 100
!*********************************************************************************************
 End
!###########################################################################################
