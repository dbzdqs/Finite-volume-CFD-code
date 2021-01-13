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
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function MoveIsAcceptable(Dim,Corn,Neib,X,Y,TElms,TEC,QElms,QEC,PreTMetrics,PostTMetrics,PreQMetrics,PostQMetrics)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,X,Y,TElms,TEC,QElms,QEC,PreTMetrics,PostTMetrics,PreQMetrics,PostQMetrics

Integer,Parameter::MyMin = 0.05

Integer::Dim,Nimp,Nwrs,Nimsig,Nwrsig,Ninv,TEC,QEC,E,I,J,P
Integer,Dimension(1:1000)::TElms,QElms
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::MoveIsAcceptable,ElementOverSpanned,TriIsInverted,QuadIsInverted
Real(8)::PI,THETA_MAX,Pre,Post,AMC,MinDisMetric,Theta
Real(8),Dimension(1:100)::PreTMetrics,PostTMetrics,PreQMetrics,PostQMetrics
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
MoveIsAcceptable = .False.

PI = 3.14159265358979d0

THETA_MAX =(PI/180)*200 !----------- Maximum Angle Allowed in Quads (200 Degrees) ----------

AMC = 0.0 !---------------------------- Average Metric Change ------------------------------

ElementOverSpanned = .False.

MinDisMetric = PreQMetrics(1)  

Nimp = 0 !------------------ Number of elements whose metric improves ----------------------
Nwrs = 0 !------------------ Number of elements whose metric gets worse --------------------

Nimsig = 0 !------------- Number of elements whose metric improves significantly -----------
Nwrsig = 0 !------------- Number of elements whose metric worsens significantly ------------

Ninv = 0 !-------------------- Number of inverted elements created ------------------------- 

!Part 1:

do I=1,QEC

	Pre = PreQMetrics(I)
	Post = PostQMetrics(I)

	!------------------------- Specifying Number of Inverted Elements ----------------------

	E = QElms(I)

	if(QuadIsInverted(Dim,Corn,X,Y,E)) then

		Ninv = Ninv + 1

	endif

	!------------------- Detecting Element Over-Spanning -----------------------------------

	Call getLargestAngleOfQuad(Dim,Corn,Neib,X,Y,E,Theta,P)

	if(Theta > THETA_MAX) then

		ElementOverSpanned = .True.

	endif

	!-------------------- Calculating Average Metric Change --------------------------------

	AMC = AMC + (Post - Pre)

	!------------------------- Specifying Nimp and Nimsig Values ---------------------------

	if(Pre < Post) then

		Nimp = Nimp + 1 

		if(Pre < 0 .Or. (Pre < MyMIN .And. Post >= MyMin)) then

			Nimsig = Nimsig + 1

		endif

	endif

	!------------------------- Specifying Nwrs and Nwrsig Values ---------------------------

	if(Pre > Post) then

		Nwrs = Nwrs + 1

		if(Post < 0 .Or. Post < MyMIN) then

			Nwrsig = Nwrsig + 1

		endif

	endif

end do

do I=1,TEC

	Pre = PreTMetrics(I)
	Post = PostTMetrics(I)

	!------------------------- Specifying Number of Inverted Elements ----------------------

	E = TElms(I)
        
	if(TriIsInverted(Dim,Corn,X,Y,E)) then

		Ninv = Ninv + 1

	endif

	!-------------------- Calculating Average Metric Change --------------------------------

	AMC = AMC + (Post - Pre)

	!------------------------- Specifying Nimp and Nimsig Values ---------------------------

	if(Pre < Post) then

		Nimp = Nimp + 1 

		if(Pre < 0 .Or. (Pre < MyMIN .And. Post >= MyMin)) then

			Nimsig = Nimsig + 1

		endif

	endif

	!------------------------- Specifying Nwrs and Nwrsig Values ---------------------------

	if(Pre > Post) then

		Nwrs = Nwrs + 1

		if(Post < 0 .Or. Post < MyMIN) then

			Nwrsig = Nwrsig + 1

		endif

	endif

end do

!Part 2:

AMC = AMC/(QEC + TEC) !---------------------- Average Metric Change ------------------------- 

if(Nwrs == (QEC+TEC) .Or. Ninv > 0 .Or. Nwrsig > Nimsig .Or. ElementOverSpanned .Or. (AMC < (MyMIN*(-1)))) then

	MoveIsAcceptable = .False.

elseif(Nimp == (QEC+TEC) .Or. (Nimsig > 0 .And. Nwrsig == 0) .Or. (Nimsig >= Nwrsig .And. (AMC > (MyMIN*(-1))))) then

	MoveIsAcceptable = .True.

endif 

!===========================================================================================
End Function MoveIsAcceptable
!*********************************************************************************************
