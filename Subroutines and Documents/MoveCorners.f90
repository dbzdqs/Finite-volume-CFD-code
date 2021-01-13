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
Function MoveCorners(Dim,NBE,BFP,Corn,X,Y,P,Q,List,Count)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,Corn,P,Q,List,Count
Intent(InOut)::X,Y

Integer::Dim,NBE,P,Q,Count,Elm,I
Integer,Dimension(1:1000)::List
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn
Logical::IsOnExteriorBoundary,MoveCorners,QuadIsInverted,TriIsInverted
Real(8)::x_value,y_value,xp,yp,xq,yq
Real(8),Dimension(1:Dim)::X,Y 
!===========================================================================================
MoveCorners = .True.
xp = X(P)
yp = Y(P)

xq = X(Q)
yq = Y(Q)

if(IsOnExteriorBoundary(Dim,NBE,BFP,P)) then !-------------- Move Q to P ------------------- 
!Part 1:
	X(Q) = X(P)
	Y(Q) = Y(P)
    
    do I=1,Count
        Elm = List(I)
        if(Corn(Elm,4) == 0) then
            if(TriIsInverted(Dim,Corn,X,Y,Elm)) then
                X(Q) = xq
	            Y(Q) = yq
                MoveCorners = .False.
                exit
            endif
        else
            if(QuadIsInverted(Dim,Corn,X,Y,Elm)) then
                X(Q) = xq
	            Y(Q) = yq
                MoveCorners = .False.
                exit    
            endif
        endif
    end do

elseif(IsOnExteriorBoundary(Dim,NBE,BFP,Q)) then !------------ Move P to Q -----------------
!Part 2:
	X(P) = X(Q)
	Y(P) = Y(Q)
    
    do I=1,Count
        Elm = List(I)
        if(Corn(Elm,4) == 0) then
            if(TriIsInverted(Dim,Corn,X,Y,Elm)) then
                X(P) = xp
	            Y(P) = yp
                MoveCorners = .False.
                exit
            endif
        else
            if(QuadIsInverted(Dim,Corn,X,Y,Elm)) then
                X(P) = xp
	            Y(P) = yp
                MoveCorners = .False.
                exit    
            endif
        endif
    end do

else !-------------------------- Move P to midpoint of PQ segment --------------------------
!Part 3:
	x_value = (X(P) + X(Q))/2
	y_value = (Y(P) + Y(Q))/2

	X(P) = x_value
	Y(P) = y_value
    
    X(Q) = x_value
    Y(Q) = y_value
    
    do I=1,Count
        Elm = List(I)
        if(Corn(Elm,4) == 0) then
            if(TriIsInverted(Dim,Corn,X,Y,Elm)) then
                X(P) = xp
	            Y(P) = yp
                X(Q) = xq
                Y(Q) = yq
                MoveCorners = .False.
                exit
            endif
        else
            if(QuadIsInverted(Dim,Corn,X,Y,Elm)) then
                X(P) = xp
	            Y(P) = yp
                X(Q) = xq
                Y(Q) = yq
                MoveCorners = .False.
                exit    
            endif
        endif
    end do

endif
!===========================================================================================
End Function MoveCorners
!*********************************************************************************************
