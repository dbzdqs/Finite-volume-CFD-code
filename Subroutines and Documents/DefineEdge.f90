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
Subroutine DefineEdge(Dim,NC,NP,NK,NV,FrontIndex,ME,ANL,ANR,Corn,Neib,FrontEdges,Fronts,X,Y,States,newQuad,MFEC,NFE,NFEC,NF,Situation)
Implicit None
!===========================================================================================
Intent(In)::Dim,NK,NV,FrontIndex,ANL,ANR,Situation
Intent(Out)::newQuad
Intent(Inout)::Corn,Neib,X,Y,States,NC,NP,MFEC,ME,NFE,NFEC,NF

Integer,Parameter::Non_Front=-1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Element=4

Integer::Dim,NC,NP,NK,NV,NM,Nn,Fronts,FrontIndex,ME,MFEC,NFE,NFEC,NF,Situation,K,tri,tri_up,tri_down,I,E,V1,V2,P,Q,FirCorner,SecCorner
Integer,Dimension(1:2)::Ei
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:1000)::TElms,QElms
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::SwapIsSafe,isOnTheBoundary,NumOfEdgesInTheLoopIsEven,isFrontEdge,TopEdgeRecoveryIsPossible,CornerFound,isSideEdge
Real(8)::A,Alpha,ANL,ANR,GetAngle,GetNorm,pi,epsilon,x_value,y_value,temp,Norm_KM,Norm_MF,Norm_NF,Alpha_Left,Alpha_Right,Beta,Beta1,Beta2,Beta_Left,Beta_Right,Theta1,Theta2,a1,b1,c1,a2,b2,c2
Real(8)::Ux,Uy,Vx,Vy,norm1,norm2,norm3,BVx,BVy,value,t,Deltax,Deltay,ExtProduct
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!pi = 4*ATAN(1.0)
pi = 3.14159265358979d0
epsilon = pi/6

Select Case(Situation)

    Case(LeftVertex) !------------ Defining Edge on Left Side of current Front -------------
        
        !------------------ Part 1: Finding Side Edge from EXISTING edge -------------------
        CornerFound = .False.
        Alpha = epsilon
        Alpha_Right = 0.0
        tri = ME
        P = MFEC
        Q = NV
!Part 1:    
        do !------------- Loop for finding suitable side edge from existing edges ----------
            
            Alpha_Right = Alpha_Right + GetAngle(Dim,NK,Q,P,X,Y)
			Theta1 = Real(ANL)/2 - Alpha_Right
            
            if(DABS(Theta1) <= Alpha) then
                Alpha = DABS(Theta1)
                SecCorner = P
                E = tri
                CornerFound = .True.
            endif
            
            do I=1,3 
                if(Corn(tri,I) == Q) then
                    tri = Neib(tri,I)
                    exit
                endif   
            end do
            
            if(tri == 0) then
                exit    
            elseif(Corn(tri,4) /= 0) then
                exit    
            endif
            
            Q = P
            
            do I=1,3
                if(Corn(tri,I) /= NK .And. Corn(tri,I) /= Q) then
                    P = Corn(tri,I)
                    exit
                endif
            end do
            
        end do
        
        if(CornerFound) then
!Part 2:            
            if(.Not. isOnTheBoundary(Dim,Fronts,SecCorner,FrontEdges,States) .And. SecCorner /= newQuad(1)) then

					print *,'EXSISTING edge for 2nd Corner'
					newQuad(2)=SecCorner                  ! Saving point number of second corner of new Quad

            else !------------------------ Closing The Front Loop --------------------------

                if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,SecCorner,1) .And. NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,SecCorner,2) .And. SecCorner /= newQuad(1)) then
                    
                    print *,'EXSISTING edge for 2nd Corner'
					newQuad(2)=SecCorner                  ! Saving point number of second corner of new Quad
                    
                else
                    
                    print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                    
                    do I=1,3
						if(Corn(E,I) /= SecCorner .And. Corn(E,I)/=NK) then
							Q = Corn(E,I)
                            tri = Neib(E,I)
							exit
						endif
					end do
                   
					do I=1,3
						if(Neib(tri,I) == E) then
							NM = Corn(tri,I)
							exit
						endif
                    end do
		
                    Call CalcBisectorIntersection(Dim,X,Y,Q,NK,SecCorner,x_value,y_value)
                    Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,Q,NK,SecCorner,NM,Nn)

					newQuad(2)=Nn                       ! Saving point number of second corner of new Quad
                    
                endif
                
            endif
            
        else !--------------- Part 2: Trying SWAP or SPLIT operations ----------------------
!Part 3:               
            Ux = X(NV) - X(NK) !---------------- Defining Vector U ---------------------
            Uy = Y(NV) - Y(NK)
            
            Vx = X(NF) - X(NK) !---------------- Defining Vector V ---------------------
            Vy = Y(NF) - Y(NK)
            
            ExtProduct = Ux*Vy - Uy*Vx
            
            if(ExtProduct /= 0) then
                
                norm1 = DSQRT(Ux*Ux + Uy*Uy) !----------- Norm of U vector -----------------
                norm2 = DSQRT(Vx*Vx + Vy*Vy) !----------- Norm of V vector -----------------
            
                BVx = norm1*Vx + norm2*Ux !----------- Defining Bisector Vector ------------ 
                BVy = norm1*Vy + norm2*Uy
            
            else
                
                Ux = X(NV) - X(NK) !---------------- Defining Vector U ---------------------
                Uy = Y(NV) - Y(NK)
                
                BVx = Uy*(-1.0)    !-------------- Defining Bisector Vector ---------------- 
                BVy = Ux 
                
            endif
            
            !--------------------------- Finding Edge Eo -----------------------------------
            tri = ME
            P = MFEC
            Q = NV
    
            do
                
               Ux = X(Q) - X(P) !------------ Defining Vector U (from P to Q) --------------
               Uy = Y(Q) - Y(P) 
               
               !--- Checking Intersection of Bisector Vector (BV) and current Vector (U) ---
               !----------- Consider BV = D0 as a vector in form of (P0 + sD0) -------------
               !----------- Consider U = D1 as a vector in form of (P1 + tD1) --------------
               !----------------- According to that P0 = NK and P1 = P ---------------------
               Deltax = X(P) - X(NK)  
               Deltay = Y(P) - Y(NK)
               
               value = BVx*Uy - BVy*Ux
               
               t = (Deltax*BVy - Deltay*BVx)/value
               
               if(t > 0 .And. t < 1) then
                   E = tri
                   exit
               endif
               
                do I=1,3 
                    if(Corn(tri,I) == Q) then
                        tri = Neib(tri,I)
                        exit
                    endif   
                end do
            
                Q = P
            
                do I=1,3
                    if(Corn(tri,I) /= NK .And. Corn(tri,I) /= Q) then
                        P = Corn(tri,I)
                        exit
                    endif
                end do
            
            end do
            
            Ei(1) = P
            Ei(2) = Q
!Part 4:            
            if(.Not. isFrontEdge(Dim,Fronts,FrontEdges,States,Ei) .And. .Not. isSideEdge(newQuad,Ei)) then
            
                do I=1,3
                    if(Corn(E,I) == NK) then
                        tri = Neib(E,I)
                        exit
                    endif
                end do
                
                do I=1,3
                    if(Neib(tri,I) == E) then
                        NM = Corn(tri,I)
                        exit
                    endif
                end do
                
                Vx = X(NM) - X(NK) !----------- Defining vector V (from NK to NM) --------------
                Vy = Y(NM) - Y(NK)
                
                norm1 = DSQRT(Vx*Vx + Vy*Vy) !------ Norm of V vector (from NK to NM) ----------
                norm2 = DSQRT(BVx*BVx + BVy*BVy) !------------ Norm of BV vector ---------------
                
                temp = Vx*BVx + Vy*BVy !------------------- Inner Product ----------------------
                Beta = DACOS(temp/(norm1*norm2))
                
                norm2 = GetNorm(Dim,NK,P,X,Y) !-------- Norm of edge from NK to P --------------
                norm3 = GetNorm(Dim,NK,Q,X,Y) !-------- Norm of edge from NK to Q --------------
                
                temp = SQRT(3.0)*(norm2 + norm3)/2
                
                if(Beta < epsilon .And. norm1 < temp .And. SwapIsSafe(Dim,Corn,X,Y,tri,E)) then !--------- SWAP operation -----------
!Part 5:                    
                    if(.Not. isOnTheBoundary(Dim,Fronts,NM,FrontEdges,States)) then
                        
                        print *,'SWAP for 2nd Corner'
                        
                        Call SafeSwap(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,tri,E)
                        
                        newQuad(2)=NM                       ! Saving point number of second corner of new Quad
                    
                    else
                        
                        if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,NM,1) .And. NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,NM,2)) then 
                            
                            print *,'SWAP for 2nd Corner'
                            
                            Call SafeSwap(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,tri,E)
                            
                            newQuad(2)=NM                   ! Saving point number of second corner of new Quad
                            
                        else
                            
                            print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                            Call CalcBisectorIntersection(Dim,X,Y,NK,P,Q,x_value,y_value)
                            Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NK,P,Q,NM,Nn)
										
							newQuad(2)=Nn                   ! Saving point number of second corner of new Quad
                            
                        endif
                        
                    endif
                
                else !---------------------------- SPLIT operation -------------------------------
!Part 6:                
                    print *,'SPLIT for 2nd Corner'
                    
                    x_value = X(P) + t*Ux 
                    y_value = Y(P) + t*Uy
                    
                    Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NK,P,Q,NM,Nn)
                    
                    newQuad(2)=Nn                   ! Saving point number of second corner of new Quad
                    
                endif
                
            else !------------------- Case when edge Eo is a Front or Side Edge --------------------------
!Part 7:                
                if(isSideEdge(newQuad,Ei)) then
                    
                    do I=1,3
                        if(Corn(E,I) == NV) then
                            tri = Neib(E,I)
                            exit
                        endif
                    end do
                    
                    do I=1,3
                        if(Neib(tri,I) == E) then
                            NM = Corn(tri,I)
                            exit
                        endif
                    end do
                    
                    if(P /= NV) then
                        
                        print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                        Call CalcBisectorIntersection(Dim,X,Y,NV,P,NK,x_value,y_value)
                        Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NV,NK,P,NM,Nn)
										
						newQuad(2)=Nn                   ! Saving point number of second corner of new Quad
                        
                    elseif(Q /= NV) then
                        
                        print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                        Call CalcBisectorIntersection(Dim,X,Y,NV,Q,NK,x_value,y_value)
                        Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NV,NK,Q,NM,Nn)
										
						newQuad(2)=Nn                   ! Saving point number of second corner of new Quad
                        
                    endif
                    
                else
!Part 8:                    
                    if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,P,2)) then
                    
                        if(P /= newQuad(1)) then
                            
                            newQuad(2) = P              ! Saving point number of second corner of new Quad
                        
                        else
                            
                            do I=1,3
                                if(Corn(E,I) == Q) then
                                    tri = Neib(E,I)
                                    exit
                                endif
                            end do
                            
                            do I=1,3
                                if(Neib(tri,I) == E) then
                                    NM = Corn(tri,I)
                                    exit
                                endif
                            end do
                            
                            print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                            Call CalcBisectorIntersection(Dim,X,Y,Q,P,NK,x_value,y_value)
                            Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,Q,NK,P,NM,Nn)
										
						    newQuad(2)=Nn                   ! Saving point number of second corner of new Quad
                            
                        endif
                        
                    else
                        
                        if(Q /= newQuad(1)) then
                            
                            newQuad(2) = Q              ! Saving point number of second corner of new Quad 
                            
                        else
                            
                            do I=1,3
                                if(Corn(E,I) == P) then
                                    tri = Neib(E,I)
                                    exit
                                endif
                            end do
                            
                            do I=1,3
                                if(Neib(tri,I) == E) then
                                    NM = Corn(tri,I)
                                    exit
                                endif
                            end do
                            
                            print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                            Call CalcBisectorIntersection(Dim,X,Y,P,Q,NK,x_value,y_value)
                            Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,P,NK,Q,NM,Nn)
										
						    newQuad(2)=Nn                   ! Saving point number of second corner of new Quad
                            
                        endif
                        
                    endif
                    
                endif
               
            endif
            
        endif
    
    Case(RightVertex) !----------- Defining Edge on Right Side of current Front ------------
    
        !------------------ Part 1: Finding Side Edge from EXISTING edge -------------------
        CornerFound = .False.
        Alpha = epsilon
        Alpha_Right = 0.0
        tri = ME
        P = MFEC
        Q = NV
!Part 1:    
        do !------------- Loop for finding suitable side edge from existing edges ----------
            
            Alpha_Right = Alpha_Right + GetAngle(Dim,NK,Q,P,X,Y)
			Theta1 = Real(ANR)/2 - Alpha_Right
            
            if(DABS(Theta1) <= Alpha) then
                Alpha = DABS(Theta1)
                FirCorner = P
                E = tri
                CornerFound = .True.
            endif
            
            do I=1,3 
                if(Corn(tri,I) == Q) then
                    tri = Neib(tri,I)
                    exit
                endif   
            end do
            
            if(tri == 0) then
                exit    
            elseif(Corn(tri,4) /= 0) then
                exit    
            endif
            
            Q = P
            
            do I=1,3
                if(Corn(tri,I) /= NK .And. Corn(tri,I) /= Q) then
                    P = Corn(tri,I)
                    exit
                endif
            end do
            
        end do
        
        if(CornerFound) then
!Part 2:            
            if(.Not. isOnTheBoundary(Dim,Fronts,FirCorner,FrontEdges,States) .And. FirCorner /= newQuad(2)) then

					print *,'EXSISTING edge for first Corner'
					newQuad(1)=FirCorner                  ! Saving point number of first corner of new Quad

            else !------------------------ Closing The Front Loop --------------------------

                if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,FirCorner,1) .And. NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,FirCorner,2) .And. FirCorner /= newQuad(2)) then
                    
                    print *,'EXSISTING edge for first Corner'
					newQuad(1)=FirCorner                  ! Saving point number of first corner of new Quad
                    
                else
                    
                    print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                    
                    do I=1,3
						if(Corn(E,I) /= FirCorner .And. Corn(E,I)/=NK) then
							Q = Corn(E,I)
                            tri = Neib(E,I)
							exit
						endif
					end do
                   
					do I=1,3
						if(Neib(tri,I) == E) then
							NM = Corn(tri,I)
							exit
						endif
                    end do
		
                    Call CalcBisectorIntersection(Dim,X,Y,Q,NK,FirCorner,x_value,y_value)
                    Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,Q,NK,FirCorner,NM,Nn)

					newQuad(1)= Nn               ! Saving point number of first corner of new Quad
                    
                endif
                
            endif
            
        else !--------------- Part 2: Trying SWAP or SPLIT operations ----------------------
!Part 3:                        
                
            Ux = X(NV) - X(NK) !---------------- Defining Vector U ---------------------
            Uy = Y(NV) - Y(NK)
            
            Vx = X(NF) - X(NK) !---------------- Defining Vector V ---------------------
            Vy = Y(NF) - Y(NK)
            
            ExtProduct = Ux*Vy - Uy*Vx
            
            if(ExtProduct /= 0) then
            
                norm1 = DSQRT(Ux*Ux + Uy*Uy) !----------- Norm of U vector -----------------
                norm2 = DSQRT(Vx*Vx + Vy*Vy) !----------- Norm of V vector -----------------
            
                BVx = norm1*Vx + norm2*Ux !---------- Defining Bisector Vector ------------- 
                BVy = norm1*Vy + norm2*Uy
            
            else
            
                Ux = X(NV) - X(NK) !---------------- Defining Vector U ---------------------
                Uy = Y(NV) - Y(NK)
                
                BVx = Uy*(-1.0)    !-------------- Defining Bisector Vector ---------------- 
                BVy = Ux
            
            endif
            
            !--------------------------- Finding Edge Eo -----------------------------------
            tri = ME
            P = MFEC
            Q = NV
    
            do
                
               Ux = X(Q) - X(P) !------------ Defining Vector U (from P to Q) --------------
               Uy = Y(Q) - Y(P) 
               
               !--- Checking Intersection of Bisector Vector (BV) and current Vector (U) ---
               !----------- Consider BV = D0 as a vector in form of (P0 + sD0) -------------
               !----------- Consider U = D1 as a vector in form of (P1 + tD1) --------------
               !----------------- According to that P0 = NK and P1 = P ---------------------
               Deltax = X(P) - X(NK)  
               Deltay = Y(P) - Y(NK)
               
               value = BVx*Uy - BVy*Ux
               
               t = (Deltax*BVy - Deltay*BVx)/value
               
               if(t > 0 .And. t < 1) then
                   E = tri
                   exit
               endif
               
                do I=1,3 
                    if(Corn(tri,I) == Q) then
                        tri = Neib(tri,I)
                        exit
                    endif   
                end do
            
                Q = P
            
                do I=1,3
                    if(Corn(tri,I) /= NK .And. Corn(tri,I) /= Q) then
                        P = Corn(tri,I)
                        exit
                    endif
                end do
            
            end do
            
            Ei(1) = P
            Ei(2) = Q
!Part 4:            
            if(.Not. isFrontEdge(Dim,Fronts,FrontEdges,States,Ei) .And. .Not. isSideEdge(newQuad,Ei)) then
            
                do I=1,3
                    if(Corn(E,I) == NK) then
                        tri = Neib(E,I)
                        exit
                    endif
                end do
                
                do I=1,3
                    if(Neib(tri,I) == E) then
                        NM = Corn(tri,I)
                        exit
                    endif
                end do
                
                Vx = X(NM) - X(NK) !----------- Defining vector V (from NK to NM) --------------
                Vy = Y(NM) - Y(NK)
                
                norm1 = DSQRT(Vx*Vx + Vy*Vy) !------ Norm of V vector (from NK to NM) ----------
                norm2 = DSQRT(BVx*BVx + BVy*BVy) !------------ Norm of BV vector ---------------
                
                temp = Vx*BVx + Vy*BVy !------------------- Inner Product ----------------------
                Beta = DACOS(temp/(norm1*norm2))
                
                norm2 = GetNorm(Dim,NK,P,X,Y) !-------- Norm of edge from NK to P --------------
                norm3 = GetNorm(Dim,NK,Q,X,Y) !-------- Norm of edge from NK to Q --------------
                
                temp = SQRT(3.0)*(norm2 + norm3)/2
                
                if(Beta < epsilon .And. norm1 < temp .And. SwapIsSafe(Dim,Corn,X,Y,tri,E)) then !--------- SWAP operation -----------
!Part 5:                    
                    if(.Not. isOnTheBoundary(Dim,Fronts,NM,FrontEdges,States)) then
                        
                        print *,'SWAP for first Corner'
                        
                        Call SafeSwap(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,tri,E)
                        
                        newQuad(1)=NM                       ! Saving point number of first corner of new Quad
                    
                    else
                        
                        if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,NM,1) .And. NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,NM,2)) then 
                            
                            print *,'SWAP for first Corner'
                            
                            Call SafeSwap(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,tri,E)
                            
                            newQuad(1)=NM                       ! Saving point number of first corner of new Quad
                            
                        else
                            
                            print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                            Call CalcBisectorIntersection(Dim,X,Y,NK,P,Q,x_value,y_value)
                            Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NK,P,Q,NM,Nn)
										
							newQuad(1)=Nn                       ! Saving point number of first corner of new Quad
                            
                        endif
                        
                    endif
                
                else !---------------------------- SPLIT operation -------------------------------
!Part 6:                
                    print *,'SPLIT for 2nd Corner'
                    
                    x_value = X(P) + t*Ux 
                    y_value = Y(P) + t*Uy
                    
                    Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NK,P,Q,NM,Nn)
                    
                    newQuad(1)=Nn                       ! Saving point number of first corner of new Quad
                    
                endif
                
            else !------------------- Case when edge Eo is a Front or Side Edge --------------------------
!Part 7:                
                if(isSideEdge(newQuad,Ei)) then
                    
                    do I=1,3
                        if(Corn(E,I) == NV) then
                            tri = Neib(E,I)
                            exit
                        endif
                    end do
                    
                    do I=1,3
                        if(Neib(tri,I) == E) then
                            NM = Corn(tri,I)
                            exit
                        endif
                    end do
                    
                    if(P /= NV) then
                        
                        print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                        Call CalcBisectorIntersection(Dim,X,Y,NV,P,NK,x_value,y_value)
                        Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NV,NK,P,NM,Nn)
										
						newQuad(1)=Nn                   ! Saving point number of first corner of new Quad
                        
                    elseif(Q /= NV) then
                        
                        print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                        Call CalcBisectorIntersection(Dim,X,Y,NV,Q,NK,x_value,y_value)
                        Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,NV,NK,Q,NM,Nn)
										
						newQuad(1)=Nn                   ! Saving point number of first corner of new Quad
                        
                    endif
                    
                else
!Part 8:                    
                    if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,NK,P,1)) then
                    
                        if(P /= newQuad(2)) then
                            
                            newQuad(1) = P              ! Saving point number of second corner of new Quad
                        
                        else
                            
                            do I=1,3
                                if(Corn(E,I) == Q) then
                                    tri = Neib(E,I)
                                    exit
                                endif
                            end do
                            
                            do I=1,3
                                if(Neib(tri,I) == E) then
                                    NM = Corn(tri,I)
                                    exit
                                endif
                            end do
                            
                            print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                            Call CalcBisectorIntersection(Dim,X,Y,Q,P,NK,x_value,y_value)
                            Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,Q,NK,P,NM,Nn)
										
						    newQuad(1)=Nn                   ! Saving point number of second corner of new Quad
                            
                        endif
                        
                    else
                        
                        if(Q /= newQuad(2)) then
                            
                            newQuad(1) = Q              ! Saving point number of second corner of new Quad 
                            
                        else
                            
                            do I=1,3
                                if(Corn(E,I) == P) then
                                    tri = Neib(E,I)
                                    exit
                                endif
                            end do
                            
                            do I=1,3
                                if(Neib(tri,I) == E) then
                                    NM = Corn(tri,I)
                                    exit
                                endif
                            end do
                            
                            print *,'<<<<<<<<<<--- Closing The Front Loop --->>>>>>>>>>'
                            
                            Call CalcBisectorIntersection(Dim,X,Y,P,Q,NK,x_value,y_value)
                            Call AddNewPoint(tri,E,x_value,y_value,Neib,Corn,FrontEdges,Fronts,States,X,Y,Dim,NC,NP,P,NK,Q,NM,Nn)
										
						    newQuad(1)=Nn                   ! Saving point number of second corner of new Quad
                            
                        endif
                        
                    endif
                    
                endif
                
            endif
            
        endif
        
End Select
!===========================================================================================
End Subroutine DefineEdge
!*********************************************************************************************
