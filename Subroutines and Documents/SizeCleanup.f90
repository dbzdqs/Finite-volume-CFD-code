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
Subroutine SizeCleanup(Dim,NC,NP,NBE,BFP,Corn,Neib,X,Y,ME,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,ME,NBE,BFP
Intent(Out)::done
Intent(InOut)::NC,NP,Corn,Neib,X,Y

Integer::Dim,NP,NC,NBE,PC,A,Ai1,Ai2,B,Bi1,Bi2,C,Ci,D,Di,E,Ei,F,Fi,I,J,ME,NE,Nbc,Ncd,Nda,Nae,Nef,Nfb,getNeibour,newNode,newElement
Integer,Dimension(1:1000)::Points
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::IsBoundaryElement,SegmentJointIsPossible,QuadIsInverted,IsInTheList,done,CE,DF
Real(8)::GetNorm,normAB,normAD,normBC,normCE,normDF,temp1,temp2
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================

if(Corn(ME,4) /= 0) then
    
    if(.Not. IsBoundaryElement(Dim,Corn,NBE,BFP,ME)) then
        
        do I=1,4
!Part 1:    
            PC = 0
            done = .False.
        
            if(I /= 4) then
                A = Corn(ME,I)
                B = Corn(ME,I+1)
                Ai1 = I
                Bi1 = I+1
            else
                A = Corn(ME,I)
                B = Corn(ME,1)
                Ai1 = I
                Bi1 = 1
            endif
        
            if(.Not. IsInTheList(Points,PC,A)) then

			    PC = PC + 1
			    Points(PC) = A

            endif
    
           Call getOppoCorner(Dim,Corn,ME,B,D,Di) 
       
           if(.Not. IsInTheList(Points,PC,D)) then

			    PC = PC + 1
			    Points(PC) = D

            endif
       
           Call getOppoCorner(Dim,Corn,ME,A,C,Ci)
   
           if(.Not. IsInTheList(Points,PC,C)) then

			    PC = PC + 1
			    Points(PC) = C

           endif
       
           if(.Not. IsInTheList(Points,PC,B)) then

			    PC = PC + 1
			    Points(PC) = B

            endif
!Part 2:       
           normAB = GetNorm(Dim,A,B,X,Y)
           normAD = GetNorm(Dim,A,D,X,Y)
           normBC = GetNorm(Dim,B,C,X,Y)
   
           if(normAB > normAD .Or. normAB > normBC) then
   
               temp1 = normAB/normAD
               temp2 = normAB/normBC
       
               if(temp1 >= 2.5 .Or. temp2 >= 2.5) then
!Part 3:                
                   NE = getNeibour(Dim,Corn,Neib,A,B,ME)
           
                   if(NE /= 0) then !---------------------- if element exists ----------------
           
                       if(Corn(NE,4) /= 0) then !---------- if element is a quad -------------
               
                           if(.Not. IsBoundaryElement(Dim,Corn,NBE,BFP,NE)) then
                           
                               do J=1,4
                       
                                   if(Corn(NE,J) == A) Ai2 = J
                                   if(Corn(NE,J) == B) Bi2 = J
                           
                               end do
                       
                               Call getOppoCorner(Dim,Corn,NE,A,F,Fi)
                       
                               if(.Not. IsInTheList(Points,PC,F)) then

			                        PC = PC + 1
			                        Points(PC) = F

                                endif
                       
                               Call getOppoCorner(Dim,Corn,NE,B,E,Ei)
                       
                               if(.Not. IsInTheList(Points,PC,E)) then

			                        PC = PC + 1
			                        Points(PC) = E

                                endif
                   
                               normCE = GetNorm(Dim,C,E,X,Y) 
                               normDF = GetNorm(Dim,D,F,X,Y)
                   
                               temp1 = normAD/normCE
                               temp2 = normBC/normCE
                       
                               if(temp1 < 2.5 .And. temp1 > 0.4) then
                                   if(temp2 < 2.5 .And. temp2 > 0.4) then
                                       CE = .True.
                                   else 
                                       CE = .False.    
                                   endif
                               else
                                   CE = .False.
                               endif
                       
                               temp1 = normAD/normDF
                               temp2 = normAD/normDF
                       
                               if(temp1 < 2.5 .And. temp1 > 0.4) then
                                   if(temp2 < 2.5 .And. temp2 > 0.4) then
                                       DF = .True.
                                   else 
                                       DF = .False.    
                                   endif
                               else
                                   DF = .False.
                               endif
                       
                               Nbc = getNeibour(Dim,Corn,Neib,B,C,ME)
                               Ncd = getNeibour(Dim,Corn,Neib,C,D,ME)
                               Nda = getNeibour(Dim,Corn,Neib,D,A,ME)
                       
                               Nae = getNeibour(Dim,Corn,Neib,A,E,NE)
                               Nef = getNeibour(Dim,Corn,Neib,E,F,NE)
                               Nfb = getNeibour(Dim,Corn,Neib,F,B,NE)
                   
                               if(PC == 6) then
!Part 4:                       
                                   if(CE .And. SegmentJointIsPossible(Dim,X,Y,C,E,Points,PC)) then !----------------- Swap AB with CE ----------------
                       
                                       !----------------- set B to E in ME (ABCD to AECD)------------------
                                       Corn(ME,Bi1) = E
                                       !----------------- set A to C in NE (ABFE to CBFE)------------------
                                       Corn(NE,Ai2) = C
                           
                                       !--------------------- Modify ME neibours --------------------------
                                       Call setNeibour(Dim,Corn,Neib,ME,A,E,Nae)
                                       Call setNeibour(Dim,Corn,Neib,ME,C,E,NE)
                           
                                       !--------------------- Modify NE neibours --------------------------
                                       Call setNeibour(Dim,Corn,Neib,NE,B,C,Nbc)
                                       Call setNeibour(Dim,Corn,Neib,NE,C,E,ME)
                           
                                       !--------------------- Modifying neibours --------------------------
                                       Call setNeibour(Dim,Corn,Neib,Nae,A,E,ME)
                                       Call setNeibour(Dim,Corn,Neib,Nbc,B,C,NE)
                           
                                       if(QuadIsInverted(Dim,Corn,X,Y,ME) .Or. QuadIsInverted(Dim,Corn,X,Y,NE)) then
                           
                                           print*,'Unable to SizeCleanup: ',ME
                               
                                           Corn(ME,Bi1) = B
                                           Corn(NE,Ai2) = A
                               
                                           Call setNeibour(Dim,Corn,Neib,ME,A,B,NE)
                                           Call setNeibour(Dim,Corn,Neib,ME,B,C,Nbc)
                                           Call setNeibour(Dim,Corn,Neib,NE,A,B,ME)
                                           Call setNeibour(Dim,Corn,Neib,NE,A,E,Nae)
                               
                                           Call setNeibour(Dim,Corn,Neib,Nae,A,E,NE)
                                           Call setNeibour(Dim,Corn,Neib,Nbc,B,C,ME)
                               
                                       else
                                        
                                           done = .True.
                               
                                       endif
                           
                                       exit
                           
                                   elseif(DF .And. SegmentJointIsPossible(Dim,X,Y,D,F,Points,PC)) then !------------- Swap AB with DF ----------------
!Part 5:                       
                                       !---------------- set A to F in ME (ABCD to FBCD) ------------------
                                       Corn(ME,Ai1) = F
                                       !---------------- set B to D in NE (ABFE to ADFE) ------------------
                                       Corn(NE,Bi2) = D
                           
                                       !--------------------- Modify ME neibours --------------------------
                                       Call setNeibour(Dim,Corn,Neib,ME,F,B,Nfb)
                                       Call setNeibour(Dim,Corn,Neib,ME,D,F,NE)
                           
                                       !--------------------- Modify NE neibours --------------------------
                                       Call setNeibour(Dim,Corn,Neib,NE,D,A,Nda)
                                       Call setNeibour(Dim,Corn,Neib,NE,F,D,ME)
                           
                                       !--------------------- Modifying neibours --------------------------
                                       Call setNeibour(Dim,Corn,Neib,Nfb,F,B,ME)
                                       Call setNeibour(Dim,Corn,Neib,Nda,D,A,NE)
                           
                                       if(QuadIsInverted(Dim,Corn,X,Y,ME) .Or. QuadIsInverted(Dim,Corn,X,Y,NE)) then
                           
                                           print*,'Unable to SizeCleanup: ',ME
                               
                                           Corn(ME,Ai1) = A
                                           Corn(NE,Bi2) = B
                               
                                           Call setNeibour(Dim,Corn,Neib,ME,A,B,NE)
                                           Call setNeibour(Dim,Corn,Neib,ME,D,A,Nda)
                                           Call setNeibour(Dim,Corn,Neib,NE,A,B,ME)
                                           Call setNeibour(Dim,Corn,Neib,NE,F,B,Nfb)
                               
                                           Call setNeibour(Dim,Corn,Neib,Nfb,F,B,NE)
                                           Call setNeibour(Dim,Corn,Neib,Nda,D,A,ME)
                               
                                       else
                                           
                                           done = .True.
                               
                                       endif
                           
                                       exit
                           
                                   else !--------------------- Replace with 3 Quads -----------------------
!Part 6:                       
                                       NP = NP + 1
                                       newNode = NP
                           
                                       X(newNode) = (X(A) + X(B))/2
                                       Y(newNode) = (Y(A) + Y(B))/2
                           
                                       !if(PolygonIsConvex(Dim,PC,Points,X,Y)) then
                                       if(SegmentJointIsPossible(Dim,X,Y,D,newNode,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,E,newNode,Points,PC)) then
                           
                                           !----------------- Modifying 'A' corner of ME,NE -------------------
                                            Corn(ME,Ai1) = newNode
                                            Corn(NE,Ai2) = newNode
                           
                                           !---------------- Adding new Element based on ME -------------------
                           
                                           NC = NC + 1
                                           newElement = NC
                           
                                           Corn(newElement,Ai1) = E
                                           Corn(newElement,Bi1) = newNode
                                           Corn(newElement,Ci) = D
                                           Corn(newElement,Di) = A
                               
                                           !-------------------- Setting newElement neibours ------------------
                               
                                           Call setNeibour(Dim,Corn,Neib,newElement,A,E,Nae)
                                           Call setNeibour(Dim,Corn,Neib,newElement,D,A,Nda)
                                           Call setNeibour(Dim,Corn,Neib,newElement,newNode,E,NE)
                                           Call setNeibour(Dim,Corn,Neib,newElement,newNode,D,ME)
                               
                                           !--------------------- Modifying ME,NE neibours --------------------
                               
                                           Call setNeibour(Dim,Corn,Neib,NE,newNode,E,newElement)
                                           Call setNeibour(Dim,Corn,Neib,ME,newNode,D,newElement)
                               
                                           !------------------------ Modifying neibours -----------------------
                           
                                           Call setNeibour(Dim,Corn,Neib,Nae,A,E,newElement)
                                           Call setNeibour(Dim,Corn,Neib,Nda,D,A,newElement)
                               
                                           if(QuadIsInverted(Dim,Corn,X,Y,ME) .Or. QuadIsInverted(Dim,Corn,X,Y,NE) .Or. QuadIsInverted(Dim,Corn,X,Y,newElement)) then
                           
                                               print*,'Unable to SizeCleanup: ',ME
                                   
                                               Corn(ME,Ai1) = A
                                               Corn(NE,Ai2) = A
                                   
                                               Call setNeibour(Dim,Corn,Neib,NE,A,E,Nae)
                                               Call setNeibour(Dim,Corn,Neib,ME,D,A,Nda)
                                               Call setNeibour(Dim,Corn,Neib,Nae,A,E,NE)
                                               Call setNeibour(Dim,Corn,Neib,Nda,D,A,ME)
                                   
                                               Call DeleteCell(Dim,newElement,Corn)
                                   
                                           else
                                    
                                               done = .True.
                               
                                           endif
                               
                                           exit
                               
                                       endif
                           
                                   endif
                       
                               endif
                           
                           endif
                   
                       endif
               
                   endif
           
               endif
       
           endif
    
        end do
    
    endif

endif
!===========================================================================================
End Subroutine SizeCleanup 
!*********************************************************************************************
