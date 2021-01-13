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
!// Date: June, 10, 2017                                                                   //!
!// Developed by: H. Morad Tabrizi, Mechanical Eng., Amirkabir University of Technology    //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine CST_PreProcc(Dim,DimU,DimL,NP,NR,NFR,IDS,UpRegion,LwRegion,X,Y,NPtCurv_Lw,NPtCurv_Up,&
                         Zita_TE,TE_Thick,IBP_Lw,IBP_uP,X_Up,Y_Up,X_Lw,Y_Lw)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,DimU,DimL,NP,NR,NFR,IDS,UpRegion,LwRegion
 Intent(InOut)::X,Y
 Intent(Out  )::NPtCurv_Lw,NPtCurv_Up,Zita_TE,TE_Thick,IBP_Lw,IBP_uP,X_Up,Y_Up,X_Lw,Y_Lw

 Integer::Dim,DimU,DimL,I,J,NP,NR,SumF,NF_Up1,NF_Up2,NF_Lw1,NF_Lw2,UpRegion,LwRegion,P1,P2,&
          NPtCurv_Lw,NPtCurv_Up,LE_Point
 Real(8)::Zita_TE,TE_Thick,Xo,Yo
 Integer,Dimension(1:DimL)::IBP_Lw
 Integer,Dimension(1:DimU)::IBP_Up
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1:DimU)::X_Up,Y_Up
 Real(8),Dimension(1:DimL)::X_Lw,Y_Lw
!*********************************************************************************************
!Part 1:
 SumF = 0
 Do J=1,NR
    IF(J==UpRegion)Then
	 NF_Up1 = SumF + 1
     NF_Up2 = SumF + NFR(UpRegion)
    EndIF

	SumF = SumF + NFR(J)
 End do

!Part 2:
 SumF = 0
 Do J=1,NR
    IF(J==LwRegion)Then
	 NF_Lw1 = SumF + 1
     NF_Lw2 = SumF + NFR(LwRegion)
    EndIF

	SumF = SumF + NFR(J)
 End do


!Part 3:
 NPtCurv_Lw = 0
 Do I=NF_Lw2,NF_Lw1,-1

    P1=IDS(3,I)
	P2=IDS(4,I)
		  
	NPtCurv_Lw = NPtCurv_Lw + 1
	IBP_Lw(NPtCurv_Lw) = P2

   
 End do
       	  
 NPtCurv_Lw = NPtCurv_Lw + 1
 IBP_Lw(NPtCurv_Lw) = P1

!Part 4:
 NPtCurv_Up = 0
 Do I=NF_Up1,NF_Up2
       
    P1=IDS(3,I)
    P2=IDS(4,I)
		  
	NPtCurv_Up = NPtCurv_Up + 1
	IBP_UP(NPtCurv_Up) = P1

 End do
	  
 NPtCurv_Up = NPtCurv_Up + 1
 IBP_UP(NPtCurv_Up) = P2

!Part 5:
 LE_Point = IBP_UP(1)
 Xo = X(LE_Point)
 Yo = Y(LE_Point)

 Do I=1,NP
   X(I) = X(I) - Xo
   Y(I) = Y(I) - Yo 
 End Do 

!Part 6:
 Do I=1,NPtCurv_Up
    P1 = IBP_UP(I)
    
	X_up(I) = X(P1)
    Y_up(I) = Y(P1) 
 End Do

!Part 6:
 Do I=1,NPtCurv_Lw
    P1 = IBP_Lw(I)
    
	X_Lw(I) = X(P1)
    Y_Lw(I) = Y(P1) 
 End Do

!Part 7:
 TE_Thick = Y_up(NPtCurv_Up) - Y_Lw(NPtCurv_Lw)
 Zita_te = 0.5*TE_Thick

!*********************************************************************************************
 End 
!###########################################################################################
