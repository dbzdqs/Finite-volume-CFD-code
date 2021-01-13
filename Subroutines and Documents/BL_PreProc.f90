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
!// Developed by: R. Amery, Mathmatical, Amirkabir university of Technology                //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine BL_PreProc(Dim,NR,NFR,IDS,NBL,Region_BL,X,Y,NCurv,NedgCurvs,NPtCurvs,NPBL,XBL,YBL,&
                       NEdgCurv,NPtCurv,EdgPt,IConectedEdg,BLPt,Dis1,OveralThick,Ratio,StrmThikRatio,Xref,DistributionType)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NR,NFR,IDS,NBL,Region_BL,X,Y
 Intent(Out  )::NCurv,NedgCurvs,NPtCurvs,NPBL,XBL,YBL,NEdgCurv,NPtCurv,EdgPt,IConectedEdg,BLPt
 Intent(InOut)::Dis1,OveralThick,Ratio,StrmThikRatio,Xref,DistributionType
 
 Integer::Dim,I,J,Sum,P1,P2,Pt1,Pt2
 Integer::NR
 Integer::NP
 Integer::NBL
 Integer::NCurv
 Integer::NedgCurvs
 Integer::NPtCurvs
 Integer::NPBL
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:100)::NEdgCurv
 Integer,Dimension(1:100)::NPtCurv
 Integer,Dimension(1:100)::Region_BL
 Integer,Dimension(1:2,1:Dim)::EdgPt
 Integer,Dimension(1:2,1:Dim)::IConectedEdg
 Integer,Dimension(1:25,1:Dim)::BLPt
 Integer,Dimension(1:Dim)::IP
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1:Dim)::XBL,YBL
 Integer,Dimension(1:100)::DistributionType
 Real(8),Dimension(1:100)::OveralThick,Dis1,Ratio,StrmThikRatio,Xref
!*********************************************************************************************
 !Part 1
 IP(:)=0
 NCurv = 0
 NPtCurvs=0
 Sum=0
 
 !Part 2
 Do I=1,NR
    IF( Region_BL(I)/=0 )Then  
     NCurv = NCurv + 1
     NPtCurv(NCurv) = 0
     
    !Part 3
     Do J=Sum+1,Sum+NFR(I)
	    P1 = IDS(3,J)
	    P2 = IDS(4,J)
       
       !Part 4
        If( IP(P1)==0 )Then
         NPtCurvs=NPtCurvs+1
         IP(P1)=NPtCurvs
         
         BLPt(1    ,NPtCurvs) = NPtCurvs
         BLPt(NBL+2,NPtCurvs) = P1
        
         XBL(NPtCurvs) = X(P1)
         YBL(NPtCurvs) = Y(P1)
         
         NPtCurv(NCurv) = NPtCurv(NCurv) + 1
        End If
        
       !Part 5
        If( IP(P2)==0 )Then
         NPtCurvs=NPtCurvs+1
         IP(P2)=NPtCurvs
                 
         BLPt(1    ,NPtCurvs) = NPtCurvs
         BLPt(NBL+2,NPtCurvs) = P2
         
         XBL(NPtCurvs) = X(P2)
         YBL(NPtCurvs) = Y(P2)
         
         NPtCurv(NCurv) = NPtCurv(NCurv) + 1
        End If      

     End Do
     
     Xref(NCurv)             = Xref(I)
     OveralThick(NCurv)      = OveralThick(I) 
     Dis1(NCurv)             = Dis1(I)
     Ratio(NCurv)            = Ratio(I)
     StrmThikRatio(NCurv)    = StrmThikRatio(I)
     DistributionType(NCurv) = DistributionType(I)
    
    EndIF
    
    Sum = Sum + NFR(I)		
 End Do

!Part 6
 NCurv = 0
 NedgCurvs = 0
 Sum=0
 Do I=1,NR
    IF( Region_BL(I)==0 )Goto 10
	NCurv = NCurv + 1
    NEdgCurv(NCurv) = NFR(I)
    
   !Part 7
    Do J=Sum+1,Sum+NFR(I)

	   Pt1 = IDS(3,J) 
	   Pt2 = IDS(4,J)	
       
       P1 = IP(Pt1)
       P2 = IP(Pt2)
       
	   NedgCurvs = NedgCurvs + 1
	   EdgPt(1,NedgCurvs) = P1 
	   EdgPt(2,NedgCurvs) = P2
	   
    End Do
		   	       
10  Sum = Sum + NFR(I)		
 End Do

 !Part 8
 IConectedEdg(:,:) = 0
 Do I=1,NedgCurvs

    P1 = EdgPt(1,I) 
	P2 = EdgPt(2,I)

    IConectedEdg(1,P2) = I
    IConectedEdg(2,P1) = I

 End Do
 NPBL = NPtCurvs
!*********************************************************************************************
 END
!###########################################################################################