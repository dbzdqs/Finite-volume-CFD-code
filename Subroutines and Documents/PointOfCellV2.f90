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
!// Supervisor:                                                                            //!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Date: Feb., 10, 2018                                                                   //!
!// Supervisor:                                                                            //!
!// Developed by: A. Hemati zadeh, Tabriz University                                       //!
!// Developed by: K. Safari, Tehran University                                             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!********************************************************************************************* 
  Subroutine PointOfCellV2(Dim,NF,NC,IDS,NEdgeOfCell,IEdgeOfCell,NPoint_Cell,IPoint_Cell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF,IDS,NEdgeOfCell,IEdgeOfCell
 Intent(Out  )::NPoint_Cell,IPoint_Cell
 
 Integer::Dim,NC,NF,I,J,K,S,Point,Face,Repeated
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:4,1:Dim)::IEdgeOfCell
 Integer,Dimension(1:4,1:Dim)::IPoint_Cell
 Integer,Dimension(1:Dim)::NEdgeOfCell,NPoint_Cell
!*********************************************************************************************
!Part 1:
 NPoint_Cell(1:NC)   = 0
 IPoint_Cell(:,1:NC) = 0
 
!Part 2:
 Do I=1,NC 
    Do J=1,NEdgeOfCell(I)
       Face = IEdgeOfCell(J,I)
       
      !Part 3:
       Do S=3,4
          Point = IDS(S,Face)
             
         !Part 4:
          repeated=.FALSE.
          Do K=1,NPoint_Cell(I)
             if(IPoint_Cell(K,I)==Point)Then
              repeated = .TRUE.
              exit
             Endif
          EndDo
          If(.NOT. repeated)Then
           NPoint_Cell(I)  =  NPoint_Cell(I) + 1
           IPoint_Cell(NPoint_Cell(I),I) = Point
          Endif
          
       EndDo
    EndDo

 EndDo
!*********************************************************************************************
 End
!###########################################################################################
