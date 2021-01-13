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
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine CalcMaxAngleDiff(Dim,NP,X,Y,Z,Point,NBoundaryFaces,Ux,Uy,Uz,Vx,Vy,Vz,MaxAngleDiff)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,I,J,NBoundaryFaces,Point
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:100)::Ux,Uy,Uz,Vx,Vy,Vz
 Real(8),Dimension(1:100)::NormalVectorLength
 Real(8)::MaxAngleDiff,IJ,CTeta,AngleDeg
 Real(8)::Vx1,Vy1,Vz1,Vx2,Vy2,Vz2,Nx,Ny,Nz
!*********************************************************************************************
!Part 1:
 Do I=1,NBoundaryFaces
    Nx = ( Uy(I)-Y(Point))*(Vz(I)-Z(Point)) - (Uz(I)-Z(Point))*(Vy(I)-Y(Point) )
    Ny = ( Uz(I)-Z(Point))*(Vx(I)-X(Point)) - (Ux(I)-X(Point))*(Vz(I)-Z(Point) )
    Nz = ( Ux(I)-X(Point))*(Vy(I)-Y(Point)) - (Uy(I)-Y(Point))*(Vx(I)-X(Point) )
    NormalVectorLength(I) = SQRT( (Nx**2)+(Ny**2)+(Nz**2) )
 EndDo
 
!Part 2:
 MaxAngleDiff = 0
 Do I=1,NBoundaryFaces
    Do J=(I+1),NBoundaryFaces
       Nx = (Uy(I)-Y(Point))*(Vz(I)-Z(Point)) - (Uz(I)-Z(Point))*(Vy(I)-Y(Point))
       Ny = (Uz(I)-Z(Point))*(Vx(I)-X(Point)) - (Ux(I)-X(Point))*(Vz(I)-Z(Point))
       Nz = (Ux(I)-X(Point))*(Vy(I)-Y(Point)) - (Uy(I)-Y(Point))*(Vx(I)-X(Point))
         
       Vx1 = Nx
       Vy1 = Ny
       Vz1 = Nz
       Nx = (Uy(J)-Y(Point))*(Vz(J)-Z(Point)) - (Uz(J)-Z(Point))*(Vy(J)-Y(Point))
       Ny = (Uz(J)-Z(Point))*(Vx(J)-X(Point)) - (Ux(J)-X(Point))*(Vz(J)-Z(Point))
       Nz = (Ux(J)-X(Point))*(Vy(J)-Y(Point)) - (Uy(J)-Y(Point))*(Vx(J)-X(Point))
        
      !Part 3:
       Vx2 = Nx
       Vy2 = Ny
       Vz2 = Nz  
       IJ    = Vx1*Vx2 + Vy1*Vy2 + Vz1*Vz2
       CTeta = IJ/(NormalVectorLength(I)*NormalVectorLength(J))
       AngleDeg = ACOS(Cteta)*57.295779513
       If(MaxAngleDiff<AngleDeg) MaxAngleDiff = AngleDeg
       
    EndDo
 EndDo
!*********************************************************************************************
 End
!###########################################################################################