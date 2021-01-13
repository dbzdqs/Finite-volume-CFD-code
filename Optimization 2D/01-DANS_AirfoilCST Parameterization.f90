!DDDDDDDDDDDDDDDDDDDDDDDDDDDAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNSSSSSSSSSSSSSSSSSSSSSSSSSS
!//             /////////////       ////////////////    ////////     //////    ////////////////        //!
!//             /////////////      ////////////////    //////////   //////    ////////////////         //!
!//            /////    /////     //////    //////    //////////// //////    /////                     //!
!//           /////    //////    ////////////////    ///////////////////    ////////////////           //!
!//          /////    //////    ////////////////    ////// ////////////               /////            //!
!//         ///////////////    //////    //////    //////   //////////    ////////////////             //!
!//       ///////////////     //////    //////    //////     ////////    ////////////////              //!
!//          Developer            Assistant    in      Numerical             Sciences                  //!
!//----------------------------------------------------------------------------------------------------//!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace Department, Amirkabir University of Technology           //!
!// Chief Developer: N. msnkre, Aerospace eng., Amirkabir University of Technology                     //!
!// Date: October, 14, 2013                                                                            //!
!//                                                                                                    //!
!// The Program is Available Through the Website: www.DANS.ir                                          //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                               //!
!//----------------------------------------------------------------------------------------------------//!
!// Description:                                                                                       //!
!// This code designed to parameterize airfoil by CST (Class-Shape-Transformation) method for          //!
!// optimization purposes.                                                                             //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_CSTAirfoilParameterization
 Implicit None
!===============================
 Integer,Parameter::Dim=100000  
!===============================

 Integer::I,J
 Integer::NP  !Number of Existing Points
 Integer::NC !Number of Cells of mesh
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
 Integer::NR   !Number Of Regions of mesh
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Integer::BSOrder_Up
 Integer::BSOrder_Lw
 Integer::NPtCurv_Up
 Integer::NPtCurv_Lw
 Integer::UpRegion
 Integer::LwRegion
 Integer::DimL
 Integer::DimU
 Real(8)::Zita_TE
 Real(8)::TE_Thick
 Real(8)::Shap 
 Real(8),Allocatable,Dimension(:,:)::BMI_Up
 Real(8),Allocatable,Dimension(:,:)::BMI_Lw
 Integer,Allocatable,Dimension(:)::IBP_Lw
 Integer,Allocatable,Dimension(:)::IBP_Up
 Real(8),Allocatable,Dimension(:)::Pol_Coeff_Up
 Real(8),Allocatable,Dimension(:)::Pol_Coeff_Lw
 Real(8),Allocatable,Dimension(:)::SBC_Up
 Real(8),Allocatable,Dimension(:)::SBC_Lw
 Real(8),Allocatable,Dimension(:)::X_Up
 Real(8),Allocatable,Dimension(:)::Y_Up
 Real(8),Allocatable,Dimension(:)::X_Lw
 Real(8),Allocatable,Dimension(:)::Y_Lw
 Real(8),Allocatable,Dimension(:)::S_Up
 Real(8),Allocatable,Dimension(:)::S_Lw
 Real(8),Allocatable,Dimension(:)::Ynew_Lw
 Real(8),Allocatable,Dimension(:)::Ynew_Up
 Real(8),Allocatable,Dimension(:)::ShapeFunc_Coe
!**************************************Main****************************************************	
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call Read_CST_Parameters(BSOrder_Up,BSOrder_Lw,UpRegion,LwRegion)
 
!Part 3:
 DimU=NFR(UpRegion)+1
 DimL=NFR(LwRegion)+1

!Part 4:
 Allocate( Pol_Coeff_Lw(0:BSOrder_Lw)            , Pol_Coeff_Up(0:BSOrder_Up),&
		   SBC_Lw(1:BSOrder_Lw+1)                , SBC_Up(1:BSOrder_Up+1),&            
           BMI_Lw(1:BSOrder_Lw+1,1:BSOrder_Lw+1) , BMI_Up(1:BSOrder_Up+1,1:BSOrder_Up+1),& 
		   IBP_Lw(1:DimL)                        , IBP_Up(1:DimU),&
		   X_Lw(1:DimL)                          , X_Up(1:DimU),&
		   Y_Lw(1:DimL)                          , Y_Up(1:DimU),&
		   S_Lw(1:DimL)                          , S_Up(1:DimU),&
		   Ynew_Lw(1:DimL)                       , Ynew_Up(1:DimU),&
           ShapeFunc_Coe(1:BSOrder_Up+1+BSOrder_Lw+1))

!Part 5:
 Call CST_PreProcc(Dim,DimU,DimL,NP,NR,NFR,IDS,UpRegion,LwRegion,X,Y,&
                   NPtCurv_Lw,NPtCurv_Up,Zita_TE,TE_Thick,IBP_Lw,IBP_Up,X_Up,Y_Up,X_Lw,Y_Lw)

!Part 6:
 Call Shape_Function(DimU,DimL,NPtCurv_Up,X_Up,Y_Up,NPtCurv_Lw,X_Lw,Y_Lw,Zita_TE,S_Up,S_Lw)

!Part 7:
 Call Bernstein_Matrix_Inverse(BSOrder_Up+1,BMI_Up)
 Call Bernstein_Matrix_Inverse(BSOrder_Lw+1,BMI_Lw)

!Part 8:
 Call Poly_Fit(NPtCurv_Up,X_Up,S_Up,BSOrder_Up,Pol_Coeff_Up)
 Call Poly_Fit(NPtCurv_Lw,X_Lw,S_Lw,BSOrder_Lw,Pol_Coeff_Lw)

!Part 9:
 Do I=1,BSOrder_Up+1

    Shap=0.0
    Do J = 1, BSOrder_Up+1
       Shap = Shap + BMI_Up(I,J)*Pol_Coeff_up(J-1)				
    End Do

  ShapeFunc_Coe(I)= Shap
 End Do

!Part 10:
 Do I=1,BSOrder_Lw+1

    Shap=0
    Do J = 1, BSOrder_Lw+1
       Shap = Shap + BMI_Lw(I,J)*Pol_Coeff_Lw(J-1)
    End Do

    ShapeFunc_Coe(I+BSOrder_Up+1)= Shap
 End Do

!Part 11:
 Call CST_InversToShap(NPtCurv_Up,NPtCurv_Lw,BSOrder_Up,BSOrder_Lw,ShapeFunc_Coe,Zita_te,X_Up,X_Lw,Ynew_Up,Ynew_Lw)

!Part 12:
 Deallocate( Pol_Coeff_Lw , Pol_Coeff_Up,SBC_Lw , SBC_Up,BMI_Lw , BMI_Up )
 
!Part 13:
 Open(10,File='ShapOut.Plt')

 Do I=1,NPtCurv_Lw
    write(10,*) X_Lw(I),Ynew_Lw(I)  
 End Do
 Do I=NPtCurv_Up,1,-1
    write(10,*) X_Up(I),Ynew_Up(I)   
 End Do

!Part 13:
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Omesh.plt")
!*********************************************************************************************
 End
!###########################################################################################

