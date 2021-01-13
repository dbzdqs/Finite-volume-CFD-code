!!!!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!!!!*******************************************************************************************
!!! Subroutine FaceValuToCell(NC,NFace_Cell,IFace_Cell,DirFace_Cell,FaceVal,CellVal)
!!! Implicit None
!!!!*******************************************************************************************
!!!
!!! integer::I,J,L,R,ifac,nFace
!!! real(8),dimension(:,:)::CellVal
!!! real(8),dimension(:,:)::FaceVal
!!! 
!!! Integer,Dimension(:,:)::IFace_Cell
!!! Integer,Dimension(:,:)::DirFace_Cell
!!! Integer,Dimension(:)::NFace_Cell
!!!!*******************************************************************************************
!!!
!!! do I=1,NC
!!!     
!!!    CellVal(:,i) = 0.0
!!!    nFace=NFace_Cell(i)
!!!    do j=1,nFace
!!!       iFac = IFace_Cell(j,i)
!!!       CellVal(:,i) = CellVal(:,i) + DirFace_Cell(:,iFac)*FaceVal(:,iFac)
!!!    enddo
!!!
!!! end do
!!! 
!!!!*******************************************************************************************
!!! End
!!!!###########################################################################################
