















































































































































































































































      MODULE CWSTWVCP

!
!svn $Id: waves_coupler.F 756 2008-09-14 20:18:28Z jcwarner $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group      Hernan G. Arango   !
!   Licensed under a MIT/X style license                               !
!   See License_ROMS.txt                                               !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between SWAN   !
!  other coupled model(s) using the Model Coupling Toolkit (MCT).      !
!                                                                      !
!=======================================================================
!
!  Componenet model registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain decompositin descriptor datatype and assocoiated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field storage data types and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_indxR => indexRA
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communitcations scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY : MCT_isend => isend
      USE m_Transfer, ONLY : MCT_irecv => irecv
      USE m_Transfer, ONLY : MCT_waitr => waitrecv
      USE m_Transfer, ONLY : MCT_waits => waitsend
!
!
!  Sparse Matrix DataType and associated methods.
!
      USE m_SparseMatrix, ONLY : SparseMatrix
      USE m_SparseMatrix, ONLY : SparseMatrix_init => init
      USE m_SparseMatrix, ONLY : SparseMatrix_importGRowInd =>          &
     &                           importGlobalRowIndices
      USE m_SparseMatrix, ONLY : SparseMatrix_importGColInd =>          &
     &                           importGlobalColumnIndices
      USE m_SparseMatrix, ONLY : SparseMatrix_importMatrixElts =>       &
     &                           importMatrixElements
      USE m_SparseMatrix, only : SparseMatrix_lsize => lsize
      USE m_SparseMatrix, only : SparseMatrix_clean => clean
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_init => init
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_clean => clean
!
!  Decompose matrix by row.
!
      USE m_SparseMatrixPlus, ONLY : Xonly
!
!  Matrix-Vector multiply methods.
!
      USE m_MatAttrVectMul, ONLY : MCT_MatVecMul => sMatAvMult

      USE W3ODATMD, ONLY: SCREEN

      implicit none
!
!     PRIVATE

      PUBLIC :: INIT_WVCP
      PUBLIC :: initialize_wav_routers
      PUBLIC :: COAWST_CPL
      PUBLIC :: wav2atm_coupling
      PUBLIC :: wavfatm_coupling
      PUBLIC :: finalize_wav_coupling
!
!  Declarations.
!
      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapSWAN         ! GloabalSegMap variables
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)

      TYPE T_AttrVect_G
        TYPE(AttrVect) :: atm2wav_AV            ! AttrVec variables
        TYPE(AttrVect) :: wav2atm_AV            ! AttrVec variables
      END TYPE T_AttrVect_G
      TYPE (T_AttrVect_G), ALLOCATABLE :: AttrVect_G(:)


      TYPE T_GSMapInterp_A
        TYPE(GlobalSegMap) :: GSMapWRF        ! GloabalSegMap variables
      END TYPE T_GSMapInterp_A
      TYPE (T_GSMapInterp_A), ALLOCATABLE :: GSMapInterp_A(:,:)

      TYPE T_Router_A
        type(Router)   :: SWANtoWRF           ! Router variables
      END TYPE T_Router_A
      TYPE (T_Router_A), ALLOCATABLE :: Router_A(:,:)

      TYPE T_AV2_A
        TYPE(AttrVect) :: wav2atm_AV2         ! AttrVect variables
        TYPE(AttrVect) :: atm2wav_AV2 
      END TYPE T_AV2_A
      TYPE (T_AV2_A), ALLOCATABLE :: AV2_A(:,:)

      TYPE(SparseMatrix) :: sMatW             ! Sparse matrix elements
      TYPE(SparseMatrix) :: sMatA             ! Sparse matrix elements
      TYPE T_SMPlus_G
        TYPE(SparseMatrixPlus) :: A2WMatPlus    ! Sparse matrix plus elements
        TYPE(SparseMatrixPlus) :: W2AMatPlus    ! Sparse matrix plus elements
      END TYPE T_SMPlus_G
      TYPE (T_SMPlus_G), ALLOCATABLE :: SMPlus_G(:,:)
 
      CONTAINS

!***********************************************************************
!                                                                      *
      SUBROUTINE COAWST_CPL (first)
!                                                                      *
!***********************************************************************
!
      USE MCT_COUPLER_PARAMS
!
      IMPLICIT NONE
!
      INTEGER,  intent(in) :: first
      INTEGER   io, iw, ia, offset, run_couple

!
!     Call to get data from atm model.
      DO iw=1,Nwav_grids
!        CALL INIT_POINTERS(iw)
!        CALL INIT_COUPLING_POINTERS(iw)
        DO ia=1,Natm_grids
          run_couple=1
!          IF ((first.eq.1).and.(iics(iw).eq.0)) run_couple=0
          IF (MOD(first, nWAVFATM(1,1)).ne.0) run_couple=0
          IF (run_couple.eq.1) THEN
            CALL WAVFATM_COUPLING (iw, ia)
          ELSE
            GOTO 55
          END IF
        END DO
      END DO
  55  CONTINUE
!
!     Send data to atm model.
      DO iw=1,Nwav_grids
!        CALL INIT_POINTERS(iw)
!        CALL INIT_COUPLING_POINTERS(iw)
        DO ia=1,Natm_grids
          run_couple=1
!          IF ((first.eq.1).and.(iics(iw).eq.0)) run_couple=0
          IF (MOD(first, nWAV2ATM(1,1)).ne.0) run_couple=0
          IF (run_couple.eq.1) THEN
            CALL WAV2ATM_COUPLING (iw, ia)
          ELSE
            GOTO 45
          END IF
        END DO
      END DO
  45  CONTINUE

      RETURN
      END SUBROUTINE COAWST_CPL

      SUBROUTINE INIT_WVCP (ng)
!
!=======================================================================
!                                                                      !
!  Initialize waves and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE MCT_COUPLER_PARAMS
      USE W3GDATMD, ONLY: NX, NY, NSEA
!     USE mod_coupler_iounits
!
      include 'mpif.h'
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: gsmsize, Nprocs
      integer :: i, j, io, ia, Isize, Jsize, Asize
      integer :: nRows, nCols, num_sparse_elems
      integer :: cid, cad
      character (len=70)  :: nc_name
      character (len=20)  :: to_add
      character (len=120) :: wostring
      character (len=120) :: owstring

      real :: cff

!     integer, dimension(2) :: src_grid_dims, dst_grid_dims
      integer, allocatable :: start(:), length(:)
!
!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (WAV_COMM_WORLD, Nprocs, MyError)
!
!  Initialize MCT coupled model registry.
!
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Nwav_grids))
        ALLOCATE(AttrVect_G(Nwav_grids))
        ALLOCATE(SMPlus_G(Nwav_grids,Natm_grids))
        ALLOCATE(AV2_A(Nwav_grids,Natm_grids))
        ALLOCATE(GSMapInterp_A(Nwav_grids,Natm_grids))
      END IF
!
      WAVid=wavids(ng)
      IF (Nwav_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      WAV_COMM_WORLD,myids=wavids)
      ELSE
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      WAV_COMM_WORLD,WAVid)
      END IF
!
!  Initialize a Global Segment Map for non-haloed transfer of data for
!  SWAN. Determine non-haloed start and length arrays for this
!  processor. For now, this will set up a tiled exchange that is 
!  not identical to thw ww3 tiling.
!
      Jsize=INT(NY/Nprocs)
      IF (MyRank.eq.Nprocs-1) THEN
        Jsize=NY-Jsize*(Nprocs-1)
      ENDIF
      IF (.not.allocated(start)) THEN
        allocate ( start(1) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(1) )
      END IF
      start=(MyRank*INT(NY/Nprocs))*NX+1
      length=Jsize*NX

      gsmsize=length(1)
!
      CALL GlobalSegMap_init (GlobalSegMap_G(ng)%GSMapSWAN, start,      &
     &                        length, 0, WAV_COMM_WORLD, WAVid)
      deallocate (start)
      deallocate (length)

!
!  If wave grid and atm grids are different sizes, then
!  develop sparse matrices for interpolation.
!
  35  FORMAT(a3,i1,a7,i1,a11)
      DO ia=1,Natm_grids
!!!!!!!!!!!!!!!!!!!!!!
! First work on atm to wave.
!!!!!!!!!!!!!!!!!!!!!!
!
        IF (MyRank.eq.0) THEN
!          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'atm',ia,'_to_wav',ng,'_weights.nc'
!          ELSE
!            nc_name=A2Wname(ia,ng)
!          END IF
          call get_sparse_matrix (ng, nc_name, num_sparse_elems,        &
     &                            src_grid_dims, dst_grid_dims)
!
! Init the sparse matrix.
!
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          nCols=src_grid_dims(1)*src_grid_dims(2)
!
! Create sparse matrix.
!
!         Sparse rows is the dst address. Multiply the interp weights
!         by the dst masking.
!
          DO i=1,num_sparse_elems
            j=sparse_rows(i)
            cff=REAL(dst_grid_imask(j),m8)
            sparse_weights(i)=sparse_weights(i)*cff
          END DO

          call SparseMatrix_init(sMatA,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatA, sparse_rows,           &
     &                                    size(sparse_rows))
          call SparseMatrix_importGColInd(sMatA, sparse_cols,           &
     &                                    size(sparse_cols))
          call SparseMatrix_importMatrixElts(sMatA, sparse_weights,     &
     &                                       size(sparse_weights))
!
! Deallocate arrays.
!
          deallocate ( sparse_rows )
          deallocate ( sparse_cols )
          deallocate ( sparse_weights )
          deallocate ( dst_grid_imask )

!!!!!!!!!!!!!!!!!!!!!!
! Second work on waves to atm.
!!!!!!!!!!!!!!!!!!!!!!
!
!          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'wav',ng,'_to_atm',ia,'_weights.nc'
!          ELSE
!            nc_name=W2Aname(ng,ia)
!          END IF
          call get_sparse_matrix (ng, nc_name, num_sparse_elems,        &
     &                            src_grid_dims, dst_grid_dims)
!
! Init the sparse matrix.
!
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          nCols=src_grid_dims(1)*src_grid_dims(2)
!
! Create sparse matrix.
!
          DO i=1,num_sparse_elems
            j=sparse_rows(i)
            cff=REAL(dst_grid_imask(j),m8)
            sparse_weights(i)=sparse_weights(i)*cff
          END DO
!
! Load the dst grid as a coupling mask.
!
          allocate(W2A_CPLMASK(ng,ia)%dst_mask(nRows))
          DO i=1,nRows
            W2A_CPLMASK(ng,ia)%dst_mask(i)=dst_grid_imask(i)
          END DO
!
          call SparseMatrix_init(sMatW,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatW, sparse_rows,           &
     &                                    size(sparse_rows))
          call SparseMatrix_importGColInd(sMatW, sparse_cols,           &
     &                                    size(sparse_cols))
          call SparseMatrix_importMatrixElts(sMatW, sparse_weights,     &
     &                                       size(sparse_weights))
!
! Deallocate arrays.
!
          deallocate ( sparse_rows )
          deallocate ( sparse_cols )
          deallocate ( sparse_weights )
          deallocate ( dst_grid_imask )
        END IF
        CALL mpi_bcast(dst_grid_dims, 2, MPI_INTEGER, 0,                &
     &                 WAV_COMM_WORLD, MyError)
!
! scatter dst_grid_imask to be used as cpl_mask
!
        IF (MyRank.ne.0) THEN
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          allocate(W2A_CPLMASK(ng,ia)%dst_mask(nRows))
        END IF
        CALL mpi_bcast(W2A_CPLMASK(ng,ia)%dst_mask, nRows,              &
     &                 MPI_INTEGER, 0,                                  &
     &                 WAV_COMM_WORLD, MyError)
!
!  Initialize a Global Segment Map for non-haloed transfer of data
!  for the atmosphere model.
!
        Isize=INT(dst_grid_dims(1)/Nprocs)
        IF (MyRank.eq.Nprocs-1) THEN
          Isize=dst_grid_dims(1)-Isize*(Nprocs-1)
        ENDIF
        IF (.not.allocated(start)) THEN
          allocate ( start(1) )
        END IF
        IF (.not.allocated(length)) THEN
          allocate ( length(1) )
        END IF
        start=(MyRank*INT(dst_grid_dims(1)/Nprocs))*dst_grid_dims(2)+1
        length=Isize*dst_grid_dims(2)
!
        CALL GlobalSegMap_init (GSMapInterp_A(ng,ia)%GSMapWRF,          &
     &                          start, length, 0, WAV_COMM_WORLD, WAVid)
        deallocate (start)
        deallocate (length)
        call mpi_barrier(WAV_COMM_WORLD, MyError)
!
! Create ATM sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
        call SparseMatrixPlus_init(SMPlus_G(ng,ia)%A2WMatPlus, sMatA,   &
     &                             GSMapInterp_A(ng,ia)%GSMapWRF,       &
     &                             GlobalSegMap_G(ng)%GSMapSWAN,        &
     &                             Xonly,0,WAV_COMM_WORLD, WAVid)
        call SparseMatrix_clean(sMatA)
!
! Create Wave sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
         call SparseMatrixPlus_init(SMPlus_G(ng,ia)%W2AMatPlus, sMatW,  &
     &                              GlobalSegMap_G(ng)%GSMapSWAN,       &
     &                              GSMapInterp_A(ng,ia)%GSMapWRF,      &
     &                              Xonly,0,WAV_COMM_WORLD, WAVid)
        call SparseMatrix_clean(sMatW)
      END DO
!
!  Initialize attribute vector holding the export data code strings of
!  the wave model.
!
      cad=LEN(wostring)
      DO i=1,cad
        wostring(i:i)=''
      END DO
      cid=1
!
!
!  Initialize attribute vector holding the export data code string of
!  the atmosphere model.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapSWAN,            &
     &                         WAV_COMM_WORLD)
      CALL AttrVect_init (AttrVect_G(ng)%atm2wav_AV, rlist="U10:V10",   &
     &                    lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%atm2wav_AV)
!
! Initialize atribute vector holding wave data to atm.
!
      CALL AttrVect_init(AttrVect_G(ng)%wav2atm_AV,                     &
     &                   rList="HSIGN:WLENP:RTP",                       &
     &                   lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%wav2atm_AV)
!
      DO ia=1,Natm_grids
!  Initialize attribute vector holding the export data code strings of
!  the atm model. The Asize is the number of grid point on this
!  processor.
!
        Asize=GlobalSegMap_lsize(GSMapInterp_A(ng,ia)%GSMapWRF,         &
     &                           WAV_COMM_WORLD)
        CALL AttrVect_init (AV2_A(ng,ia)%atm2wav_AV2, rlist="U10:V10",  &
     &                      lsize=Asize)
        CALL AttrVect_zero (AV2_A(ng,ia)%atm2wav_AV2)
!
!  Initialize attribute vector holding the export data code string of
!  the wave model.
!
        CALL AttrVect_init(AV2_A(ng,ia)%wav2atm_AV2,                    &
     &                     rList="HSIGN:WLENP:RTP:CPL_MASK",            &
     &                     lsize=Asize)
        CALL AttrVect_zero(AV2_A(ng,ia)%wav2atm_AV2)
      END DO

      RETURN
      END SUBROUTINE INIT_WVCP

      SUBROUTINE INITIALIZE_WAV_ROUTERS
!
!=======================================================================
!                                                                      !
!  Initialize waves routers for wave model.                            !
!                                                                      !
!=======================================================================
!
      USE MCT_COUPLER_PARAMS
!      USE ww3_iounits
!     USE M_PARALL
!
!      include 'mpif.h'
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: ng, iw, ia
!
!  Initialize MCT Routers.
!
      ALLOCATE(Router_A(Nwav_grids,Natm_grids))
!
!  Initialize a router to the atmosphere model component.
!
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          ATMid=atmids(ia)
          CALL Router_init (ATMid, GSMapInterp_A(iw,ia)%GSMapWRF,       &
     &                      WAV_COMM_WORLD, Router_A(iw,ia)%SWANtoWRF)
        END DO
      END DO

      RETURN
      END SUBROUTINE INITIALIZE_WAV_ROUTERS

      SUBROUTINE WAV2ATM_COUPLING (iw, ia)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  atm and wave models. Currently, the following data streams are      !
!  processed:                                                          !
!                                                                      !
!  Fields exported to the ATM model:                                   !
!                                                                      !
!     * Significant wave height (m)                                    !
!     * Surface wave relative peak period (s)                          !
!     * Surface wave length (m)                                        !
!                                                                      !
!  Fields imported from the ATM Model:                                 !
!                                                                      !
!     * Wind Speed (m/s)                                               !
!                                                                      !
!=======================================================================
!
      USE CONSTANTS, ONLY: PI
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF
      USE MCT_COUPLER_PARAMS
      USE W3ADATMD, ONLY: HS, PHIBBL, PHIOC, FP0, T0M1, UBA
      USE W3ADATMD, ONLY: WLM
!
      implicit none
!/MPI      INCLUDE "mpif.h"
!
!  Imported variable declarations.
!
      integer :: Numcouple, iw, ia
      integer :: IP, IX, IY
!
!  Local variable declarations.
!
      integer :: MyStatus, MyError, MySize, MyRank
      integer :: i, id, j, gsmsize, ierr, indx, Tag
      integer :: Istr, Iend, Jstr, Jend, Asize
      integer :: Isize, Jsize, INDXG, Nprocs, OFFSET
      integer :: start, length, grdsize
      integer, pointer :: points(:)

      real, pointer :: SND_BUF(:), RCV_BUF(:)
      real(m8) :: fac
      real(m8), pointer :: avdata(:)
      integer, pointer :: indices(:)
      real(m8), pointer :: Amask(:)
!
!-----------------------------------------------------------------------
!  Send wave fields to WRF.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, Nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(iw)%GSMapSWAN,          &
     &                           WAV_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      avdata=0.0_m8
!
      grdsize=NX*NY
      allocate ( SND_BUF(grdsize),stat=ierr )
      SND_BUF=0.0
      allocate ( RCV_BUF(grdsize),stat=ierr )
      RCV_BUF=0.0
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GlobalSegMap_G(iw)%GSMapSWAN,          &
     &                           MyRank, points)
!
!  Determine grid tiling for exchanges.
!
      Jsize=INT(NY/Nprocs)
      IF (MyRank.eq.Nprocs-1) THEN
        Jsize=NY-Jsize*(Nprocs-1)
      ENDIF
      start=(MyRank*INT(NY/Nprocs))*NX+1
      length=Jsize*NX
!
!  Load WW3 data into MCT storage buffers.
!  The data is exported using WRF definition for real kind m8=r8.
!-------------------------------------------------------------------
!  HS: Signfiicant wave height
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=HS(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2atm_AV,             &
     &                             "HSIGN",avdata)
!-------------------------------------------------------------------
!  RTP: Peak surface period
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=FP0(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      fac=2.0_m8*PI
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=fac/MAX(REAL(RCV_BUF(i),m8),0.001_m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2atm_AV,             &
     &                             "RTP",avdata)
!-------------------------------------------------------------------
!  WLEN: mean wave length
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=WLM(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2atm_AV,             &
     &                             "WLENP",avdata)
!-------------------------------------------------------------------
!  Send wave fields bundle to atm model, WRF.
!-----------------------------------------------------------------------
!
 35     FORMAT (a14,i2,a23,i5)
!
!  Send fields to atmosphere model.
!
        Tag=ia*10+iw
        CALL MCT_MatVecMul(AttrVect_G(iw)%wav2atm_AV,                   &
     &                     SMPlus_G(iw,ia)%W2AMatPlus,                  &
     &                     AV2_A(iw,ia)%wav2atm_AV2)
!
!  Now add in the CPL_MASK before we send it over to wrf.
!  Get the number of grid points on this processor.
!
        Asize=GlobalSegMap_lsize (GSMapInterp_A(iw,ia)%GSMapWRF,        &
     &                            WAV_COMM_WORLD)
        allocate (Amask(Asize))
        Amask=0.0
!
!  Ask for points in this tile.
!
        CALL GlobalSegMap_Ordpnts (GSMapInterp_A(iw,ia)%GSMapWRF,       &
     &                             MyRank, points)
!
!  Load the dst grid cpl mask into the attr vect.
!
        DO i=1,Asize
          Amask(i)=REAL(W2A_CPLMASK(iw,ia)%dst_mask(points(i)),m8)
        END DO
        deallocate (points)
        CALL AttrVect_importRAttr (AV2_A(iw,ia)%wav2atm_AV2, "CPL_MASK",&
     &                             Amask, Asize)
        CALL MCT_isend (AV2_A(iw,ia)%wav2atm_AV2,                       &
     &                  Router_A(iw,ia)%SWANtoWRF, Tag)
        CALL MCT_waits (Router_A(iw,ia)%SWANtoWRF)
        IF (MyRank.EQ.0) THEN
          WRITE (SCREEN,35) '== WW3 grid ',iw,' sent data to WRF grid ' &
     &                      ,ia
        END IF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupl fail swancplr, MyStatus= ', MyError
          CALL FINALIZE_WAV_COUPLING(iw)
        END IF
        deallocate (Amask)
        if (associated (indices)) then
          deallocate (indices)
        endif
      deallocate (avdata)
      deallocate (SND_BUF, RCV_BUF)
!
      RETURN
      END SUBROUTINE WAV2ATM_COUPLING
      SUBROUTINE WAVFATM_COUPLING (iw, ia)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!                                                                      !
!  Fields imported from the ATM Model:                                 !
!                                                                      !
!     * Wind E and N (m/s)                                             !
!                                                                      !
!=======================================================================
!
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF
      USE W3SERVMD
      USE W3IDATMD
      USE MCT_COUPLER_PARAMS
      USE W3ADATMD, ONLY: UA, UD
!
      implicit none
!/MPI      INCLUDE "mpif.h"
!
!  Imported variable declarations.
!
      integer :: iw, ia
!
!  Local variable declarations.
!
      integer :: IP, IX, IY, grdsize
      integer :: MyStatus, MyError, MySize, MyRank, Nprocs
      integer :: i, id, j, gsmsize, ierr, indx, Tag
      integer :: Istr, Iend, Jstr, Jend, start, length
      integer :: Isize, Jsize, INDXG, OFFSET, IAPROC
      integer, pointer :: indices(:)
      real :: cff, cffmin, cffmax
      real, parameter    ::  Large = 1.0E+20
      real, dimension(2) :: range
      real, pointer      :: SND_BUF(:), RCV_BUF(:)
      real, pointer      :: WND_U10(:), WND_V10(:)
      real(m8), pointer  :: avdata(:)
!
!-----------------------------------------------------------------------
!  Get wind data from atm.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, Nprocs, MyError)
      IAPROC=MyRank+1
      grdsize=NX*NY
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(iw)%GSMapSWAN,          &
     &                           WAV_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      avdata=0.0_m8
      MyError=0
      allocate ( SND_BUF(grdsize),stat=ierr )
      SND_BUF=0.
      allocate ( RCV_BUF(grdsize),stat=ierr )
      RCV_BUF=0.
      allocate ( WND_U10(grdsize),stat=ierr )
      WND_U10=0.
      allocate ( WND_V10(grdsize),stat=ierr )
      WND_V10=0.
!
!-----------------------------------------------------------------------
!  RCV the data from WRF.
!-----------------------------------------------------------------------
!
 35   FORMAT (a14,i2,a24,i2)
!
!  Receive fields from atmosphere model.
!
      Tag=0*100+ia*10+iw
      CALL MCT_irecv (AV2_A(iw,ia)%atm2wav_AV2,                         &
     &                Router_A(iw,ia)%SWANtoWRF, Tag)
!     Wait to make sure the WRF data has arrived.
      CALL MCT_waitr (AV2_A(iw,ia)%atm2wav_AV2,                         &
     &                Router_A(iw,ia)%SWANtoWRF)
      CALL MCT_MatVecMul(AV2_A(iw,ia)%atm2wav_AV2,                      &
     &                   SMPlus_G(iw,ia)%A2WMatPlus,                    &
     &                   AttrVect_G(iw)%atm2wav_AV)
      IF (MyRank.EQ.0) THEN
        WRITE (SCREEN,35)'== WW3 grid ',iw,' recv data from WRF grid'   &
     &                     ,ia
      END IF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling fail swancplr, MyStatus= ', MyError
        CALL FINALIZE_WAV_COUPLING(iw)
      END IF
!
! Compute local non-halo data size.
!
      Jsize=INT(NY/Nprocs)
      IF (MyRank.eq.Nprocs-1) THEN
        Jsize=NY-Jsize*(Nprocs-1)
      ENDIF
      start=(MyRank*INT(NY/Nprocs))*NX+1
      length=Jsize*NX
!
!
!  U-wind.
!
      CALL AttrVect_exportRAttr(AttrVect_G(iw)%atm2wav_AV,"U10",       &
     &                          avdata,gsmsize)
      range(1)= Large
      range(2)=-Large
      SND_BUF=0.0
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        range(1)=MIN(range(1),REAL(avdata(IP)))
        range(2)=MAX(range(2),REAL(avdata(IP)))
        SND_BUF(i)=REAL(avdata(IP))
      END DO
      CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,                 &
                         MPI_SUM, WAV_COMM_WORLD, MyError)
      CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,                 &
                         MPI_SUM, WAV_COMM_WORLD, MyError)
      IF (MyRank.eq.0) THEN
        write(SCREEN,40) 'WRFtoWW3 Min/Max U10     (ms-1):     ',       &
     &                      cffmin, cffmax
      END IF
!
!  now scatter data to all nodes.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array WND_U10 as temporary for now.
!
      DO i=1,grdsize
        IF (ia.eq.1) THEN
          WND_U10(i)=RCV_BUF(i)
        ELSE
          WND_U10(i)=WND_U10(i)+RCV_BUF(i)
        END IF
      END DO
!     DO i=1,NSEA
!       IX     = MAPSF(i,1)
!       IY     = MAPSF(i,2)
!       IP=(IY-1)*NX+IX
!       WND_SPD(i)=RCV_BUF(IP)
!     END DO
!
!  V-wind.
!
      CALL AttrVect_exportRAttr(AttrVect_G(iw)%atm2wav_AV,"V10",        &
     &                          avdata,gsmsize)
      range(1)= Large
      range(2)=-Large
      SND_BUF=0.0
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        range(1)=MIN(range(1),REAL(avdata(IP)))
        range(2)=MAX(range(2),REAL(avdata(IP)))
       SND_BUF(i)=REAL(avdata(IP))
      END DO
      CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,                 &
                         MPI_SUM, WAV_COMM_WORLD, MyError)
      CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,                 &
                         MPI_SUM, WAV_COMM_WORLD, MyError)
      IF (MyRank.eq.0) THEN
        write(SCREEN,40) 'WRFtoWW3 Min/Max V10     (ms-1):     ',       &
     &                    cffmin, cffmax
      END IF
!
!  now scatter data to all nodes.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array WND_V10 as temporary for now.
!
      DO i=1,grdsize
        IF (ia.eq.1) THEN
          WND_V10(i)=RCV_BUF(i)
        ELSE
          WND_V10(i)=WND_V10(i)+RCV_BUF(i)
        END IF
      END DO
!
!  Now we need to combine wnd speed and dir and scatter.
!
      DO i=1,NSEA
        IX     = MAPSF(i,1)
        IY     = MAPSF(i,2)
        IP=(IY-1)*NX+IX
        IF (ia.eq.1) THEN
          UA(i)=SQRT(WND_U10(IP)**2+WND_V10(IP)**2+0.000001)
          UD(i)=ATAN2(WND_V10(IP),WND_U10(IP))
        ELSE
          UA(i)=UA(i)+SQRT(WND_U10(IP)**2+WND_V10(IP)**2+0.000001)
          UD(i)=UD(i)+ATAN2(WND_V10(IP),WND_U10(IP))
        END IF
      END DO
 40     FORMAT (a36,1x,2(1pe14.6))
!
      deallocate (avdata)
      deallocate (SND_BUF, RCV_BUF)
      deallocate (WND_U10, WND_V10)
      if (associated (indices)) then
        deallocate (indices)
      endif
!
      RETURN
      END SUBROUTINE WAVFATM_COUPLING

      SUBROUTINE FINALIZE_WAV_COUPLING(iw)
!
!=======================================================================
!                                                                    ===
!  This routines terminates execution during coupling error.         ===
!                                                                    ===
!=======================================================================
      USE MCT_COUPLER_PARAMS
!
!  Local variable declarations.
!
      integer :: iw, io, ia, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      DO ia=1,Natm_grids
        CALL Router_clean (Router_A(iw,ia)%SWANtoWRF, MyError)
      END DO
      CALL AttrVect_clean (AttrVect_G(iw)%atm2wav_AV, MyError)
      CALL GlobalSegMap_clean (GlobalSegMap_G(iw)%GSMapSWAN, MyError)

      END SUBROUTINE FINALIZE_WAV_COUPLING
      END MODULE CWSTWVCP
