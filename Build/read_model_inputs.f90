      SUBROUTINE read_model_inputs
!
!=======================================================================
!                                                                      !
!  This routine reads in model input parameters of dt and              !
!  number of grids for each model.                                     !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE ww3_iounits
      USE mod_coupler_iounits
      implicit none
!
      include 'mpif.h'
!
!  Imported variable declarations.
!
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, isval, iw, ia, inp, out, status
      integer :: MyRank, MyError, MyMaster, DT, num, den
      integer :: cdecode_line, cload_i, cload_r, indx, indx2, test
      integer :: Ivalue
      integer :: sstupdate
!     integer, allocatable :: parentid(:)
      real(m8), dimension(100) :: Rval
      real(m8) :: FAC
      character (len=1 ), parameter :: blank = ' '
      character (len=1 ) :: KEY
      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160) :: aline
      character (len=160) :: saveline1, saveline2, saveline3
      character (len=160), dimension(100) :: Cval
!
      inp=1
      out=stdout
      MyMaster=0
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
!
!     Read WW3 input file
!
      iw=1
      OPEN (inp, FILE=TRIM(Wname(iw)), FORM='formatted', STATUS='old',  &
     &      ERR=110)
      GO TO 130
 110  WRITE (out,120) Wname(iw)
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Wname(iw))
      RETURN
 120  FORMAT (/,' READ MODEL INPUTS - Unable to open ww3 input file.',  &
     &        /,a80)
 130  CONTINUE
!
!  set number of WW3 grids to = 1 for now.
!
      Nwav_grids=1
      allocate (dtwav(Nwav_grids))
!
!  now read the ww3_grid.inp file to get the dt
!
      REWIND (inp)
      test=1
      DO WHILE (test.eq.1)
! Name of grid :
        READ (inp,'(a)',ERR=116,END=140) line
        IF (line(1:1).ne.'$') test=0
      END DO
!
      test=1
      DO WHILE (test.eq.1)
! Define spectrum
        READ (inp,'(a)',ERR=116,END=140) line
        IF (line(1:1).ne.'$') test=0
      END DO
!
      test=1
      DO WHILE (test.eq.1)
! Define model run flags
        READ (inp,'(a)',ERR=116,END=140) line
        IF (line(1:1).ne.'$') test=0
      END DO
!
      test=1
      DO WHILE (test.eq.1)
! Define model time steps
        READ (inp,'(a)',ERR=116,END=140) line
        IF (line(1:1).ne.'$') test=0
      END DO
      aline=ADJUSTL(line)
      indx=INDEX(aline,blank)
      aline=aline(1:indx-1)
      read(aline,*) dtwav(iw)
! Set minimum 1 sec dt
      dtwav(iw)=MAX(dtwav(iw),1.0_m8)
!
!     write(*,*) 'ww3 dt is ', dtwav(iw)
!
      test=1
      DO WHILE (test.eq.1)
! Read rest of file
        READ (inp,'(a)',ERR=116,END=140) line
      END DO
 116  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 140  CLOSE (inp)
!
!
!     Read WRF input file
!
      OPEN (inp, FILE=TRIM(Aname), FORM='formatted', STATUS='old',      &
     &      ERR=210)
      GO TO 230
 210  WRITE (out,220) Aname
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Aname)
!     exit_flag=4
      RETURN
 220  FORMAT (/,' READ MODEL INPUTS - Unable to open wrf input file.',  &
     &        /,a80)
 230  CONTINUE
!
      sstupdate=0
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=215,END=240) line
        aline=ADJUSTL(line)
        IF(aline(1:10).eq.'time_step ') THEN
          saveline1=aline
        END IF
        IF(aline(1:19).eq.'time_step_fract_num') THEN
          saveline2=aline
        END IF
        IF(aline(1:19).eq.'time_step_fract_den') THEN
          saveline3=aline
        END IF
        IF(aline(1:7).eq.'max_dom') THEN
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          indx=INDEX(aline,',')
          read(aline(1:indx-1),'(i5)') Natm_grids
          allocate (dtatm(Natm_grids))
          allocate (parentid(Natm_grids))
          allocate (wrf_e_we(Natm_grids))
!
! Process the time steps
!
!  get DT
          aline=saveline1
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          indx=INDEX(aline,',')
          read(aline(1:indx-1),'(i5)') DT
!  get num
          aline=saveline2
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          indx=INDEX(aline,',')
          read(aline(1:indx-1),'(i5)') num
!  get den
          aline=saveline3
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          indx=INDEX(aline,',')
          read(aline(1:indx-1),'(i5)') den
!  compute dt = DT + num/den
          IF (den.eq.0) THEN
            dtatm(1)=REAL(DT,m8)
          ELSE
            dtatm(1)=REAL(DT,m8)+REAL(num,m8)/REAL(den,m8)
          END IF
        END IF
        IF(aline(1:9).eq.'parent_id') THEN
          indx=INDEX(aline,'=')
          DO ia=1,Natm_grids
            aline=aline(indx+1:LEN(aline))
            indx=INDEX(aline,',')
            saveline1=TRIM(ADJUSTL(aline(1:indx-1)))
            read(saveline1,'(i5)') parentid(ia)
            IF (parentid(ia).EQ.0) parentid(ia) = 1
          END DO
        END IF
        IF(aline(1:22).eq.'parent_time_step_ratio') THEN
          indx=INDEX(aline,'=')
          DO ia=1,Natm_grids
            aline=aline(indx+1:LEN(aline))
            indx=INDEX(aline,',')
            saveline1=TRIM(ADJUSTL(aline(1:indx-1)))
            read(saveline1,'(i5)') den
!           dtatm(ia)=dtatm(1)/REAL(den,m8)
            dtatm(ia)=dtatm(parentid(ia))/REAL(den,m8)
          END DO
        END IF
        IF(aline(1:10).eq.'sst_update') THEN
          indx=INDEX(aline,'=')
          aline=ADJUSTL(aline(indx+1:LEN(aline)))
          indx=MAX(INDEX(aline,','),LEN(aline))
          read(aline(1:indx-1),'(i1)') sstupdate
        END IF
        IF(aline(1:4).eq.'e_we') THEN
          indx=INDEX(aline,'=')
          DO ia=1,Natm_grids
            aline=aline(indx+1:LEN(aline))
            indx=INDEX(aline,',')
            saveline1=TRIM(ADJUSTL(aline(1:indx-1)))
            read(saveline1,'(i5)') wrf_e_we(ia)
          END DO
        END IF
      END DO
!
 215  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 240  CLOSE (inp)
      IF (sstupdate.eq.0) THEN
        WRITE (stdout,65) sstupdate
 65     FORMAT (/,' Recommend that sst_update be set to = 1 in the '    &
     &            'namelist.input for model coupling, not = ',i5)
!       STOP
      END IF
  60  FORMAT (/,'read model inputs - Error while processing line: ',/,a)
      RETURN
      END SUBROUTINE read_model_inputs
