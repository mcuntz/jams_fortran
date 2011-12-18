MODULE mo_nml

  ! Adapted from Echam5, (C) MPI-MET, Hamburg, Germany

  ! Reading and positioning namelist file
  ! Author:
  !     L. Kornblueh, MPI, March 2001, original source

  ! Modified Jan 2011, Matthias Cuntz - compatible with gfortran <= version 4.3
  !                                     all integer(i4)
  !                                     quiet

  USE mo_kind,         ONLY: i4
  USE mo_string_utils, ONLY: tolower
  USE mo_message,      ONLY: message, message_text
  USE mo_finish,       ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: open_nml                                      ! open namelist file
  PUBLIC :: close_nml                                     ! close namelist file
  PUBLIC :: position_nml                                  ! position namelist file
  PUBLIC :: nunitnml                                      ! namelist unit
  PUBLIC :: POSITIONED, MISSING, LENGTH_ERROR, READ_ERROR ! return values from position_nml

  ! return values of function 'position_nml'
  INTEGER(i4), PARAMETER :: POSITIONED   =  0 ! file pointer set to namelist group
  INTEGER(i4), PARAMETER :: MISSING      =  1 ! namelist group is missing
  INTEGER(i4), PARAMETER :: LENGTH_ERROR =  2 !
  INTEGER(i4), PARAMETER :: READ_ERROR   =  3 !
  ! default namelist unit
  INTEGER,     SAVE      :: nunitnml         = -1

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         open_nml

  !     PURPOSE
  !         Open a namelist file.

  !     CALLING SEQUENCE
  !         call open_nml(file, unit, quiet=quiet)
  
  !     INDENT(IN)
  !         character(len=*) :: file         Namelist filename
  !         integer          :: unit         namelist unit

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         logical :: quiet                 Be verbose or not (default: .false.)
  !                                          .true.:  no messages
  !                                          .false.: write out messages

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         call open_nml('namelist.txt',nnml)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5, (C) MPI-MET, Hamburg, Germany

  SUBROUTINE open_nml(file, unit, quiet)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: file
    INTEGER         , INTENT(IN) :: unit
    INTEGER         , INTENT(IN), OPTIONAL :: quiet
    INTEGER :: istat

    nunitnml = unit
    if (.not. present(quiet)) then
       !message_text = '    This is namelist '//file
       write(message_text,'(A,A)')'    This is namelist ', trim(file)
       CALL message(trim(message_text))
    end if
    OPEN (nunitnml, file=file, iostat=istat, status='old', action='read', &
         delim='apostrophe')

    IF (istat /= 0) THEN
       !message_text = 'Could not open namelist file '//TRIM(file)
       write(message_text,'(A,A)')'Could not open namelist file ', trim(file)
      CALL finish('OPEN_NML',trim(message_text))
    END IF

  END SUBROUTINE open_nml  

  ! ------------------------------------------------------------------

  !     NAME
  !         close_nml

  !     PURPOSE
  !         Closes the namelist file.

  !     CALLING SEQUENCE
  !         call close_nml()
  
  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         open_nml remembers the namelist unit in the public, save variable nunitnml.
  !         close_nml uses nunitnml.

  !     EXAMPLE
  !         call close_nml()
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5, (C) MPI-MET, Hamburg, Germany

  SUBROUTINE close_nml()

    IMPLICIT NONE

    INTEGER :: istat

    IF (nunitnml == -1) THEN
      CALL finish('CLOSE_NML','No namelist file opened.')
    END IF

    CLOSE(nunitnml, IOSTAT=istat)

    IF (istat /= 0) THEN
      CALL finish('CLOSE_NML','Could not close namelist file.')
    END IF

    nunitnml = -1

  END SUBROUTINE close_nml

  ! ------------------------------------------------------------------

  !     NAME
  !         position_nml

  !     PURPOSE
  !         Position namelist file pointer for reading a new namelist next.
  !         It positions namelist file at correct place for reading
  !         namelist /name/ (case independent).

  !     CALLING SEQUENCE
  !         call position_nml(name, unit=unit, first=first, status=status)
  
  !     INDENT(IN)
  !         character(len=*) :: name         namelist name (case independent)

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         integer     :: unit               namelist unit (default: nunitnml)
  !         logical     :: first              start search at beginning,
  !                                           i.e. rewind the namelist first (default: .true.)
  !                                          .true.:  rewind
  !                                          .false.: continue search from current file pointer

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         integer(i4) :: status            Set on output to either of
  !                                          POSITIONED (0)   - correct
  !                                          MISSING (1)      - name not found
  !                                          LENGTH_ERROR (2) - namelist length longer then 256 characters
  !                                          READ_ERROR (3)   - error while reading namelist file

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         call position_nml('myname',nnml)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5, (C) MPI-MET, Hamburg, Germany

  SUBROUTINE position_nml(name, unit, first, status)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in)            :: name   ! namelist group name
    INTEGER,          INTENT(in)  ,OPTIONAL :: unit   ! file unit number
    LOGICAL,          INTENT(in)  ,OPTIONAL :: first ! default: true
    INTEGER(i4),      INTENT(out) ,OPTIONAL :: status ! error return value

    CHARACTER(len=256) :: yline    ! line read
    CHARACTER(len=256) :: test     ! uppercase namelist group name
    INTEGER(i4)        :: stat     ! local copy of status variable
    INTEGER            :: ios      ! status variable from read operation
    LOGICAL            :: lrew     ! local copy of rewind flag
    INTEGER(i4)        :: iunit    ! local copy of unit number
    INTEGER(i4)        :: len_name ! length of requested namelist group name
    CHARACTER          :: ytest    ! character to test for delimiter
    CHARACTER(len=12)  :: code     ! error code printed
    INTEGER(i4)        :: ind      ! index from index routine
    INTEGER(i4)        :: indc     ! index of comment character (!)

    lrew  = .TRUE.
    IF (PRESENT(first)) lrew  = first
    iunit =  nunitnml
    IF (PRESENT(unit)) iunit = unit   
    stat  =  MISSING
    code  = 'MISSING'

    len_name = LEN_TRIM(name)

    IF (len_name > LEN(test)) THEN
       stat =  LENGTH_ERROR
       code = 'LENGTH_ERROR'
    END IF

    !test = '&'//tolower(name)
    write(test,'(A,A)') '&', tolower(name)

    ! Reposition file at beginning:
    IF (lrew) REWIND(iunit)

    ! Search start of namelist
    DO
       IF (stat /= MISSING) EXIT

       yline = ' '

       READ (iunit,'(a)',IOSTAT=ios) yline
       IF (ios < 0) THEN
          EXIT  ! MISSING
       ELSE IF (ios > 0) THEN
          stat =  READ_ERROR
          code = 'READ_ERROR'
          EXIT
       END IF

       yline = tolower(yline)

       ind = INDEX(yline,TRIM(test))

       IF (ind == 0) CYCLE

       indc = INDEX(yline,'!')

       IF (indc > 0 .AND. indc < ind) CYCLE

       ! test for delimiter
       ytest = yline(ind+len_name+1:ind+len_name+1)

       IF ( (LGE(ytest,'0') .AND. LLE(ytest,'9')) .OR. &
            (LGE(ytest,'a') .AND. LLE(ytest,'z')) .OR. &
            ytest == '_'                         .OR. &
            (LGE(ytest,'A') .AND. LLE(ytest,'Z'))) THEN
          CYCLE
       ELSE 
          stat = POSITIONED
          BACKSPACE(iunit)
          EXIT
       END IF
    END DO

    IF (PRESENT(status)) status = stat
    SELECT CASE (stat)
    CASE (POSITIONED)
       RETURN
    CASE (MISSING)
       IF (PRESENT(status)) RETURN
    END SELECT

    ! Error if it reaches here
    !message_text = 'namelist /'//TRIM(name)//'/ '//code
    write(message_text,'(A,A,A,A)') 'namelist /', trim(name), '/ ', trim(code)
    CALL finish('POSITION_NML',message_text)

  END SUBROUTINE position_nml

END MODULE mo_nml
