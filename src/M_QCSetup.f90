module M_QCSetup
  use M_Kinds
  use M_Variables
  use M_ReadPar
  use M_Hashtable
  use H_PedRead
  use M_StrEdit
  implicit none

  private
  public :: setup_parameters_mod
  public :: resolve_animal_key_fields_mod
  public :: load_pedigree_info_mod
  public :: set_qcsetup_debug_mod

  logical, save :: qcsetup_debug = .false.

contains

  subroutine set_qcsetup_debug_mod(in_debug)
    logical, intent(in) :: in_debug
    qcsetup_debug = in_debug
  end subroutine set_qcsetup_debug_mod

  subroutine resolve_animal_key_fields_mod(out_ped_key_field, out_geno_key_field)
    character(len=LEN_STR), intent(out) :: out_ped_key_field
    character(len=LEN_STR), intent(out) :: out_geno_key_field

    integer :: j

    out_ped_key_field = ''
    out_geno_key_field = ''

    do j = 1, PEDFile%NVAR
      if (index(trim(to_upper(PEDFile%FieldName(j))), 'ANIMAL_') == 1) then
        out_ped_key_field = trim(PEDFile%FieldName(j))
        exit
      end if
    end do

    do j = 1, GENOFile%NVAR
      if (index(trim(to_upper(GENOFile%FieldName(j))), 'ANIMAL_') == 1) then
        out_geno_key_field = trim(GENOFile%FieldName(j))
        exit
      end if
    end do

    if (len_trim(out_ped_key_field) == 0) then
      print *, 'ERROR: PEDFILE does not define an ANIMAL_ key field.'
      print *, '       Please define a field name starting with ANIMAL_ (e.g. ANIMAL_ID).'
      stop
    end if

    if (len_trim(out_geno_key_field) == 0) then
      print *, 'ERROR: GENOFILE does not define an ANIMAL_ key field.'
      print *, '       Please define a field name starting with ANIMAL_ (e.g. ANIMAL_ID).'
      stop
    end if

    print *, '  Matching key fields resolved from parameter file:'
    print *, '    PED key  : ', trim(out_ped_key_field)
    print *, '    GENO key : ', trim(out_geno_key_field)
  end subroutine resolve_animal_key_fields_mod

  subroutine load_pedigree_info_mod(io_pht, out_ped_file, out_n_animals, in_ped_key_field)
    type(PEDHashTable), intent(inout) :: io_pht
    character(len=256), intent(out) :: out_ped_file
    integer, intent(out) :: out_n_animals
    character(len=*), intent(in) :: in_ped_key_field

    integer :: j, unitF, local_n_animals
    logical :: found

    if (len_trim(PEDFile%FileName) == 0) then
      print *, 'ERROR: PEDFile not initialized (read_parameters not called)'
      out_n_animals = 0
      return
    end if

    out_ped_file = trim(PEDFile%FileName)

    print *, '  Loading PED into hashtable: ', trim(out_ped_file)
    print *, '  PED key field used for matching: ', trim(in_ped_key_field)
    call pht_load_ped_file(io_pht, local_n_animals, unitF, trim(in_ped_key_field))
    out_n_animals = local_n_animals

    if (qcsetup_debug) then
      print *, '  --- PED HashTable Load Test ---'
      print *, '  Total records in Hashtable: ', io_pht%count

      found = .false.
      do j = 0, io_pht%size - 1
        if (associated(io_pht%table(j)%next)) then
          print *, '  Sample record from HashTable bucket ', j
          print *, '    ID:   ', trim(io_pht%table(j)%next%ped_data%ID)
          print *, '    ARN:  ', trim(io_pht%table(j)%next%ped_data%ARN)
          print *, '    Sex:  ', io_pht%table(j)%next%ped_data%Sex
          print *, '    Sire: ', trim(io_pht%table(j)%next%ped_data%Sire)
          print *, '    Dam:  ', trim(io_pht%table(j)%next%ped_data%Dam)
          found = .true.
          exit
        end if
      end do
      if (.not. found) print *, '  WARNING: HashTable seems empty after load.'
      print *, '  ---------------------------------------'
    end if

    close(unitF)
    print *, '  Finished loading pedigree info into HashTable.'
  end subroutine load_pedigree_info_mod

  subroutine setup_parameters_mod(param_file)
    character(len=*), intent(in) :: param_file
    logical :: exists

    print *, 'Loading parameters from: ', trim(param_file)

    inquire(file=trim(param_file), exist=exists)
    if (.not. exists) then
      print *, 'ERROR: Parameter file not found: ', trim(param_file)
      stop
    end if

    call read_parameters(trim(param_file))

    if (len_trim(PEDFile%FileName) == 0) then
      print *, 'ERROR: PED_FILE not specified in parameter file'
      stop
    end if

    if (len_trim(MAPFile%FileName) == 0) then
      print *, 'ERROR: MAP_FILE not specified in parameter file'
      stop
    end if

    print *, 'Parameter file loaded successfully'
    if (qcsetup_debug) then
      print *, ''
      print *, 'File Information from Parameter File:'
      print *, '  PED file: ', trim(PEDFile%FileName)
      if (PEDFile%NVAR > 0) then
        print *, '    Fields: ', PEDFile%NVAR
        print *, '    Header lines: ', PEDFile%Header
      end if

      print *, '  MAP file: ', trim(MAPFile%FileName)
      if (MAPFile%NVAR > 0) then
        print *, '    Fields: ', MAPFile%NVAR
        print *, '    Header lines: ', MAPFile%Header
      end if

      if (len_trim(GENOFile%FileName) > 0) then
        print *, '  GENO file: ', trim(GENOFile%FileName)
        if (GENOFile%NVAR > 0) then
          print *, '    Fields: ', GENOFile%NVAR
          print *, '    Header lines: ', GENOFile%Header
        end if
      end if
    end if
  end subroutine setup_parameters_mod

end module M_QCSetup
