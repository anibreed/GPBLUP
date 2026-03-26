module renPED
  use M_ped
  use M_readfile
  use libsearch, only: Qsort 
  ! inbreed module is deprecated; Colleau-based routines used instead
  use colleau_mod
  use omp_lib
  implicit none
  type(ffh_t) :: h

  contains

  subroutine renumped(PEDFile, PEDinfo, REFFile, REFinfo, inb, rel, use_gpu, device_id)
  character(LEN=*),intent(in):: PEDFile
  integer,intent(in):: PEDinfo(:)
  character(LEN=*),intent(in),optional:: REFFile
  integer,intent(in),optional:: REFinfo(:)
  integer,intent(in),optional:: inb, rel
  logical, intent(in), optional :: use_gpu
  integer, intent(in), optional :: device_id
  type(CPED):: NX(3), OX(3)
  type(VPED),allocatable:: REF(:)
  integer,allocatable:: order(:),NPED(:,:)
  real(8),allocatable:: F(:)
  real(kind=r8), allocatable :: F_tab(:)
  character(LEN=LEN_STR),allocatable:: REFAnim(:)
  character(LEN=LEN_STR):: SV(MAX_VAR), string
  character(LEN=MAX_STR):: ofile, ofmt
  logical:: missBYR, inb_cal=.false., rel_cal=.false.
  integer:: IV(MAX_VAR), strlen, ianim
  integer(kind=ki4) :: NREC, NRECTot, unt, NA, n, i, j, ii, ix, jx, ios, stats(3), status
  integer(kind=ki4) :: k
  integer(kind=ki4) :: mapi, sire_idx, dam_idx
  logical :: found
  integer(kind=ki4), allocatable :: sire_arr(:), dam_arr(:)
  logical :: use_gpu_l
  integer :: di
  integer(kind=ki4), allocatable :: anc_s(:), anc_d(:)
  real(kind=r8), allocatable :: val_s(:), val_d(:)
  integer(kind=ki4) :: n_s, n_d, maxanc
  real(kind=r8) :: a_sd_re
  integer :: omp_nprocs, omp_threads
  integer :: cnt_gt1e6, cnt_gt1e8, cnt_gt1e12
  real(kind=r8) :: minpos, maxF, sumF
      integer :: outunit, l
      character(len=128) :: fname
      character(len=32) :: itmp, jtmp
      integer :: ns_nonzero, nd_nonzero, common_count
      real(kind=r8) :: sum_s, sum_d, min_s, max_s, min_d, max_d, sumprod
      type(CPED) :: tmpX
      integer :: idx_lookup

  if(present(REFFile)) then
   if(.not.present(REFinfo)) stop "please enter REFFile column"
   if(REFinfo(1) == missI) stop "REFFile column denoted by 0"
  endif
  if(present(inb)) then
  if(inb == 1) inb_cal=.true.
  endif
  if(present(rel)) then
  if(rel == 1) rel_cal=.true.
  endif  
  NREC=N_recf(PEDFile, MAX_RECL)
  unt = fopen(PEDFile,MAX_RECL)
  strlen=0
  NA=0

  do
   SV=missC
   IV=missI
   call initX(OX)
   call readrec(unt,MAX_RECL,n,SV,IV)
   if(n <= 0) exit
   NA=NA+1
   call checkPED(SV(PEDinfo(1:3)),NX)
   if(PEDinfo(4)/=0) NX(1)%BYR=IV(PEDinfo(4))
   if(PEDinfo(5)/=0) NX(1)%SEX=IV(PEDinfo(5))
   do i=1,3
    if( .not. Valid(NX(i)%ASD(1))) cycle
    call h%get_value(NX(i)%ASD(1), OX(i), stats(i))
   enddo
   call update(OX, NX, stats)
   do i=1,3
    if( .not. Valid(NX(i)%ASD(1))) cycle
    call h%ustore_value(NX(i)%ASD(1), NX(i))
   enddo           
  enddo
  close(unt)
  NREC=NA
  print*,"end read PEDFile data", NREC
  print*,"Ped File: Hash  size = ", h%hash_mask,"    Hash filled(No. of animal on PED) = ", h%n_keys_stored
  if(present(REFFile)) then
   if(rel_cal) then
    NREC= N_recf(REFFile)
    allocate(REFAnim(NREC))
    REFAnim=missC
   endif   
   unt = fopen(REFFile,MAX_RECL)
   NA=0
   do
     SV=missC
     IV=missI
     call initX(NX)
     call readrec(unt,MAX_RECL,n,SV,IV)
     if(n <= 0) exit
     NA=NA+1
     if(.not.Valid(SV(REFinfo(1)))) cycle
     call h%get_value(SV(REFinfo(1)), NX(1), status)
     if(status == -1) then
       print*,SV(REFinfo(1)),"not included on pedigree file"
       NX(1)%ASD(1)=SV(REFinfo(1))
       NX(1)%SEX=IV(REFinfo(2))
     endif
     NX(1)%GEN = 0
     NX(1)%OBS=NX(1)%OBS+1
     call h%ustore_value(NX(1)%ASD(1),NX(1))
     if(rel_cal) REFAnim(NA)=SV(REFinfo(1))
   enddo
   close(unt)
   NREC=NA
   print*,"end read REFFile data", NREC
 else
   do i = 0, h%n_buckets-1
     if (h%valid_index(i)) then
        h%vals(i)%GEN= 0
     endif
   enddo
 endif

 print*,"After reading Refence File: Hash  size = ", h%hash_mask,"    Hash filled(No. of animal on PED) = ", h%n_keys_stored
      
 call trace_PED()

 if(present(REFFile)) then
   do i = 0, h%n_buckets-1
      if (h%valid_index(i)) then
        if(h%vals(i)%GEN < 0) then
          call h%udelete_key(h%vals(i)%ASD(1))
        endif
      endif
   enddo
   print*,"After tracing PED: Hash  size = ", h%hash_mask,"    Hash filled(No. of animal on PED) = ", h%n_keys_stored
 !    call sleep(10)
 endif

 NREC=h%n_keys_stored
 
 allocate(REF(NREC),order(NREC))

 REF(:)%ID=missC
 REF(:)%NID=missI
 order(:)=missI
 NA=0
 missBYR=.false.
 do i = 0, h%n_buckets-1
    if(h%valid_index(i)) then
      NA=NA+1
      REF(NA)%ID=h%vals(i)%ASD(1)
      if(len_trim(REF(NA)%ID) > strlen) strlen=len_trim(REF(NA)%ID)
      REF(NA)%NID=h%vals(i)%BYR*100+(100 - h%vals(i)%GEN)
      order(NA)=i
    endif
 enddo
 if(NREC/=NA) stop "Valid_index error"
 NREC=NA
 call Qsort(REF(:)%NID,order)
 do i=1, NREC
  h%vals(order(i))%NASD(1)=i
 enddo
  ! After sorting `REF(:)%NID` with `order`, the REF(:)%ID array
  ! is no longer aligned with `order`. Rebuild REF(:)%ID from
  ! the hash slots referenced by `order` so IDs and order match.
  do i = 1, NREC
    REF(i)%ID = h%vals(order(i))%ASD(1)
  end do
 do i=1,NREC
    ix = h%get_index(h%vals(order(i))%ASD(2))
    if(ix /= -1) then
      ! Always map sire index if found (avoid leaving parent index as missing)
      h%vals(order(i))%NASD(2) = h%vals(ix)%NASD(1)
    endif

  ! Configure OpenMP threads: use single thread for deterministic diagnostics
  omp_nprocs = omp_get_num_procs()
  omp_threads = 1
  call omp_set_num_threads(omp_threads)
  write(*,'(A,I4,2X,A,I4)') 'OpenMP: using threads=', omp_threads, ' (procs=', omp_nprocs,')'
    ix = h%get_index(h%vals(order(i))%ASD(3))
    if(ix /= -1) then
      ! Always map dam index if found
      h%vals(order(i))%NASD(3) = h%vals(ix)%NASD(1)
    endif
 enddo

 if(inb_cal) then
  ! Use Colleau inbreeding routine for efficiency and optional triplet output
  allocate(F(NREC))
  allocate(sire_arr(NREC))
  allocate(dam_arr(NREC))
  use_gpu_l = .false.
  if (present(use_gpu)) use_gpu_l = use_gpu
  do i = 1, NREC
    if (h%vals(order(i))%NASD(2) == missI) then
      sire_arr(i) = 0
    else
      sire_arr(i) = h%vals(order(i))%NASD(2)
    end if
    if (h%vals(order(i))%NASD(3) == missI) then
      dam_arr(i) = 0
    else
      dam_arr(i) = h%vals(order(i))%NASD(3)
    end if
  end do

    ! Summary: prepared parent index arrays (verbose debug suppressed)
    write(*,*) 'Prepared sire/dam arrays for', NREC, 'records'

  ! Call the Colleau inbreeding implementation (compute_inbreeding)
  write(*,*) '>>> renumped: calling compute_inbreeding (colleau_mod)'
  call timestamp()
  call set_colleau_verbose(.true.)
  ! Per-target diagnostics removed; proceed to compute inbreeding
  maxanc = min(NREC, MAXANC_LIMIT)
  ! Diagnostic: sample mapping of indices -> IDs and founder counts
  do i = 1, min(50, NREC)
    write(*,*) 'MAP index=', i, 'ID=', trim(REF(i)%ID)
    if (sire_arr(i) /= 0) then
      write(*,*) '  SIRE_IDX=', sire_arr(i), 'SIRE_ID=', trim(REF(sire_arr(i))%ID)
    end if
    if (dam_arr(i) /= 0) then
      write(*,*) '  DAM_IDX=', dam_arr(i), 'DAM_ID=', trim(REF(dam_arr(i))%ID)
    end if
  end do
  ! founder counts
  ii = 0
  do i = 1, NREC
    if (sire_arr(i) == 0 .or. dam_arr(i) == 0) ii = ii + 1
  end do
  write(*,'(A,I8)') 'Number of records with at least one missing parent:', ii
  ! (previous per-target presence checks removed)
  ! Compute inbreeding using Colleau sparse routine (production)
  call compute_inbreeding(NREC, sire_arr, dam_arr, F)

  ! Tabular fallback and comparison removed for production build
  ! quick check: print first few F values returned
  do i = 1, min(10,NREC)
    write(*,'(A,I6,2X,ES18.12)') 'SAMPLE F[', i, F(i)
  end do
  ! Print count of non-zero F entries
  ii = 0
  do i = 1, NREC
    if (abs(F(i)) > 1.0d-12) ii = ii + 1
  end do
  write(*,'(A,I8)') 'Number of non-zero F values:', ii
  ! Additional diagnostics: F distribution and min/max
  cnt_gt1e6 = 0; cnt_gt1e8 = 0; cnt_gt1e12 = 0
  minpos = 1.0e30_r8; maxF = 0.0_r8; sumF = 0.0_r8
  do i = 1, NREC
    if (F(i) > 1.0d-6) cnt_gt1e6 = cnt_gt1e6 + 1
    if (F(i) > 1.0d-8) cnt_gt1e8 = cnt_gt1e8 + 1
    if (F(i) > 1.0d-12) cnt_gt1e12 = cnt_gt1e12 + 1
    if (F(i) > 0.0_r8) then
      if (F(i) < minpos) minpos = F(i)
      if (F(i) > maxF) maxF = F(i)
      sumF = sumF + F(i)
    end if
  end do
  if (minpos > 1.0e29_r8) minpos = 0.0_r8
  write(*,'(A,I8)') 'Count F>1e-6:', cnt_gt1e6
  write(*,'(A,I8)') 'Count F>1e-8:', cnt_gt1e8
  write(*,'(A,I8)') 'Count F>1e-12:', cnt_gt1e12
  write(*,'(A,2X,A,ES16.8)') 'min positive F=', minpos
  write(*,'(A,2X,A,ES16.8)') 'max F=', maxF
  write(*,'(A,2X,A,ES16.8)') 'sum F=', sumF
  if (rel_cal) then
    write(*,*) '>>> renumped: writing triplets A_triplets.txt'
    call write_triplets(NREC, sire_arr, dam_arr, F, 'A_triplets.txt')
    write(*,*) '<<< renumped: triplets written'
  endif
  write(*,*) '<<< renumped: returned from compute_inbreeding'
  call timestamp()

  ! Computed inbreeding F stored for each record (detailed output suppressed)

  ! store F (double precision) without scaling
  do i = 1, NREC
    ! Prefer current hash lookup index to ensure we write to the live slot
    idx_lookup = h%get_index(trim(REF(i)%ID))
    if (idx_lookup /= -1) then
      h%vals(idx_lookup)%F = F(i)
    else
      h%vals(order(i))%F = F(i)
    end if
  end do
  
  deallocate(F, sire_arr, dam_arr)

endif

 if(present(REFFile)) then
   ofile=trim(REFFile)//"_ren"
 else
   ofile=trim(PEDFile)//"_ren"
 endif
 unt=fopen(ofile,MAX_RECL)
  call writenum(strlen,string,"i10")
  ofmt="),i2,i9,i3,i5,f10.6,3(1x,a"//trim(string)//"))"
 strlen=floor(log10(real(NREC)) + 1)
 call writenum(strlen,string,"i10")
 ofmt="(3(1x,i"//trim(string)//ofmt

 do i=1,NREC
   write(unt,fmt=ofmt) h%vals(order(i))%NASD(:), &
           h%vals(order(i))%SEX, &
           h%vals(order(i))%BYR, &
           h%vals(order(i))%GEN, &
           h%vals(order(i))%OBS, &
           h%vals(order(i))%F, &
           h%vals(order(i))%ASD(:)
 enddo
 close(unit=unt)
 end subroutine

 subroutine checkPED(SV,X)
 character(len=*),intent(in):: SV(:)
 type(CPED),intent(out):: X(:)
 character(len=len(SV(1))):: TV(size(SV))
 integer:: i
 call initX(X)
  if(.not.Valid(SV(1))) return
  TV=SV
  ! check same ID between Sire and dam
  if(SV(2)==SV(3)) then
     TV(2)=missA
     TV(3)=missA
  endif
  ! check same ID between animal and parents
  do i = 2, 3
       if(SV(1) == SV(i))  TV(i)=missA
  enddo
  do i =1, 3
      X(1)%ASD(i)=TV(i)
  enddo
  X(2)%ASD(1)=TV(2)
  X(2)%SEX=Male
  X(3)%ASD(1)=TV(3)
  X(3)%SEX=Female
 end subroutine

  subroutine update(A,B, stat)  ! modify A to B
  type(CPED),intent(inout):: A(:), B(:)
  integer,intent(in):: stat(:)
  integer:: ii, jj, status
  do  ii = 1, 3
     if(stat(ii) == -1) cycle 
    do jj=2, 3
       if(A(ii)%ASD(jj) /= missA .and. B(ii)%ASD(jj) == missA) then
         B(ii)%ASD(jj) = A(ii)%ASD(jj)
       elseif(A(ii)%ASD(jj) /= missA .and. B(ii)%ASD(jj) /= missA) then
         if(A(ii)%ASD(jj) /= B(ii)%ASD(jj)) B(ii)%ASD(jj)= missA
       endif
    enddo
    if(A(ii)%BYR /= missI .and. B(ii)%BYR == missI) B(ii)%BYR=A(ii)%BYR
    if(A(ii)%SEX /= missI) then
      if(B(ii)%SEX == missI) then
        B(ii)%SEX=A(ii)%SEX
      elseif(B(ii)%SEX /= missI) then
        if(A(ii)%SEX /= B(ii)%SEX ) then
          call h%delete_key(A(ii)%ASD(1),status)
          B(ii)%ASD(1)=missA
          if(ii > 1) B(1)%ASD(ii)=missA
       endif
      endif
    endif
  enddo
  return
  end subroutine

  subroutine trace_PED()
  type(CPED):: valA,valS,valD
  integer:: NC, newNC, round, i, j,status
  integer,allocatable:: check(:)
  round=0
  NC= huge(ki4)
  allocate(check(h%n_keys_stored))
  do 
  round=round+1
  newNC=0
  check=0

  do i = 0, h%n_buckets-1
    status = h%get_index(h%vals(i)%ASD(1))
    if (status == -1) cycle
    valA = h%vals(status)
    if(valA%ASD(2) /= missA) then
      call initX(valS)
      call h%get_value(valA%ASD(2), valS,status)
      if(status /= -1) then
        if(valS%GEN <= valA%GEN) then
           valS%GEN = valA%GEN+1
           call h%ustore_value(valS%ASD(1), valS)
           newNC=newNC+1
        endif   
      endif
    endif
    if(valA%ASD(3) /= missA) then
      call initX(valD)
      call h%get_value(valA%ASD(3), valD,status)
      if(status /= -1) then
        if(valD%GEN <= valA%GEN) then
           valD%GEN = valA%GEN+1
           call h%ustore_value(valD%ASD(1), valD)
           newNC=newNC+1
        endif   
      endif
    endif
   enddo

   if(newNC <= missI)  then
      print*,round," th trace back finished, all Animals traced!!"
      exit
   elseif(newNC < NC) then
      NC=newNC
      print*,round," th trace back",NC,"  Animals updated!!"

   elseif(newNC == NC) then
      do i =1, newNC
        if(check_IDD(check(i)))  exit
      enddo
      do i = 0, h%n_buckets-1
         if (h%valid_index(i)) then
          if(h%vals(i)%GEN > 0 ) h%vals(i)%GEN = 0 
         endif
      enddo
      NC=huge(ki4)
   endif
  enddo
 end subroutine


logical function check_IDD(step) result(check)
 integer,intent(in):: step
 integer:: first = 1
 check=.false.
 if(check_ID(h%vals(step)%ASD(1), first)) then
  check=.true.
 endif
end function

recursive logical function check_ID(ARN,firstA) result(tangled)
 character(len=*),intent(in):: ARN
 integer,intent(in),optional::firstA
 character(len=LEN_STR),save:: ARN0
 integer,save:: ngen
 integer:: stepA
 tangled=.false.
 if (ARN == missA) return
 if(present(firstA)) then
   ARN0=ARN
   ngen=0
 elseif(ARN == ARN0) then
   stepA = h%get_index(ARN)
   if(stepA == -1) return
   print*," Animal Tangled!! delete Sire and Dam    ",h%vals(stepA)%ASD(:)
   h%vals(stepA)%ASD(2:3) = missA    !    call delete_stack(step)
   tangled=.true.
   return
 endif
 ngen=ngen+1
 stepA = h%get_index(ARN)
 if(stepA == -1) return
 if(h%vals(stepA)%ASD(2) /= missA) then
     tangled=check_ID(h%vals(stepA)%ASD(2))
     if(tangled) return
 endif
 if(h%vals(stepA)%ASD(3) /= missA) then
     tangled=check_ID(h%vals(stepA)%ASD(3))
     if(tangled) return
 endif
end function check_ID

elemental subroutine initX(X)
 TYPE(CPED), intent(inout) :: X
 integer:: ii
 do ii=1,3
    X%ASD(ii)=missA
    X%NASD(ii)=missI
 enddo
 X%SEX=missI
 X%GEN= -1
 X%BYR=missI
 X%OBS=missI
 X%F=missR
end subroutine

   logical function Valid(XC)
   character(len=*),intent(in):: XC
   character(len=len(XC)):: XX
   integer:: N
   Valid = .true.
   XX=trim(adjustl(XC))
   if(iachar(XX(1:1)) < 33 .or. iachar(XX(1:1)) > 126 ) then
    Valid = .false.
   else
    N = len_trim(XX)
    select case (N) 
     case( :0 )
       Valid = .false.
     case( 1 )
       if(XX(1:1) == '0') Valid = .false.
     case( 2: )
       Valid = .true.
    end select
   endif
   end function

end module
