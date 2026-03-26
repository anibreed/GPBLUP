program relped
 use renPED
 implicit none
 character(len=MAX_STR):: PARAM_File, PEDFile, REFFile
 integer:: PEDinfo(5), REFinfo(2), inb, rel
 logical:: reference
 logical :: use_gpu
 integer :: device_id
 reference=.false. 
 PARAM_File=missC
 PEDFile=missC
 REFFile=missC
 PEDinfo=0
 REFinfo=0
 inb=0
 rel=0
 use_gpu = .false.
 device_id = 0

 call getarg(1,PARAM_File)
 PARAM_File=trim(adjustl(PARAM_File))
 if((len_trim(PARAM_File)).lt.1) stop 'Input parameter file'
 print*,"Finish to read parameter file=  ",trim(PARAM_File)

 call read_parameter(PARAM_File)

 if(reference) then
     call renumped(PEDFile=PEDFile, PEDinfo=PEDinfo, REFFile=REFFile, REFinfo=REFinfo, inb=inb, rel=rel, use_gpu=use_gpu, device_id=device_id)
 else
     call renumped(PEDFile=PEDFile, PEDinfo=PEDinfo, inb=inb, rel=rel, use_gpu=use_gpu, device_id=device_id)
 endif

 contains
 subroutine read_parameter(paramfile)
 character(len=*),intent(in):: paramfile
 character(len=MAX_STR):: XC(MAX_VAR)
 integer:: X(MAX_VAR), n, unt

 unt=fopen(paramfile,MAX_RECL)
 do
        call readline(unt,n,XC,X)
        if(n<= -1) exit
        XC(1)=upper(XC(1))
        select case(XC(1)(1:8))
                 case("PEDFILE:")
                     PEDFile=XC(2)
                     print*,"Pedigree file=  ", trim(PEDFile)
                 case("PEDINFO:")
                     if(n<3) then
                         print*,"Pedigree information was not enough"
                         stop
                     else
                         PEDinfo(1:n-1)=X(2:n)
                         print*,"Column info on pedigree file(ARN, SRN, DRN, BYR, SEX)=",PEDinfo
                     endif
                 case("REFFILE:")
                     REFFile=XC(2)
                     reference=.true.
                     print*,"Reference file name for renumbering=  ",trim(REFFile)
                 case("REFINFO:")
                     if(n<3) then
                         print*,"Reference file information was not enough"
                         stop
                     else
                         REFinfo(1)=X(2)  ! Column no. for Animal ID
                         REFinfo(2)=X(3)  ! Column no. for Animal Sex(1=Female, 2=Make)
                         print*,"Column infos for animal ID and Sex on reference file=",REFinfo

                     endif
                 case("INBREED:")
                     if(n<2) then
                         print*,"Information for calculating inb. was not enough"
                         stop
                     else
                         inb=X(2)
                         if(inb<0.or.inb>1) then
                              print*,"Information for calculating inb. should be 0 or 1",inb
                              stop
                         else
                              print*,"Calculating inb. (No=0, Yes=1)",inb
                         endif    
                     endif
                 case("RELATION")
                     if(n<2) then
                         print*,"Information for calculating A matrix was not enough"
                         stop
                     else
                         rel=X(2)
                         if(rel<0.or.rel>1) then
                              print*,"Information for calculating A matrix should be 0 or 1",rel
                              stop
                         else
                              print*,"Calculating A matrix (No=0, Yes=1)",rel
                         endif    
                     endif
                 case default 
                     continue
        end select
 enddo
 close(unit=unt)
 end subroutine
 end program
