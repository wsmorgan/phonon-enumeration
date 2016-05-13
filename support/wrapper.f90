!!<summary>Subroutines needed to provide easy access to random selections of
!!enumerate structures under the enum4 standard.</summary>
module wrapper
  use io_utils, only: write_rotperms_list
  use enumeration_types, only: RotPermList

  public write_matrices
  private
contains
  !!<summary>Returns lowest i/o unit number not in use.</summary>
  !!<parameter name="unit">Out parameter that will contain the lowest i/o number.</parameter>
  integer function fpy_newunit(unit)
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
    integer :: n

    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          fpy_newunit = n
          return
       end if
    end do
    stop "newunit ERROR: available unit not found."
  end function fpy_newunit

  !!<summary>Writes out the HNF, SNF and Left Transform matrices to file.</summary>
  !!<parameter name="HNFs,SNFs,Ls">List of all unique HNFs and SNFs with corresponding left transform matrices.</parameter>
  subroutine write_matrices(HNFs, SNFs, Ls, iperm, rotperms, ivol)
    integer, intent(in), dimension(:,:,:) :: HNFs, SNFs, Ls
    integer, intent(in) :: iperm(:), ivol
    type(RotPermList), intent(in) :: rotperms(:)
    !!<local name="nHNF">Number of HNFs in the list that match the current block index, HNFi.</local>
    !!<local name="vsH">Vector subscript for matching the HNF index to the HNF list.</local>
    !!<local name="funit">I/O unit for the file to write the output to.</local>
    integer :: nHNF, HNFi
    integer, allocatable :: vsH(:)
    logical :: exist
    integer :: ioerr, iHNF, jHNF
    integer :: funit
    character(len=:), allocatable :: setname
    character(80) :: groupname
    character(10) :: cellstr, groupstr
    
    !Write the size of the cell to string so we can concat it to the folder name.
    write(cellstr, "(I6)") ivol
    call system("mkdir -p cells."//trim(adjustl(cellstr)))
    setname = "cells."//trim(adjustl(cellstr))//"/matrices"
    
    inquire(file=setname, exist=exist)
    if (exist) then
       open(fpy_newunit(funit), file=setname, status="old", position="append",  &
            action="write", iostat=ioerr)
    else
       open(fpy_newunit(funit), file=setname, status="new", action="write", iostat=ioerr)
    end if
    
    do HNFi=1, maxval(iperm)
       !Write the symmetry group for this set of HNFs, SNFs and Ls. We adjust the file
       !name to have the HNFi index so that we don't have to duplicate the group for
       !multiple HNFs having the same group
       write(groupstr, "(I6)") HNFi
       groupname = "cells."//trim(adjustl(cellstr))//"/group."//trim(adjustl(groupstr))
       call write_rotperms_list(rotperms(HNFi), groupname)

       !Next, we loop over those HNFs and SNFs that match the permutation index we are
       !looking at currently
       nHNF = count(iperm==HNFi)
       allocate(vsH(nHNF))
       vsH = pack((/(i,i=1,size(HNFs, 3))/), HNFi==iperm)
              
       do iHNF = 1, nHNF ! Write this labeling for each corresponding HNF
          jHNF = vsH(iHNF) ! Index of matching HNFs
          write (funit, '((i4,1x),3(i2,1x),2x,6(i2,1x),2x,9(i4,1x))') &
               HNFi, SNFs(1,1,jHNF),SNFs(2,2,jHNF),SNFs(3,3,jHNF), &
               HNFs(1,1,jHNF),HNFs(2,1,jHNF),HNFs(2,2,jHNF),HNFs(3,1,jHNF),HNFs(3,2,jHNF),HNFs(3,3,jHNF), &
               transpose(Ls(:,:,jHNF))
       end do
       deallocate(vsH)
    end do
    close(funit)
  end subroutine write_matrices
end module wrapper
