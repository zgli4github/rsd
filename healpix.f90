module healpix_tool
  
  implicit none
  contains
  
  subroutine Healpixmap2fits_simple(nside,fname,map,ttype,tunit,extname,pixtype,coordsys,ordering)
    implicit none
    character(len=*) :: fname,ttype,tunit,extname,pixtype,coordsys,ordering
    real :: map(:)
    integer :: nside,npix,i,j,k,ipnest,firstpix,lastpix
    integer :: status,rwmode,blocksize,unit,pcount=0,gcount=1,bitpix=32
    integer :: tfields,naxis,nrow,ncol,varidat=0
    character*68 :: tform="E"
    character*80 :: card

    status = 0
    rwmode = 1
    blocksize = 1
    tfields = 1 ! number of columns.
    naxis = 1   ! dimension of data
    !ttype = "temprature"
    !tunit = "uK"
    !extname = "COMP_map"
    !pixtype = 'HEALPIX'
    !coordsys = 'GALACTIC'
    !ordering = 'NESTED'

    npix=12*nside**2
    k=size(map)
    if(k.ne.npix) then 
      write(*,*) "Error! Size of input map is not equal to number of Healpix pixels."
      write(*,*) k,npix,nside
      !stop  ! uncomment it just for arbitrary input array
    end if
    nrow = k
    firstpix=0
    lastpix=nrow-1
    unit=12
    !==== delete fname if it exists. Attention: It can only delete fits file !
    !==== It will do nothing if fname is not a fits file, 
    call deletefile(fname) 
    !==== create an empty fits file.
    call ftinit(unit,fname,blocksize,status)
    if(status.ne.0) call printerr(status)
    write(*,*) "hello2"
    !==== create a null primary hdu,necessary !!!
    call ftiimg(unit,bitpix,0,0,status)
    if(status.ne.0) call printerr(status)
    !==== create a binary hdu and write keywods in it.
    write(*,*) "hello3"
    call ftibin(unit,nrow,tfields,ttype,tform,tunit,extname,varidat,status)
    if(status.ne.0) call printerr(status)
    write(*,*) "create a binary hdu and write keywods in it."
    !==== write data in table.
    write(*,*) "writing data in binary table."
    call ftpcle(unit,1,1,1,nrow,map,status)
    if(status.ne.0) call printerr(status)
    call ftgcrd(unit,"extname",card,status) 
      ! Get the 80-character header record for the named keyword, "extname".
      ! This is necessary because subroutine 'ftiky[jks]' will insert a keyword 
      !      immediately following the last keyword that has been read from chu.
      ! Without this call, ftikyj will insert a keyword in the first line,which 
      ! is supposed to be the 'xtension' keyword.
    write(*,*) "writing keyword in header."
    call ftikys(unit,"PIXTYPE",pixtype," ",status)
    call ftikys(unit,"COORDSYS",coordsys,"Coordinate system",status)
    call ftikys(unit,"ORDERING",ordering,"Healpix Ordering",status)
    call ftikyj(unit,"NSIDE",nside,"Healpix Nside",status)
    call ftikyj(unit,"FIRSTPIX",firstpix," ",status)
    call ftikyj(unit,"LASTPIX",lastpix," ",status)
    !==== close file.
    write(*,*) "close file."
    call ftclos(unit,status)
    if(status.ne.0) call printerr(status)
    return
  end subroutine

  subroutine Healpixmap2fits_simpleD(nside,fname,map,ttype,tunit,extname,pixtype,coordsys,ordering)
    implicit none
    character(len=*) :: fname,ttype,tunit,extname,pixtype,coordsys,ordering
    real*8 :: map(:)
    integer :: nside,npix,i,j,k,ipnest,firstpix,lastpix
    integer :: status,rwmode,blocksize,unit,pcount=0,gcount=1,bitpix=32
    integer :: tfields,naxis,nrow,ncol,varidat=0
    character*68 :: tform="D"
    character*80 :: card

    status = 0
    rwmode = 1
    blocksize = 1
    tfields = 1 ! number of columns.
    naxis = 1   ! dimension of data
    npix=12*nside**2
    k=size(map)
    if(k.ne.npix) then 
      write(*,*) "Error! Size of input map is not equal to number of Healpix pixels."
      write(*,*) k,npix,nside
      stop
    end if
    nrow = npix
    firstpix=0
    lastpix=nrow-1
    unit=12
    call deletefile(fname) 
    call ftinit(unit,fname,blocksize,status)
    call ftiimg(unit,bitpix,0,0,status)
    write(*,*) "create a binary hdu and write keywods in it."
    call ftibin(unit,nrow,tfields,ttype,tform,tunit,extname,varidat,status)
    if(status.ne.0) call printerr(status)
    write(*,*) "writing data in binary table."
    call ftpcld(unit,1,1,1,nrow,map(0:npix-1),status)
    if(status.ne.0) call printerr(status)
    call ftgcrd(unit,"extname",card,status) 
    write(*,*) "writing keyword in header."
    call ftikys(unit,"PIXTYPE",pixtype," ",status)
    call ftikys(unit,"COORDSYS",coordsys,"Coordinate system",status)
    call ftikys(unit,"ORDERING",ordering,"Healpix Ordering",status)
    call ftikyj(unit,"NSIDE",nside,"Healpix Nside",status)
    call ftikyj(unit,"FIRSTPIX",firstpix," ",status)
    call ftikyj(unit,"LASTPIX",lastpix," ",status)
    write(*,*) "close file."
    call ftclos(unit,status)
    if(status.ne.0) call printerr(status)
    return
  end subroutine

  subroutine Healpixmap2fits_simpleI(nside,fname,map,ttype,tunit,extname,pixtype,coordsys,ordering)
    implicit none
    character(len=*) :: fname,ttype,tunit,extname,pixtype,coordsys,ordering
    integer :: map(:)
    integer :: nside,npix,i,j,k,ipnest,firstpix,lastpix
    integer :: status,rwmode,blocksize,unit,pcount=0,gcount=1,bitpix=32
    integer :: tfields,naxis,nrow,ncol,varidat=0
    character*68 :: tform="J"
    character*80 :: card

    status = 0
    rwmode = 1
    blocksize = 1
    tfields = 1 ! number of columns.
    naxis = 1   ! dimension of data
    npix=12*nside**2
    k=size(map)
    if(k.ne.npix) then 
      write(*,*) "Error! Size of input map is not equal to number of Healpix pixels."
      write(*,*) k,npix,nside
      stop
    end if
    nrow = npix
    firstpix=0
    lastpix=nrow-1
    unit=12
    call deletefile(fname) 
    call ftinit(unit,fname,blocksize,status)
    call ftiimg(unit,bitpix,0,0,status)
    write(*,*) "create a binary hdu and write keywods in it."
    call ftibin(unit,nrow,tfields,ttype,tform,tunit,extname,varidat,status)
    if(status.ne.0) call printerr(status)
    write(*,*) "writing data in binary table."
    call ftpclj(unit,1,1,1,nrow,map(0:npix-1),status)
    if(status.ne.0) call printerr(status)
    call ftgcrd(unit,"extname",card,status) 
    write(*,*) "writing keyword in header."
    call ftikys(unit,"PIXTYPE",pixtype," ",status)
    call ftikys(unit,"COORDSYS",coordsys,"Coordinate system",status)
    call ftikys(unit,"ORDERING",ordering,"Healpix Ordering",status)
    call ftikyj(unit,"NSIDE",nside,"Healpix Nside",status)
    call ftikyj(unit,"FIRSTPIX",firstpix," ",status)
    call ftikyj(unit,"LASTPIX",lastpix," ",status)
    write(*,*) "close file."
    call ftclos(unit,status)
    if(status.ne.0) call printerr(status)
    return
  end subroutine



  
  subroutine printerr(sta)
  !======================================================
  !  Printerr subroutine widely used in above code.
  !======================================================
    implicit none 
    integer,intent(in) :: sta
    character*30 :: errtxt
    character*80 :: errmsg
    if(sta.eq.0) return
    call ftgerr(sta,errtxt)
    write(*,*) "FITSIO error status =",sta,":",errtxt
    call ftgmsg(errmsg)
    do while(errmsg.ne.'')
      write(*,*) errmsg
      call ftgmsg(errmsg)
    end do
    call ftcmsg
    if(sta.ne.107) stop 
  end subroutine

  subroutine deletefile(filename)
    implicit none
    character(len=*) :: filename
    character(len=200) :: fname
    logical :: alive
    integer :: unit
    unit=99
    fname=trim(filename)
    inquire(file=fname,exist=alive)
    if(alive) then
      write(*,*) "Input file exists. Now delete it."
      write(*,*) trim(fname)
      open(12,file=fname)
      close(12,status="delete")
      inquire(file=fname,exist=alive)
      if(alive) then
        write(*,*) "Error: Delete file failed."
        stop
      else
        write(*,*) "File deleted."
      end if
    else
      write(*,*) "Input file does not exist."
      write(*,*) trim(fname)
    end if
  end subroutine

!  subroutine deletefile(filename,status)
!  !======================================================
!  ! Delete FITS file if it exists. 
!  !======================================================
!    implicit none
!    integer status,unit,blocksize
!    character*(*) :: filename
!
!    if(status.gt.0) return
!    call ftgiou(unit,status)
!    call ftopen(unit,filename,1,blocksize,status)
!    if(status.eq.0) then
!      call ftdelt(unit,status)
!    else if(status.eq.103) then
!      status = 0
!      call ftcmsg
!    else
!      status = 0
!      call ftcmsg
!      call ftdelt(unit,status)
!    end if
!    call ftfiou(unit,status)
!  end subroutine

end module
