module horizonrun
  use pix_tools,only: pix2ang_nest,ang2pix_nest,pix2vec_nest,vec2pix_nest,ang2vec,vec2ang
  use fitstools,only: output_map,input_map
  use head_fits,only: add_card
  use healpix_tool

  implicit none
  type :: hrmock
    integer :: ngalaxy
    real*8 :: rmin,rmax
    real*8,allocatable :: mass(:),x(:),y(:),z(:),vx(:),vy(:),vz(:)
    real*8,allocatable :: map2048_healpix(:)
  end type
  type(hrmock) :: hr3mock(0:26),hr2mock(0:7)
  character(len=200) :: HR3Home="/home/zg/work/galdata/simulation/HorizonRun/hr3/"
  type :: nzparam
    integer :: nbin,islog
    real*8 :: rmin,rmax
    real*8,allocatable :: rmean(:),vol(:),Nz(:)
  end type
  type(nzparam) :: hr3nz(0:26)
  real*8,parameter :: Om_HR3=0.26d0,OmL_HR3=0.74d0,h_HR3=0.72d0,Omb_HR3=0.044d0

  contains

  subroutine hr3_readmass(imock)
    implicit none
    integer :: imock,i,j,k
    real,allocatable :: mass(:)
    character(len=200) :: massf,cmock
    if(imock.eq.0) then
      cmock="00"
    else if(imock.lt.10) then
      write(cmock,"(i1)") imock
      cmock="0"//trim(adjustl(cmock))
    else
      write(cmock,"(i2)") imock
      cmock=trim(adjustl(cmock))
    end if
    massf=trim(HR3Home)//"mock3"//trim(cmock)//"/HR3_mass_sdss3_"//trim(cmock)//".dat"
    open(10,file=massf,form="binary",status="old")
    read(10) k
    hr3mock(imock)%ngalaxy=k
    allocate(mass(k),hr3mock(imock)%mass(k))
    read(10) mass
    close(10)
    hr3mock(imock)%mass=mass
    return
  end subroutine
end module
