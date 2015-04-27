use const
use distances
use coordtransform
use horizonrun
use healpix_tool
implicit none
real*8 :: om,oml,h
real*8 :: ra,dec,r,x,y,z,vx,vy,vz,vin(3),vout(3),rarange(2),decrange(2)
real*8 :: ramin,ramax,decmin,decmax,theta,rat,dect,zrange(2),rmin,rmax
integer :: imock,ngalaxy,i,j,k,n0,n1,n2,n3,n4
integer,allocatable :: icat(:,:)
character(len=200) :: ofile,ofile1,ofile2,ofile3,cmock
rarange=(/105d0,265d0/)
decrange=(/-5d0,60d0/)
zrange=(/0.43d0,0.70d0/)
ramin=rarange(1)
ramax=rarange(2)
decmin=decrange(1)
decmax=decrange(2)

om=Om_HR3
oml=OmL_HR3
h=h_HR3
rmin=comovingdistance(om,omL,zrange(1))
rmax=comovingdistance(om,omL,zrange(2))

imock=0
if(imock.eq.0) then
  cmock="00"
else if(imock.lt.10) then
  write(cmock,"(i1)") imock
  cmock="0"//trim(adjustl(cmock))
else
  write(cmock,"(i2)") imock
  cmock=trim(adjustl(cmock))
end if
ofile=trim(HR3Home)//"mock3"//trim(cmock)//"/radec.dat"
ofile1=trim(HR3Home)//"mock3"//trim(cmock)//"/radec1.dat"
ofile2=trim(HR3Home)//"mock3"//trim(cmock)//"/radec2.dat"
ofile3=trim(HR3Home)//"mock3"//trim(cmock)//"/radec3.dat"

call hr3_readpos(imock)
call hr3_readvel(imock)
ngalaxy=hr3mock(imock)%ngalaxy
allocate(icat(4,ngalaxy))
icat=0

open(10,file=ofile,status="replace")
open(11,file=ofile1,status="replace")
open(12,file=ofile2,status="replace")
do i=1,ngalaxy
   x=hr3mock(imock)%x(i)
   y=hr3mock(imock)%y(i)
   z=hr3mock(imock)%z(i)
   vx=hr3mock(imock)%vx(i)
   vy=hr3mock(imock)%vy(i)
   vz=hr3mock(imock)%vz(i)
   call Euclide2Sphere(x,y,z,r,ra,dec)
   if (r.ge.rmin.and.r.le.rmax) then
      vin(1)=x
      vin(2)=y
      vin(3)=z
      if (dec.ge.decmin.and.dec.le.decmax) then
         if (ra.ge.ramin.and.ra.le.ramax) then
            write(10,*) ra,dec
            icat(1,i)=1
         end if
         ra=ra+180d0
         if (ra.ge.360d0) ra=ra-360d0
         if (ra.ge.ramin.and.ra.le.ramax) then
            write(11,*) ra,dec
            icat(2,i)=1
         end if
      end if
      
      vin(3)=vin(3)*(-1d0)
      call rotation_z(vin,vout,-10d0)
      vin=vout
      call rotation_y(vin,vout,-40d0)
      call Euclide2Sphere(vout(1),vout(2),vout(3),r,ra,dec)
      if (ra.ge.ramin.and.ra.le.ramax.and.dec.ge.decmin.and.dec.le.decmax) then
         write(12,*) ra,dec
         icat(3,i)=1
      end if

      !call rotation_z(vin,vout,180d0)
      !vin=vout
      !call rotation_y(vin,vout,-30d0)
      !call Euclide2Sphere(vout(1),vout(2),vout(3),r,ra,dec)
      !if (ra.ge.ramin.and.ra.le.ramax.and.dec.ge.decmin.and.dec.le.decmax) then
      !   write(13,*) ra,dec
      !   icat(4,i)=1
      !end if
   end if
end do
close(10)
close(11)
close(12)

n0=0
n1=0
n2=0
n3=0
n4=0
do i=1,ngalaxy
   k=sum(icat(:,i))
   if (k.eq.0) n0=n0+1
   if (k.eq.1) n1=n1+1
   if (k.eq.2) n2=n2+1
   if (k.eq.3) n3=n3+1
   if (k.eq.4) n4=n4+1
end do
write(*,*) "n0=",n0
write(*,*) "n1=",n1
write(*,*) "n2=",n2
write(*,*) "n3=",n3
write(*,*) "n4=",n4


end
