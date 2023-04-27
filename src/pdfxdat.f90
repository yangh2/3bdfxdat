!!c     Calculates partial pair distribution functions
!!c     Requires VASP5-format XDATCAR.
!!c     Run with cat XDATCAR.* | pdfxdat
implicit none

character*75 species,line1,line2
character ( len = 8) filename
character*2 Chem(10)
character*240 lmnts
data lmnts /"H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRgCn"/
real*8 BB(96)/-3.739,3.26,-1.9,7.79,5.3,6.646,9.36,5.803, &       !!  1- 8
       5.654,4.566,3.63,5.375,3.449,4.1491,5.13,2.847,9.577, &    !!  9-17
       1.909,3.67,4.7,12.29,-3.438,-0.3824,3.635,-3.73,9.45, &    !! 18-26
       2.49,10.3,7.718,5.68,7.288,8.185,6.58,7.97,6.795,7.81, &   !! 27-36
       7.09,7.02,7.75,7.16,7.054,6.715,0,7.03,5.88,5.91,5.922, &  !! 37-47
       4.87,4.065,6.225,5.57,5.8,5.28,4.92,5.42,5.07,8.24,4.84, & !! 48-58
       4.58,7.69,0,0.8,7.22,6.5,7.38,16.9,8.01,7.79,7.07,12.43, & !! 59-70
       7.21,7.7,6.91,4.86,9.2,10.7,10.6,9.6,7.63,12.692,8.776, &  !! 71-81
       9.405,8.532,0,0,0,0,0,0,10.31,0,8.417,0,0,0,0/             !! 82-96

integer natoms,nsteps,ntype,numb(106)
integer i,j,k
integer isep,hist(2000,10,10)
integer alim,blim,clim,idr
integer filenumber,N12(10,10),ns,t1,t2,la,lb,lc
integer vlabc(3)
integer typeZ(10)
real*8 abc(3,3)
real*8 xyz(10,10000,3), dxyz(3)  !! 10 species 10000 atoms 3 coordinates
real*8 dij,pdf,g,sigma,volume,dx,dy,dz,r
real*8 pdfX(2000),pdfN(2000)
real*8 Xfactor(10,10),Nfactor(10,10),C(10),Zbar,Bbar

integer Count_Items
real*8 gauss,triple

filename="pdf.xyzw"
!!$ bin width for histogram.
sigma=0.025
read(*,'(a)') line1
read(*,*)
read(*,*) abc(1,:)
read(*,*) abc(2,:)
read(*,*) abc(3,:)
!!$ Try to guess limits for repeated images.
!!$ Improved efficiency possible from minimum image method
!!$ if size is greater than 20 Angstrom.
alim=1+20./sqrt(dot_product(abc(1,:),abc(1,:)))
blim=1+20./sqrt(dot_product(abc(2,:),abc(2,:)))
clim=1+20./sqrt(dot_product(abc(3,:),abc(3,:)))

!!$ Read header of first file
volume=triple(abc(1,:),abc(2,:),abc(3,:))
read(*,'(A)') species
ntype=Count_Items(species)
if(ntype.gt.10) then
   write(*,*) 'Too many types ',ntype
   stop
endif
read(*,*) (numb(i),i=1,ntype)
natoms=sum(numb)
if(natoms.gt.10000) then
   write(*,*) 'Too many atoms ',natoms
   stop
endif

!!$ Read in types and determine x-ray/neutron scale factors
do t1=1,ntype
   call s_word_extract_first(species,Chem(t1))
   typeZ(t1)=(1+index(lmnts,Chem(t1)))/2
   C(t1)=float(numb(t1))/float(natoms)
   Zbar=Zbar+C(t1)*float(typeZ(t1))
   Bbar=Bbar+C(t1)*BB(typeZ(t1))
enddo
do t1=1,ntype
   do t2=1,ntype
      N12(t1,t2)=numb(t1)*numb(t2)
      if(t1.eq.t2) N12(t1,t2)=numb(t1)*(numb(t2)-1)
      Xfactor(t1,t2) = C(t1)*C(t2)*float(typeZ(t1))*float(typeZ(t2))
      Nfactor(t1,t2) = C(t1)*C(t2)*BB(typeZ(t1))*BB(typeZ(t2))
   enddo
enddo

nsteps=0
!!$ Start loop over XDATCAR files
do
!!$ Detect start of new file
   read(*,'(a)',end=1) line2
   if(line1 == line2) then
      do k=1,7
         read(*,*)
      enddo
   endif
   nsteps=nsteps+1

!!$ Read structure from file
   do t1=1,ntype
      do i=1,numb(t1)
         read(*,*) xyz(t1,i,:)
      enddo
   enddo

!!$ Loop over type1, type2, repeated images
   do t1=1,ntype
      do t2=1,ntype
         do i=1,numb(t1)
            do j=1,numb(t2)
               do la=-alim,alim
                  do lb=-blim,blim
                     do lc=-clim,clim
                        dxyz=xyz(t1,i,:)-xyz(t2,j,:)+(/la,lb,lc/)
                        dxyz=matmul(dxyz,abc)
                        dij=sqrt(dot_product(dxyz,dxyz))
                        isep=100.0*dij
                        if(isep.le.2000) then
                           hist(isep,t1,t2)=hist(isep,t1,t2)+1
                           if((t1.eq.t2).and.(i.eq.j)) then
                              hist(isep,t1,t2)=hist(isep,t1,t2)-1
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
enddo

1 continue
!!$ Done with loop over XDATCAR files

!!$ Post-processing smearing, x-ray/neutron weighting and file output
do t1=1,ntype
   do t2=1,ntype
      filename(5:6)=Chem(t1)
      filename(5+len_trim(Chem(t1)):)=Chem(t2)
      filenumber=10*t1+t2
      open(filenumber,file=filename)
      open(101,file="pdf.xray")
      open(102,file="pdf.neutron")
      do idr=1,2000
         r=float(idr)/100.0
         pdf=0.0
         do isep=1,2000
            g=gauss((r-(0.5+float(isep))/100.0)/(sigma))/sigma
            pdf=pdf+float(hist(isep,t1,t2))*g/float(nsteps)
         enddo
         pdf=volume*pdf/(4.0*3.141592654*float(N12(t1,t2))*r*r)
         pdfX(idr)=pdfX(idr)+Xfactor(t1,t2)*pdf/Zbar/Zbar
         pdfN(idr)=pdfN(idr)+Nfactor(t1,t2)*pdf/Bbar/Bbar
         write(filenumber,900) r,pdf
      enddo
   enddo
enddo
do idr=1,2000
   r=float(idr)/100.0
   write(101,900) r,pdfX(idr)
   write(102,900) r,pdfN(idr)
enddo
stop
900 format(2f12.6)
end program

! --------------------
INTEGER FUNCTION Count_Items(s)  ! in string or C string that are blank or comma separated
  CHARACTER(*) :: s
  INTEGER :: i, k

  k = 0  ; IF (s /= ' ') k = 1      ! string has at least 1 item
  DO i = 1,LEN_TRIM(s)-1
     IF (s(i:i) /= ' '.AND.s(i:i) /= ',' &
          .AND.s(i+1:i+1) == ' '.OR.s(i+1:i+1) == ',') k = k+1
  END DO
  Count_Items = k
END FUNCTION Count_Items


subroutine s_word_extract_first ( s, w )
  implicit none
  integer   ( kind = 4 ) get1
  integer   ( kind = 4 ) get2
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character ( len = * )  w
  w = ' '
  s_length = len_trim ( s )
  if ( s_length < 1 ) then
    return
  end if
!  Find the first nonblank.
  get1 = 0
  do
    get1 = get1 + 1
    if ( s_length < get1 ) then
      return
    end if
    if ( s(get1:get1) /= ' ' ) then
      exit
    end if
  end do
!  Look for the last contiguous nonblank.
  get2 = get1
  do
    if ( s_length <= get2 ) then
      exit
    end if
    if ( s(get2+1:get2+1) == ' ' ) then
      exit
    end if
    get2 = get2 + 1
  end do
!  Copy the word.
  w = s(get1:get2)
!  Shift the string.
  s(1:get2) = ' '
  s = adjustl ( s )
  return
end

real*8 function gauss(x)
  real*8 x
  if(abs(x).gt.4.) then
     gauss=0.
     return
  else
     gauss=(1.0/sqrt(2.0*3.141592654))*exp(-x*x/2.0)
  endif
  return
end function gauss

real*8 function triple(a,b,c)
  real*8 a(3),b(3),c(3),d(3)
  d(1)=a(2)*b(3)-a(3)*b(2)
  d(2)=a(3)*b(1)-a(1)*b(3)
  d(3)=a(1)*b(2)-a(2)*b(1)
  triple=dot_product(c,d)
  return
end function triple
