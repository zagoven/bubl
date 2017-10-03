module methan
 real,dimension(1000):: z,t,s,d,zm,cm !height,temperature,salinity 
 real,parameter:: g=981,rgaz=82.0575,rc=0.0584 !gravity constant, individual gas constant, critical radius
 real,parameter:: pi=3.14256,dif=1.5e-05,hg=7.9e+05 !pi, diffusion, Henry law constant
 real,parameter:: pa=1.
 integer:: nz,nm
 real:: kb,bb,tt ! k mass transfer
end module
!
program methane
 use methan
 implicit none
 real:: fr,mt,vb
 real:: rstp,time,aa
 real,dimension(500000):: frr,mtt,ft
 real:: rn,rp,dt,zr,bot,mn,mp,rnn,mnn
 real:: k1,k2,k3,k4,fnn
 integer:: n,i
! Reading z,t,s
 open(20,file='mar.dat')
 open(21,file='BlackSea.dat')
 open(22,file='mar_3500.dat')
 nz=0
 do while(.not.eof(20))
  nz=nz+1
! Reading depth(horizont)(meters), temperature (Celsium degrees)
  read(20,*) z(nz),t(nz),s(nz)
  d(nz)=rstp(s(nz),t(nz),0.,aa)/1000. ! in gr/cm^3
 end do
 z=z*1.e2 ! in cm 
 
! Counter for aqua concentrarion of methane 
! It is better to add the same for oxygen data
 nm=0 
 do while(.not.eof(21))
  nm=nm+1
  read(21,*) zm(nm),cm(nm) ! Nanomol ! reading the horizont and value of aqua concentrarion on this horizont
 end do
 zm=zm*1.e2   ! depth (horizonts) to cm
 cm=cm/1.e9   ! concentration to mol
 
!Initial radius of bubble
 rp=3500. ! r in micrometer
!Initial depth
 bot=1.e04 ! in cm
!The step in sec 
 dt=0.1 
 zr=bot 
 
!Initial volume of methane ( where r=0)
mp=4.8e-04  ! for rp=3500. 
!mp=2.4e-03 !for rp=6000. 
!mp=2.4e-05  
!mp=1.2e-05 !for rp=1000.

! Time cycle
 n=0; time=0.
 do while(.true.)
  n=n+1
  time=time+dt
  if(n<4) then
      
! Radius calculation 
! With Adams–Bashforth/Moulton method 4 order

   frr(n)=fr(zr,rp)
   zr=zr-dt*vb(rp)   
   rn=rp+dt*frr(n)
   mn=mp+4.*pi*rp*rp*bb*1.e-8*dt  
   rp=rn
   mp=mn
   
  else

   zr=zr-dt*vb(rp)
   frr(n)=fr(zr,rp)
   rnn=rp+dt*(55.*frr(n)-59.*frr(n-1)+37.*frr(n-2)-9.*frr(n-3))/24.
   zr=zr-dt*vb(rnn)
   fnn=fr(zr,rnn)
   rn=rp+dt*(9.*fnn+19*frr(n)-5.*frr(n-1)+frr(n-2))/24.
   mn=mp+4.*pi*rp*rp*bb*1.e-8*dt
   rp=rn
   mp=mn
  end if

  if(zr<100.) exit
  write(22,'(4e13.5,\)') time,zr/100.,rn/1000.,mn, kb
 end do    

    end program
    
real function fr(zr,r)
 use methan
 implicit none
 real:: r,re,v,rcm,sur,zr !radius,reinolds, rcm - radius in cm
 real:: dvisc,surften,vb
 integer:: i
 real:: ss,dd,aa,cc,vbb
 real:: a1,a2,a3,a4,pb,dm
! Interpolation t,s to zr
 do i=1,nz-1
  if(zr>=z(i).and.zr<z(i+1)) then
   aa=(zr-z(i))/(z(i+1)-z(i))
   tt=t(i)+(t(i+1)-t(i))*aa
   ss=s(i)+(s(i+1)-s(i))*aa
   dd=d(i)+(d(i+1)-d(i))*aa
   exit
  end if
 end do
!
 do i=1,nm-1
  if(zr>=zm(i).and.zr<zm(i+1)) then
   cc=cm(i)+(cm(i+1)-cm(i))*(zr-zm(i))/(zm(i+1)-zm(i))
   exit
  end if
 end do 
! Velocity of bubble rising
 vbb=vb(r)
! 
 v=dvisc(ss,tt,0.)/(dd*1000.)*1.e04 ! cm^2/sec
 rcm=r*1.e-04
 re=2.*rcm*vbb/v ! Re, d/less
 kb=sqrt(2/pi*(1-2.89*re**-0.5)*dif*vbb/rcm) ! cm/sec
! mean density
 dm=0.
 do i=1,nz
  if(zr>z(i)) then
   dm=dm+(d(i+1)+d(i))*0.5*(z(i+1)-z(i))
  else
   exit
  end if
 end do
 dm=dm/zr
 sur=surften(ss,tt)/(r*1.e-06)
 pb=pa+(0.1*g*dm*zr+2.*sur)*9.87167e-6
 a1=3.*(pa*101325.+dm*1000.*g/100.*zr/100.)+4.*sur
 a2=r*1.e-06*dm*1000.*g/100.*vbb/100.
 bb=kb*(cc-pb/hg)
 a3=3*rgaz*(tt+273.15)*bb/100.*101325.
 fr=(a2+a3)/a1*1.e6
 return
end function
!
real function mt(r)
 use methan
 implicit none
 real:: r,rcm
 rcm=r*1.e-04
 mt=bb*4.*3.14256*rcm*rcm
 return
    end function mt

    ! Selection and calculation of bubble velocity
    
real function vb(r)
 use methan
 implicit none
 real:: r
 real:: m1,m2,j1,j2,vbmin 
  if (r<4.e03)  then
  m1 = -0.849
  m2 = -0.815
  j1 = 0.733
  j2 = 4.792E-4
  vbmin = 22.16
 else if (r>=4.e03.and.r<4.e04) then 
  m1 = 0.0
  m2 = 0.0
  j1 = 11.05
  j2 = 0.0
  vbmin = 19.15
 else
  print *, 'Incorrect radius'
  stop
 end if
 vb=(vbmin+j1*(r-rc)**m1)*exp(j2*tt*(r-rc)**m2) ! cm/sec
 return
end function