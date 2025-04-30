subroutine filer_final()
use variables
implicit none

   
   open (11,file='information.dat',position='append')
   write(11,*) ' ' 
   write(11,*) '======================================================='
   write(11,*) 'total cost time = ',totalcosttime,'sec'
   write(11,*) 'final time = ',time,'sec'
   write(11,*) '======================================================='
   write(11,*) ' ' 
   close(11)

   close(21)
   close(31)
  
end subroutine filer_final

subroutine filerInfo()
use variables
implicit none
   open (1,file='information.dat',position='append')

   write(1,*) '!-----------------------------------------------------------!'
   write(1,*) '!  Finite Volume Method by projection method with mpi       !'
   write(1,*) '!-----------------------------------------------------------!'
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'number of processor(MPI) = ',nproc
   write(1,*) 'number of threads(OpenMP) = ',nthreads
   write(1,*) 'number of grid = ',nx*ny*nz/1000000.d0,'M'
   write(1,*) 'AOA = ',REAL(AOA)
   write(1,*) '                '
   write(1,*) 'x-component length = ',REAL(lx)
   write(1,*) 'y-component length = ',REAL(ly)
   write(1,*) 'z-component length = ',REAL(lz)
   write(1,*) '                '
   write(1,*) 'nx = ',nx
   write(1,*) 'ny = ',ny  
   write(1,*) 'nz = ',nz
   write(1,*) '                '
   write(1,*) 'dt = ',REAL(dt)
   write(1,*) 'CFL x = ', REAL(dt / dzsml), '  , CFL y = ', REAL(dt / dySml)
   write(1,*) 'Fo = ', dt/Re/dySml/dzSml
   write(1,*) '                ' 
   write(1,*) 'DFIB x_c = ',REAL(xc)
   write(1,*) 'DFIB yc = ',REAL(yc)
   write(1,*) 'DFIB zc = ',REAL(zc)
   write(1,*) '                '
   write(1,*) 'Re = ',REAL(Re)
   write(1,*) 'Free stream velocity = ', REAL(Freestream) ,' m/s '
   write(1,*) '                '
   if(Gridder=='non-uniform')then
   write(1,*) '====================non-uniform========================'
   write(1,*) 'dx = ',dx
   write(1,*) 'Large dy = ',dy
   write(1,*) 'Large dz = ',dz
   write(1,*) '                '
   write(1,*) 'nySml = ',nySml
   write(1,*) 'nzSml = ',nzSml
   write(1,*) '                '
   write(1,*) 'lyMid = ',lyMid
   write(1,*) 'lzMid = ',lzMid
   write(1,*) '                '
   write(1,*) 'Middle dy = ',dyMid
   write(1,*) 'Middle dz = ',dzMid
   write(1,*) '                '
   write(1,*) 'Small dy = ',dySml
   write(1,*) 'Small dz = ',dzSml
   write(1,*) '                '
   write(1,*) 'Small grid yc = ',GridderYc
   write(1,*) 'Small grid zc = ',GridderZc
   write(1,*) '                '
   write(1,*) 'lySml  = ',lySml
   write(1,*) 'lzSml  = ',lzSml
   write(1,*) '====================non-uniform========================'
   write(1,*) '                '
   end if
      if(Gridder=='non-uniform-sin2')then
   write(1,*) '====================non-uniform-sin2========================'
   write(1,*) 'dx = ',REAL(dx)
   write(1,*) 'Large dy = ',REAL(dy)
   write(1,*) 'Large dz = ',REAL(dz)
   write(1,*) '                '
   write(1,*) 'Middle dy = ',REAL(dyMid)
   write(1,*) 'Middle dz = ',REAL(dzMid)
   write(1,*) '                '
   write(1,*) 'Small dy = ',REAL(dySml)
   write(1,*) 'Small dz = ',REAL(dzSml)
   write(1,*) '                '
   write(1,*) 'Small grid yc = ',REAL(GridderYc)
   write(1,*) 'Small grid zc = ',REAL(GridderZc)
   write(1,*) '                '
   write(1,*) 'lyMid  = ',REAL(lyMid)
   write(1,*) 'nyMid = ',nyMid
   write(1,*) '                '
   write(1,*) 'lySml  = ',REAL(lySml)
   write(1,*) 'nySml = ',nySml
   write(1,*) '                '
   write(1,*) 'lzMid  = ',REAL(lzMid)
   write(1,*) 'nzMid = ',nzMid
   write(1,*) '                '
   write(1,*) 'lzSml  = ',REAL(lzSml)
   write(1,*) 'nzSml = ',nzSml
   write(1,*) '====================non-uniform-sin2========================'
   write(1,*) '                '
   end if
      if(Gridder=='non-uniform-sin')then
   write(1,*) '====================non-uniform-sin========================'
   write(1,*) 'Small grid yc = ',REAL(GridderYc)
   write(1,*) 'Small grid zc = ',REAL(GridderZc)
   write(1,*) '                '
   write(1,*) 'dx = ',REAL(dx)
   write(1,*) 'Large dy = ',REAL(dy)
   write(1,*) 'Large dz = ',REAL(dz)
   write(1,*) '                '
   write(1,*) 'Middle dy = ',REAL(dyMid)
   write(1,*) 'Middle dz = ',REAL(dzMid)
   write(1,*) '                '
   write(1,*) 'Small dy = ',REAL(dySml)
   write(1,*) 'Small dz = ',REAL(dzSml)
   write(1,*) '                '
   write(1,*) 'lyMid  = ',REAL(lyMid)
   write(1,*) 'nyMid = ',nyMid
   write(1,*) '                '
   write(1,*) 'lySml  = ',REAL(lySml)
   write(1,*) 'nySml = ',nySml
   write(1,*) '                '
   write(1,*) 'lzMid  = ',REAL(lzMid)
   write(1,*) 'nzMid = ',nzMid
   write(1,*) '                '
   write(1,*) 'lzSml  = ',REAL(lzSml)
   write(1,*) 'nzSml = ',nzSml
   write(1,*) '====================non-uniform-sin2========================'
   write(1,*) '                '
      end if
      if(Gridder=='uniform')then
   write(1,*) '=====================uniform==========================='
   write(1,*) 'dx = ',dx
   write(1,*) 'dy = ',dy
   write(1,*) 'dz = ',dz
   write(1,*) '=====================uniform==========================='
   write(1,*) '                '
   endif
   if (steadiness == 1) then; write(1,*) 'steadiness = steady'
   else if(steadiness == 2) then; write(1,*) 'steadiness = unsteady'
   end if
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   if (DBD == 1) then; write(1,*) 'DBD actuator on'
   else if(DBD == 0) then; write(1,*) 'DBD actuator off'
   end if
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   if (LES == 1) then; write(1,*) 'LES mode on'
   else if(LES == 0) then; write(1,*) 'LES mode off'
   end if
   write(1,*) 'LES Cs = ',REAL(Cs)
   write(1,*) 'Reduced Frequency k = ',REAL(reduce_frequency_k)
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'dt = ',dt
   write(1,*) 'max time step = ',nstep*dt,'sec'
   write(1,*) 'each time step = ',isto*dt,'sec'
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'p Residual =',REAL(zeta)
   write(1,*) 'velocity Residual =',REAL(zeta_vel)
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'isto =',isto
   write(1,*) 'istea =',istea
   write(1,*) '======================================================='
   close(1)



   open (21,file='CD_Time.dat',position='append')
   write(21,*) ' TITLE     = "" '
   write(21,*) ' VARIABLES = t*,AOA(<math>0</math>),C<sub>D</sub>,C<sub>L</sub> '
  

   !open (31,file='CD_AOA.dat',position='append')
   !write(31,*) ' TITLE     = "" '
   !write(31,*) ' VARIABLES = AOA,C<sub>D</sub>,C<sub>L</sub> '


   open (51,file='Cp.dat',position='append')
   write(51,*) ' TITLE     = "" '
   write(51,*) ' VARIABLES = x/c,C<sub>p</sub>'



end subroutine filerInfo



subroutine filerProcess()
use variables
implicit none

   write(21,*) REAL(time), REAL(AOA), REAL(cDrag), REAL(cLift)

   !write(31,*) REAL(AOA), REAL(cDrag), REAL(cLift)
   
end subroutine filerProcess



subroutine filerProcess_cp()
use variables
implicit none
real*8 :: x_c, p_avg, p_ref

   p_ref=0.d0
   do j=1,ny; do i=1,nx
      p_ref=p_ref+p(i,j,nz)*iDx(i)*iDy(j)
   enddo; enddo
   p_ref=p_ref/(lx*ly)

   write(51,*) 'ZONE solutiontime=',REAL(time)
   do k=1,nz
      p_avg=0.d0
      do j=1,ny
         
         if (ETA(5,j,k)>=0.5d0 .AND. ETA(5,j+1,k)<0.5d0) then
         
            ratio=(ETA(1,j+1,k)-0.5d0)/(ETA(1,j+1,k)-ETA(1,j,k))

            x_c= COS((AOA)*PI/180.d0)*(Zs(k)-zc-0.25d0) - SIN((AOA)*PI/180.d0)*(Ys(j)*ratio+Ys(j+1)*(1.d0-ratio)-yc) + 0.25d0  

            do i=1,nx
               p_avg=p_avg+p(i,j,k)*ratio+p(i,j+1,k)*(1.d0-ratio)
            enddo

            p_avg=p_avg/(nx*1.d0)-p_ref

            write(51,*) REAL(x_c), REAL(2.d0*p_avg)
                  
            exit

         endif
      enddo
   enddo

   do k=nz,1,-1
      p_avg=0.d0
      do j=1,ny
         
         if (ETA(5,j,k)<=0.5d0 .AND. ETA(5,j+1,k)>0.5d0) then
         
            ratio=(ETA(1,j+1,k)-0.5d0)/(ETA(1,j+1,k)-ETA(1,j,k))

            x_c= COS((AOA)*PI/180.d0)*(Zs(k)-zc-0.25d0) - SIN((AOA)*PI/180.d0)*(Ys(j)*ratio+Ys(j+1)*(1.d0-ratio)-yc) + 0.25d0  

            do i=1,nx
               p_avg=p_avg+p(i,j,k)*ratio+p(i,j+1,k)*(1.d0-ratio)
            enddo

            p_avg=p_avg/(nx*1.d0)-p_ref

            write(51,*) REAL(x_c), REAL(2.d0*p_avg)
                  
            exit

         endif
      enddo
   enddo

end subroutine filerProcess_cp



subroutine filereachtime()
use variables
implicit none

   !$OMP PARALLEL DO PRIVATE(i,j)  
   do k=1,nz; do j=1,ny; do i=1,nx
      uc(i,j,k) = 0.5d0*( u(i,j,k)+u(i-1,j,k) )
      vc(i,j,k) = 0.5d0*( v(i,j,k)+v(i,j-1,k) )
      wc(i,j,k) = 0.5d0*( w(i,j,k)+w(i,j,k-1) )
   enddo; enddo; enddo
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE(i,j)  
   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1) = p(i,j,k)
      Qout(i,j,k,2) = uc(i,j,k)
      Qout(i,j,k,3) = vc(i,j,k)
      Qout(i,j,k,4) = wc(i,j,k)
      Qout(i,j,k,5) = ETA(i,j,k)
   enddo; enddo; enddo
   !$OMP END PARALLEL DO


	write(filename,'(I4.4)')num
 	fileformat = '.q'

   open (17,file=TRIM(filename)//fileformat,position='append',form='unformatted')
   write(17) nblocks
   write(17) nx, ny, nz
   write(17) temp, temp, temp, REAL(time)
   write(17) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )

   close(17)

   num = num + 1



   !$OMP PARALLEL DO PRIVATE(i,j)  
   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1) = p(i,j,k)
      Qout(i,j,k,2) = u(i,j,k)
      Qout(i,j,k,3) = v(i,j,k)
      Qout(i,j,k,4) = w(i,j,k)
      Qout(i,j,k,5) = ETA(i,j,k)
   enddo; enddo; enddo
   !$OMP END PARALLEL DO


   if( istep < StartDynamic )then
      open (unit=18,form='unformatted',file='Static0.Q')
   else 
      open (unit=18,form='unformatted',file='OOXX.Q')
   endif
   write(18) nblocks
   write(18) nx, ny, nz
   write(18) temp, temp, temp, REAL(time)
   write(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )

   close(18)

end subroutine filereachtime




subroutine filer_bodyforce()
use variables
implicit none

   !$OMP PARALLEL DO PRIVATE(i,j)  
   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1)=edelta(i,j,k)
      Qout(i,j,k,2)=0.d0
      Qout(i,j,k,3)=F_tavey(i,j,k)
      Qout(i,j,k,4)=F_tavex(i,j,k)
      Qout(i,j,k,5)=ETA(i,j,k)
   enddo; enddo; enddo
   !$OMP END PARALLEL DO


   !open (unit=18,form='unformatted',file='PLASMA.Q')
   !write(filename,'(I4.4)')num
 	!fileformat = 'plasma.q'

   !open (18,file=TRIM(filename)//fileformat,position='append',form='unformatted')
   open (unit=18,form='unformatted',file='PLASMA.q')
   write(18) nblocks
   write(18) nx, ny, nz
   write(18) temp, temp, temp, REAL(time)
   write(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )

   close(18)


end subroutine filer_bodyforce


