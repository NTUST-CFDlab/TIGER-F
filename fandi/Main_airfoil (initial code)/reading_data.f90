subroutine reading_data()
   use variables
   implicit none

    
   if(Gridder=='non-uniform')then
      !Small interval
      dySml = lySml/nySml
      dzSml = lzSml/nzSml

      !Middle intervel
      dyMid = (lyMid-lySml)/(nyMid-nySml)
      dzMid = (lzMid-lzSml)/(nzMid-nzSml)

      !Large interval
      dy = ( ly-lyMid ) / ( ny - nyMid )
      dz = ( lz-lzMid ) / ( nz - nzMid )

       

      dx = lx / nx 
    

   else if(Gridder=='uniform')then

      dx = lx / nx 
      dy = ly / ny
      dz = lz / nz

   else if(Gridder=='non-uniform-sin')then
      !Small interval
      dySml = lySml/nySml
      dzSml = lzSml/nzSml

      dyMid = (lyMid-lySml)/(nyMid-nySml)
      dzMid = (lzMid-lzSml)/(nzMid-nzSml)
      !Large interval
      dy = ( ly-lySml ) / ( ny - nySml )
      dz = ( lz-lzMid ) / ( nz - nzMid )

      dx = lx / nx 

   else if(Gridder=='non-uniform-sin2')then
      !Small interval
      dySml = lySml/nySml
      dzSml = lzSml/nzSml

      dyMid = (lyMid-lySml)/(nyMid-nySml)
      dzMid = (lzMid-lzSml)/(nzMid-nzSml)
      !Large interval
      dy = ( ly-lyMid ) / ( ny - nyMid )
      dz = ( lz-lzMid ) / ( nz - nzMid )

      dx = lx / nx 
        
   end if



   nu = DBLE(1.d0/Re)*1.d0

end subroutine reading_data






subroutine reading_variables()
use variables
implicit none

open (18,file=inputfile,form='unformatted')
read(18) inblocks
read(18) inx, iny, inz
read(18) temp, temp, temp, temp
read(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
close(18)


   do k=1,nz; do j=1,ny; do i=1,nx  

      p(i,j,k) = Qout(i,j,k,1)
      p_no(i,j,k) = Qout(i,j,k,1)
      p_pre(i,j,k,1) = Qout(i,j,k,1)   

   end do; end do; end do  
   

   do k=1,nz; do j=1,ny; do i=1,nx  

      u(i,j,k) = Qout(i,j,k,2)
      v(i,j,k) = Qout(i,j,k,3)
      w(i,j,k) = Qout(i,j,k,4)  

   enddo; enddo; enddo 


   temp=0.

end subroutine reading_variables


