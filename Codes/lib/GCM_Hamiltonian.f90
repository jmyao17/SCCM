       subroutine HN_Matrix(jmn,jmx,jdf,iq1,iq2,icase,lpr)
       use VAPHFB_Par
       implicit real*8 (a-h,o-z)
       logical lpr

   30 format('nqqjkk(',2i2,2x,i2,2x,2i2,')=',f15.8,2x,'hqqjkk=',f15.8)

      if(icase.eq.0) then
!--- loop over total angular momentum 
      k1  = 0
      k2  = 0
      k1m = 0 
      k2m = 0
      do iis = jmn, jmx, jdf
         if(iis.eq.1) goto 119
         GCM%nmaxdi(iis)  = GCM%NOQ 
         GCM%hjkkqq(iis,iq1,iq2) = dreal(Kernel%hjkk(iq1,iq2,iis,k1,k2))
         GCM%njkkqq(iis,iq1,iq2) = dreal(Kernel%njkk(iq1,iq2,iis,k1,k2))
!------ symmetry in Kernels 
         GCM%hjkkqq(iis,iq2,iq1) = GCM%hjkkqq(iis,iq1,iq2)
         GCM%njkkqq(iis,iq2,iq1) = GCM%njkkqq(iis,iq1,iq2)

         if(lpr) then
            write(*,30) iq1,iq2,iis,k1m,k2m,GCM%njkkqq(iis,iq1,iq2),GCM%hjkkqq(iis,iq1,iq2)
         endif
  119 enddo ! iis         
      return
      endif 
!--------------------------------------------- 
 222   if(icase.eq.1) then
       do iis= jmn, jmx, jdf
          if(iis.eq.1) goto 120
!----- initialization  
          fac05 = one
          !fac3 = fac05*(2*iis+1)/(8*pi1**2)
                if(iv(iis).gt.zero) then       ! for even spin
                 GCM%nmaxdi(iis)  = (iis/2+1)*GCM%NOQ
                     iki0     = 0
              endif
                if(iv(iis).lt.zero) then       ! for odd spin
                 GCM%nmaxdi(iis)  = (iis-1)/2*GCM%NOQ
                      iki0    = 2
              endif
!     ..................................................................................
!      K1, K2 start from 0 or 2 and only the matrix elements with even K are non-zero
!      For even spin J, K1 and K2 start from iki0 =0
!      For odd  spin J, K1 and K2 start from iki0 =2
!      the mesh point i is given by: i = iq + (K-iki0)/2*maxmp  
!     ..................................................................................
!      For example:
!**************************************************************************************************
!      1)  J=4, iki0=0
!          K=   iki0,                     2,                               4
!          i= 1,2,...,maxmp,  maxmp+1,maxmp+2,...,2*maxmp,  2*maxmp+1,2*maxmp+2,...,2*maxmp+maxmp
!      2)  J=7, iki0=2
!          K=   iki0,                     4,                               6
!          i= 1,2,...,maxmp,  maxmp+1,maxmp+2,...,2*maxmp,  2*maxmp+1,2*maxmp+2,...,2*maxmp+maxmp
!**************************************************************************************************
              do k1=iki0,iis,2   ! K-value
              do k2=iki0,iis,2
                   k1m = k1-iki0   ! mesh point
                   k2m = k2-iki0
                 GCM%hjkkqq(iis,iq1+k1m/2*GCM%NOQ,iq2+k2m/2*GCM%NOQ)  &
     &           = dreal(Kernel%hjkk(iq1,iq2,iis,k1,k2)) !
                 GCM%njkkqq(iis,iq1+k1m/2*GCM%NOQ,iq2+k2m/2*GCM%NOQ)  &
     &           = dreal(Kernel%njkk(iq1,iq2,iis,k1,k2)) !nkk(iis,k1,k2)
              enddo ! k2
              enddo ! k1 
!---- symmetry in kernels 
             do k1=iki0,iis,2
             do k2=iki0,iis,2
                  k1m = k1-iki0   ! mesh point
                  k2m = k2-iki0
                GCM%hjkkqq(iis,iq2+k2m/2*GCM%NOQ,iq1+k1m/2*GCM%NOQ)    &
     &           = dreal(Kernel%hjkk(iq1,iq2,iis,k1,k2)) !hkk(iis,k1,k2)
                GCM%njkkqq(iis,iq2+k2m/2*GCM%NOQ,iq1+k1m/2*GCM%NOQ)    &
     &           = dreal(Kernel%njkk(iq1,iq2,iis,k1,k2)) !nkk(iis,k1,k2)
              enddo !
              enddo ! 
!---------- output
       if(lpr) then
           do k1=iki0,iis,2
           do k2=iki0,iis,2
                k1m = k1-iki0   ! mesh point
                k2m = k2-iki0
                write(*,30) iq1,iq2,iis,k1m,k2m,  &
     &          GCM%njkkqq(iis,iq2+k2m/2*GCM%NOQ,iq1+k1m/2*GCM%NOQ), &
     &          GCM%hjkkqq(iis,iq2+k2m/2*GCM%NOQ,iq1+k1m/2*GCM%NOQ)     

           enddo !
           enddo !
       endif
  120 enddo ! iis
      endif

       return
       end
 




