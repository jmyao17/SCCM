       subroutine HN_Matrix_general(jmn,jmx,jdf,iq1,iq2,lpr)
       ! for a general case 
       use VAPHFB_Par
       implicit real*8 (a-h,o-z)
       logical lpr
   30 format('nqqjkk(',i2,2x,2i4,2x,2i4,')=',f15.8,2x,'hqqjkk=',f15.8)

       jqk(j,iq,k) = k+j+1+(2*j+1)*(iq-1) 

      do iis = jmn, jmx, jdf
         GCM%qkmax(iis) = (2*iis+1)*GCM%NOQ 
         if(iis.eq.1)  cycle
         do k1=-iis,iis, 1 !2
            iqk1 = jqk(iis,iq1,k1)  
            GCM%ik(iqk1,iis) = k1 
            GCM%iq(iqk1,iis) = iq1 
         do k2=-iis,iis, 1 !2
            iqk2 = jqk(iis,iq2,k2)  
            GCM%nmaxdi(iis)  = GCM%NOQ*(2*iis+1) 
            GCM%hjkkqq(iis,iqk1,iqk2) = dreal(Kernel%hjkk(iq1,iq2,iis,k1,k2))
            GCM%njkkqq(iis,iqk1,iqk2) = dreal(Kernel%njkk(iq1,iq2,iis,k1,k2))

            GCM%hjkkqq(iis,iqk2,iqk1) = GCM%hjkkqq(iis,iqk1,iqk2)
            GCM%njkkqq(iis,iqk2,iqk1) = GCM%njkkqq(iis,iqk1,iqk2)

         if(lpr .and. abs(GCM%njkkqq(iis,iqk1,iqk2)).gt.CHOP) then
            write(*,30) iis,k1,k2,iq1,iq2,GCM%njkkqq(iis,iqk1,iqk2),&
            GCM%hjkkqq(iis,iqk1,iqk2)/GCM%njkkqq(iis,iqk1,iqk2)
         endif
         enddo          
         enddo          
      enddo ! iis         
      return
       end
 




