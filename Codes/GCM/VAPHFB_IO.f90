        subroutine print_matrix(LDA,A_matrix)
        implicit none
        integer LDA,ii,jj
        real*8  A_matrix(LDA,LDA)

        do ii=1,LDA
        do jj=1,LDA
           write(911,'(2i6,f12.8)') jj,ii,A_matrix(jj,ii)
        enddo
        enddo
        stop
        END SUBROUTINE

        subroutine print_zmatrix(LDA,ZA_matrix)
        implicit none
        integer     LDA,ii,jj
        complex*16  ZA_matrix(LDA,LDA)

        do ii=1,LDA
        do jj=1,LDA
           write(911,'(2i6,2e12.3)') jj,ii,ZA_matrix(jj,ii)
        enddo
        enddo
        stop
        END SUBROUTINE

        subroutine print_array(LDA,V,istop)
        implicit none
        integer     LDA,ii,istop
        real*8 V(LDA)

        do ii=1,LDA
           write(911,'(i6,e12.3)') ii,V(ii)
        enddo
        if(istop .ne.0) stop

        END SUBROUTINE

      subroutine ReadCoordinators()
      use vaphfb_par
      implicit none
      integer NQ,iq
!    .......... read coordinates
         open(15,file='betgam.dat',status='old')
         read(15,*)
         read(15,*)
         read(15,*)
         read(15,*)
!     ............
         read(15,*) NQ
         print *, NQ
         stop
         if(.NOT. ALLOCATED(Const%beta2t_mesh))   ALLOCATE(Const%beta2t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%gamma2t_mesh))  ALLOCATE(Const%gamma2t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%P00_mesh))      ALLOCATE(Const%P00_mesh(1:NQ))
         do iq=1,NQ
            read(15,*) Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),Const%P00_mesh(iq)
         enddo !
!         Input%NGCM=NQ
         GCM%NOQ   =NQ

      end
!________________________________________________________________________
      subroutine Filename4wfs(betat,gammat,p00)
!.......................................................................
      USE VAPHFB_PAR
      implicit none
      character*1 sign1,sign2
      character*6 name0
      real*8      betat,gammat,p00,ab2c1
      integer     nuc1,nuc2,name1,name2,name3,&
     &            name4,name5,namep1,namep2,namep3
!---------------------------
              name0 = '../../'
              if(betat.ge.0.d0) sign1='+'
              if(betat.lt.0.d0) sign1='-'
              ab2c1  = abs(betat)
              name1 = ab2c1+48
              name2 = mod(ab2c1*10,10.d0)  +48
              name3 = mod(ab2c1*100,10.d0) +48
              name4 = mod(gammat/10,10.d0)+48
              name5 = mod(gammat,10.d0)   +48
              nuc1 = mod(Nucl%nucleon(2)/10,10) +48
              nuc2 = mod(Nucl%nucleon(2),10)    +48
              namep1 = p00+48
              namep2 = mod(p00*10,10.d0)  +48
              namep3 = mod(p00*100,10.d0) +48
!---------------------------  
!      if(PNP%NFOM.gt.1) then
!      File%wf =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)    &
!     &           //'_VAP_'//Input%cIntID       &
!     &           //'_beta'//sign1//char(name1)                       &
!     &          //char(name2)//char(name3)                           &
!     &          //'gam'//char(name4)//char(name5)                    &
!     &          //'np'//char(namep1)//char(namep2)//'.wf' !      
!
!      else
      File%wf =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)    &
     &           //'_HFB_'//Input%cIntID       &
     &           //'_beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &          //'np'//char(namep1)//char(namep2)//'.wf' !      
!      endif

!---------------------------  
      if(PNP%NFOM.gt.1) then
      File%out =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)   &
     &           //'_VAP_'//Input%cIntID       &
     &           //'_beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &       //'np'//char(namep1)//char(namep2)//'.out' !      

      else
      File%out =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)   &
     &           //'_HFB_'//Input%cIntID       &
     &           //'_beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &        //'gam'//char(name4)//char(name5)                    &
     &       //'np'//char(namep1)//char(namep2)//'.out' !      
      endif

      return
      end





!=======================================================================
      subroutine nucleus(is,npro,te)

!=======================================================================
!
!     is = 1 determines the symbol for a given proton number npro
!          2 determines the proton number for a given symbol te
!
!-----------------------------------------------------------------------
!
      PARAMETER (MAXZ=140)
!
      CHARACTER TE*2,T*(2*MAXZ+2)
!
      T(  1: 40) = '  _HHeLiBe_B_C_N_O_FNeNaMgAlSi_P_SClAr_K'
      T( 41: 80) = 'CaSsTi_VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr_Y'
      T( 81:120) = 'ZrNbMoTcRuRhPdAgCdInSnSbTe_IXeCsBaLaCePr'
      T(121:160) = 'NdPmSmEuGdTbDyHoErTmYbLuHfTa_WReOsIrPtAu'
      T(161:200) = 'HgTlPbBiPoAtRnFrRaAcThPa_UNpPuAmCmBkCfEs'
      T(201:240) = 'FmMdNoLrRfHaSgNsHsMr10111213141516171819'
      T(241:280) = '2021222324252627282930313233343536373839'
      T(281:282) = '40'
!
! ... Rf is called also as Ku (kurchatovium)
! ... Ha: IUPAC calls it as dubnium (Db). J.Chem.Educ. 1997, 74, 1258
! ... Ha is called also as Db (Dubnium)
!
      if (is.eq.1) then
         if (npro.lt.0.or.npro.gt.maxz) stop 'in NUCLEUS: npro wrong'
         te = t(2*npro+1:2*npro+2)
         return
      else
!
         do np = 0,maxz
            if (te.eq.t(2*np+1:2*np+2)) then
               npro = np
               if (npro.gt.maxz) write(6,100) TE
               return
            endif
         enddo
!
         write(6,100) TE
  100    format(//,' NUCLEUS ',A2,'  UNKNOWN')
      endif
!
      stop
      END
!
