*-----------------------------------------------------------------------
*.specification

      subroutine b2txcorrft (nCv, ns, dtim, switch, geo, mpg, ne, te,
     &  ti, sig0)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      use b2mod_b2cmpa
      use b2mod_switches
      use b2mod_user_namelist
     & , only : ft_omp, icsepomp
      use b2mod_numerics_namelist
     & , only : corr_core_safety, corr_core_ratio
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nCv, ns
      type (switches), intent(in) :: switch
      type (geometry), intent(in) :: geo
      type (mapping), intent(in) :: mpg
      real (kind=R8), intent(in) :: dtim, ne(nCv), te(nCv), ti(nCv),
     &  sig0(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TXCORRFT computes, for every flux tube lying on a closed flux
*     surface, the ratio of the stability limit of the artificial
*     poloidal averaging (used by b2upco and b2upht for the speed-up of
*     the pressure correction) to the time step,
*
*       delta_t_max = sigma_AN*(B*R)**2/(n_e*(T_e+T_i))
*       ratio       = delta_t_max/b2mndr_dtim
*
*     so that the weight of the local (non averaged) part of the
*     correction applied by b2upco/b2upht is
*
*       alpha(is) = min(1,corr_core_safety(is)*ratio)
*
*     with alpha=1 corresponding to the unmodified correction and
*     alpha -> 0 to the flux surface averaged one.  sigma_AN is the
*     artificial anomalous electrical conductivity of the potential
*     equation (sig0, as built by b2tqna with its barrier and SOL
*     enhancement), so that the margin is a local quantity and no
*     arbitrary radial profile is needed.
*     See Pressure_correction_speed-up.pdf in $SOLPSTOP/doc.
*
*     The result is stored in corr_core_ratio(nCv); a large value means
*     "no speed-up", i.e. alpha=1 for any safety factor.
*
*     All the input quantities are given by the caller, i.e. from the
*     driver (b2news_/b2news_m/b2mndt), right after the transport
*     coefficients, where co%sig0 is available.  The speed-up is active
*     only when at least one corr_core_safety is larger than zero;
*     otherwise the routine does nothing and corr_core_dn/corr_core_dt
*     are used as before.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFt, is
      integer, save :: ncall = 0
      real (kind=R8), parameter :: big = 1.0e30_R8
      real (kind=R8) ::
     &  vv(mpg%nFt), pp(mpg%nFt), aa(mpg%nFt), nn(mpg%nFt),
     &  bb(mpg%nFt), rr(mpg%nFt), ratio_ft(mpg%nFt)
      real (kind=R8) :: tl, al, almin, almax
*   ..procedures
      external xertst
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2txcorrft')
*   ..test nCv, ns
      call xertst (0.lt.nCv.and.0.lt.ns, 'faulty argument nCv, ns')
*   ..test that the module workspace is allocated
      if (.not.allocated(corr_core_ratio)) then
        call subend ()
        return
      endif

* ..no speed-up unless a safety factor is given
      corr_core_ratio = big
      if (maxval(corr_core_safety(0:ns-1)).le.0.0_R8) then
        ncall = ncall+1
        call subend ()
        return
      endif
*   ..test the time step
      if (dtim.le.0.0_R8) then
        if (ncall.eq.0) write(*,*) 'b2txcorrft: no time step '//
     &   'available, no speed-up of the pressure correction'
        ncall = ncall+1
        call subend ()
        return
      endif

* ..volume weighted sums over the cells on closed flux surfaces
      vv = 0.0_R8
      pp = 0.0_R8
      aa = 0.0_R8
      nn = 0.0_R8
      bb = 0.0_R8
      rr = 0.0_R8
      do iCv = 1, mpg%nCi
        if (.not.mpg%cvOnClosedSurface(iCv)) cycle
        iFt = mpg%cvFt(iCv)
        if (iFt.le.0.or.iFt.gt.mpg%nFt) cycle
        vv(iFt) = vv(iFt) + geo%cvVol(iCv)
        nn(iFt) = nn(iFt) + geo%cvVol(iCv)*ne(iCv)
*     ..te and ti are in Joule in the B2.5 state, so p = ne*(te+ti) is in Pa
        pp(iFt) = pp(iFt) + geo%cvVol(iCv)*ne(iCv)*
     &    (te(iCv)+ti(iCv))
        bb(iFt) = bb(iFt) + geo%cvVol(iCv)*
     &    (geo%cvBb(iCv,3)*geo%cvX(iCv))**2
        aa(iFt) = aa(iFt) + geo%cvVol(iCv)*sig0(iCv)*
     &    (geo%cvBb(iCv,3)*geo%cvX(iCv))**2
        rr(iFt) = rr(iFt) + geo%cvVol(iCv)*geo%cvX(iCv)
      enddo

* ..stability margin of the artificial poloidal averaging per tube
      do iFt = 1, mpg%nFt
        if (aa(iFt).gt.0.0_R8.and.pp(iFt).gt.0.0_R8) then
          ratio_ft(iFt) = (aa(iFt)/pp(iFt))/dtim
        else
          ratio_ft(iFt) = big
        endif
*   ..the outermost nrings surfaces are left without speed-up
        if (ft_omp(iFt).ge.icsepomp - switch%nrings)
     &    ratio_ft(iFt) = big
      enddo

* ..distribute the per-flux-tube margin over the cells
      do iCv = 1, mpg%nCi
        if (.not.mpg%cvOnClosedSurface(iCv)) cycle
        iFt = mpg%cvFt(iCv)
        if (iFt.gt.0.and.iFt.le.mpg%nFt)
     &    corr_core_ratio(iCv) = ratio_ft(iFt)
      enddo

* ..report the margins (switch corr_core_iout)
      if (switch%corr_core_iout.gt.0) then
        write(*,'(/a)') ' b2txcorrft: local stability limit of the'//
     &   ' pressure correction'
        write(*,'(a,1p,e11.3,a)') '  dtim=',dtim,' s'
        write(*,'(a)') '  CORR_CORE_SAFETY, species:'
        do is = 0, ns-1
          if (is_neutral(is)) cycle
          write(*,'(a,i3,a,1p,e10.2)') '    species ',is,' : ',
     &      corr_core_safety(is)
        enddo
        write(*,'(a)') '  tube    R_av   Te+Ti[eV]  ne[1e19]'//
     &   '   sig_AN[S/m]  B*R[Tm]  dtmax[s]    ratio  al_min  al_max'
        do iFt = 1, mpg%nFt
          if (vv(iFt).le.0.0_R8) cycle
          tl = pp(iFt)/(qe*nn(iFt))
          almin = big
          almax = 0.0_R8
          do is = 0, ns-1
            if (is_neutral(is)) cycle
            al = min(1.0_R8,corr_core_safety(is)*ratio_ft(iFt))
            almin = min(almin,al)
            almax = max(almax,al)
          enddo
          write(*,'(i5,f8.3,f11.1,f9.3,1p,4e11.3,2e9.2)')
     &      iFt, rr(iFt)/vv(iFt), tl, nn(iFt)/vv(iFt)*1.0e-19,
     &      aa(iFt)/bb(iFt), sqrt(bb(iFt)/vv(iFt)),
     &      aa(iFt)/pp(iFt), ratio_ft(iFt), almin, almax
        enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2txcorrft

      end subroutine b2txcorrft
