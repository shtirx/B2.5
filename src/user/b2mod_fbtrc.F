      module b2mod_fbtrc
      use b2mod_types
      use b2mod_diag
      use b2mod_time
      use b2us_feedback, only: fb_current, fb_target, fb_rescale,
     &  saved_fb_actuator, fb_type, fb_rescale_option, fb_actuator,
     &  fb_alpha, fb_beta, fb_const, nFb, fb_species
      implicit none

c========================================================================
c*** Writes per-iteration feedback control data to
c*** "tracing/fbcontrol.trc" (added 2026-06-03)
c***
c*** This trace captures the state of the feedback control system
c*** (b2.feedback_control.parameters) at each outer iteration:
c***   - fb_curr_   : measured quantity (e.g. separatrix density)
c***   - fb_targ_   : target value
c***   - fb_resc_   : rescale factor (or absolute puff for option 9/10)
c***   - fb_actv_   : applied actuator value (puff rate or BC flux)
c***   - fb_type_   : feedback type number (constant from namelist)
c***   - fb_rsclopt : rescale option number (constant)
c***   - fb_actu_   : actuator number (constant)
c***   - fb_alpha_  : alpha parameter (constant)
c***   - fb_beta_   : beta parameter (constant)
c***   - fb_const_  : const/plateau parameter (constant)
c***
c*** One line per feedback per iteration. If NFB>1 in
c*** b2.feedback_control.parameters, all columns repeat for
c*** feedback #2, #3, etc. (appended to the same line).
c***
c*** Controlled by ank_tracing and filedata(12) in b2.user.parameters.
c*** If nFb=0 (no feedback defined), the file is not created.
c-----------------------------------------------------------------------
c
c*** lv    :  version number to appear in the file
c
      integer lv
      parameter (lv=20260603)
      integer ltrc
      parameter (ltrc=ntrc_fbc)
      integer, save :: lun

c-----------------------------------------------------------------------
c*** Feedback quantities (per active feedback scheme):
c***
c***    fb_curr_   fb_current         measured value
c***    fb_targ_   fb_target          target value
c***    fb_resc_   fb_rescale         rescale factor (or absolute puff)
c***    fb_actv_   saved_fb_actuator  applied puff / BC value
c***    fb_type_   fb_type            feedback type (constant)
c***    fb_rsclopt fb_rescale_option  rescale option (constant)
c***    fb_actu_   fb_actuator        actuator type (constant)
c***    fb_alpha_  fb_alpha           rescaling alpha (constant)
c***    fb_beta_   fb_beta            rescaling beta (constant)
c***    fb_const_  fb_const           rescaling const / plateau (constant)
c-----------------------------------------------------------------------

      integer, save :: nfb_tmp
      integer nq_tmp
      parameter (nq_tmp=10)
      character, save :: f*(lenfnm), hlp_frm*512
      character, save :: hlp_hdr*256

      character, save :: hlp_hq(nq_tmp)*12
      data hlp_hq / 'fb_curr_', 'fb_targ_', 'fb_resc_', 'fb_actv_',
     &              'fb_type_', 'fb_rsclopt','fb_actu_',
     &              'fb_alpha_','fb_beta_', 'fb_const_' /

      contains

      subroutine b2fbtrc()
      implicit none
c-----------------------------------------------------------------------
c*** Write one line of feedback data
c-----------------------------------------------------------------------
      integer ifb

      if(filedata(ltrc)) then
        if(nfb_tmp.gt.0) then
          write(lun,hlp_frm,err=94) tim,
     &       (fb_current(ifb),         ifb=1,nfb_tmp),
     &       (fb_target(ifb),          ifb=1,nfb_tmp),
     &       (fb_rescale(ifb),         ifb=1,nfb_tmp),
     &       (saved_fb_actuator(ifb),  ifb=1,nfb_tmp),
     &       (dble(fb_type(ifb)),            ifb=1,nfb_tmp),
     &       (dble(fb_rescale_option(ifb)),  ifb=1,nfb_tmp),
     &       (dble(fb_actuator(ifb)),        ifb=1,nfb_tmp),
     &       (fb_alpha(ifb),           ifb=1,nfb_tmp),
     &       (fb_beta(ifb),            ifb=1,nfb_tmp),
     &       (fb_const(ifb),           ifb=1,nfb_tmp)
        end if
      end if
      return

   94 write(*,*) 'b2fbtrc: cannot write in the file ',trim(f)
      filedata(ltrc)=.false.
      return
      end subroutine b2fbtrc
c=======================================================================
c
      subroutine b2fbtrci(switch)
      use b2mod_switches
      implicit none
c-----------------------------------------------------------------------
c*** Initialise the file for feedback tracing data
c-----------------------------------------------------------------------
      type(switches), intent(in) :: switch
      integer i, j, ifb, nv, nfb_old, nq_old
      character hlp_c
      logical ex, warn_tmp

      if (nFb.le.0) then
        filedata(ltrc) = .false.
        return
      end if

      nfb_tmp = nFb

      lun = ntrc(ltrc)
      f = fn_ntrc(ltrc)

      call find_file(f, ex)
      open(lun, file=f, err=90)
      rewind(lun)
      if (ex) read(lun, *, end=91) hlp_c
      goto 93
   91 ex = .false.
   93 continue
      rewind(lun)

      if (ex .and. switch%b2mndr_stim .lt. 0.0_R8) then
c***   Old file exists and restart is requested — check consistency
        read(lun, *, err=92) hlp_c, nv, nfb_old, nq_old
        warn_tmp = .false.

        if (nv .ne. lv) then
          write(*,'(a,i8,3a,i8)')
     &      '!!! b2fbtrci: version mismatch, current=', lv,
     &      ' file=', trim(f), ' version=', nv
          filedata(ltrc) = .false.
          close(lun)
          return
        end if

        if (nfb_old .ne. nfb_tmp) then
          write(*,'(a,i3,3a,i3)')
     &      '!!! b2fbtrci: nFb mismatch, current=', nfb_tmp,
     &      ' file=', trim(f), ' nFb=', nfb_old
          filedata(ltrc) = .false.
          close(lun)
          return
        end if

        if (nq_old .ne. nq_tmp) then
          write(*,'(a,i3,3a,i3)')
     &      '!!! b2fbtrci: nq mismatch, current=', nq_tmp,
     &      ' file=', trim(f), ' nq=', nq_old
          filedata(ltrc) = .false.
          close(lun)
          return
        end if

        if (warn_tmp) write(*,*) ' You should remove the file!'

c***   Read to end of file for append
   50   read(lun, '(a1)', end=55) hlp_c
        go to 50
   55   backspace lun

c***   Build format strings (already computed for writing)
        if (nfb_tmp .gt. 0) then
          write(hlp_frm, '(a,i3,a)') '(1p,e16.7,',
     &      nq_tmp * nfb_tmp,
     &      '(2x,e12.4e3))'
        end if

      else
c***   New file (or fresh start with stim >= 0)
c***   Line 1: version, nFb, nq
        write(lun, '(a,3i8)') '#   ', lv, nfb_tmp, nq_tmp

c***   Line 2: run label
        write(lun, '(a,a)') '#  ', trim(label)

c***   Build column header format string
        if (nfb_tmp .gt. 0) then
          write(hlp_frm, '(a,i3,a)')
     &      '(a1,6x,a4,4x,',
     &      nq_tmp * nfb_tmp,
     &      '(4x,a8,i2.2))'
        end if

c*** Write column header line
      if (nfb_tmp .gt. 0) then
        write(lun, hlp_frm) '#', 'time',
     &    ((hlp_hq(j), ifb, ifb=1,nfb_tmp), j=1,nq_tmp)
      end if

c***   Build data format string
        if (nfb_tmp .gt. 0) then
          write(hlp_frm, '(a,i3,a)') '(1p,e16.7,',
     &      nq_tmp * nfb_tmp,
     &      '(2x,e12.4e3))'
        end if
      end if

      return

   90 write(*,*) 'b2fbtrci: cannot open file ', trim(f)
      filedata(ltrc) = .false.
      return

   92 write(*,*) 'b2fbtrci: wrong format of the file ', trim(f)
      filedata(ltrc) = .false.
      close(lun)
      return
      end subroutine b2fbtrci
c=======================================================================
c
      subroutine b2fbtrcc(iline)
      implicit none
      integer, intent(out) :: iline
      real (kind=R8) :: tdummy
      character hlp_c
      logical ex
c-----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c-----------------------------------------------------------------------
      if (.not. filedata(ntrc_fbc)) return
      lun = ntrc(ltrc)
      f = fn_ntrc(ltrc)
      call find_file(f, ex)
      open(lun, file=f, err=190)
      rewind(lun)
      iline = 0
      read(lun, '(a1)', err=155, end=155) hlp_c
      do while (hlp_c .eq. '#')
        iline = iline + 1
        read(lun, '(a1)', err=155, end=155) hlp_c
      end do
      backspace lun
      read(lun, *, err=155, end=155) tdummy
      do while (tdummy .le. tim)
        iline = iline + 1
        read(lun, *, err=155, end=155) tdummy
      end do
  155 continue
      close(lun)
      return
  190 write(*,*) 'b2fbtrcc: cannot open file ', trim(f)
      filedata(ltrc) = .false.
      iline = -1
      return
      end subroutine b2fbtrcc
c=======================================================================
c
      subroutine b2fbtrcs
c-----------------------------------------------------------------------
c*** Save the values (not used for this set)
c-----------------------------------------------------------------------
      return
      end subroutine b2fbtrcs
c=======================================================================
c
      subroutine b2fbtrcf
c-----------------------------------------------------------------------
c*** Finalization step
c-----------------------------------------------------------------------
      logical file_open

      lun = ntrc(ltrc)
      inquire(unit=lun, opened=file_open)
      if (file_open) close(lun)
      return
      end subroutine b2fbtrcf
c=======================================================================
      end module b2mod_fbtrc
