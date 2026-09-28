      subroutine b2agbs_st (nx, ny, nncutmax,
     1 leftix, leftiy, rightix, rightiy, region, crx, cry,
     1 ffbz, bb, wbbl, gs, qz, qc, qac, hx, hy, vol,
     1 pbs, cell_width, cell_height, cflags)
      use b2mod_types
      use b2mod_b2cmfs
      use b2mod_constants
      use carre_constants
      use tradui_constants
      use b2mod_connectivity
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nncutmax
      integer leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny)
      integer rightix(-1:nx,-1:ny), rightiy(-1:nx,-1:ny)
      integer region(-1:nx,-1:ny,0:2)
      integer cflags(-1:nx,-1:ny,GRID_N_CELLFLAGS)
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     *  ffbz(-1:nx,-1:ny,0:3),
     *  bb(-1:nx,-1:ny,0:3), wbbl(-1:nx,-1:ny,0:3),
     *  gs(-1:nx,-1:ny,0:2), qz(-1:nx,-1:ny,0:1),
     *  qc(-1:nx,-1:ny,0:1), qac(-1:nx,-1:ny,0:1), vol(-1:nx,-1:ny)
*   ..input/output arguments (defined on entry, changed upon exit)
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1),
     *  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
*   ..local variables
      integer i, ix, iy, ixl, ixr, iyl, iyr, jxi, jxa, jsep
      integer ixpref, ixsel, ncel, l, lpref, ixw, iyw
      real (kind=R8) :: dist, dist_best
      logical :: okw
      integer ntubes(-1:ny), ixref(nncutmax,-1:ny), iseg(-1:nx,-1:ny),
     *  ixleft(nncutmax,-1:ny),
     *  isnum(nncutmax,-1:ny)
      real (kind=R8) :: hz(-1:nx,-1:ny),
     *                  hzl(-1:nx,-1:ny),    ! hz at the left cell face
     *                  pol_flux, width, wrk
      real (kind=R8) :: cutoff_width
      data width /1.0_R8/
      data cutoff_width /0.2_R8/
      save width, cutoff_width
*     ------------------------------------------------------------------
*     B2AGBS computes a correction to the pbsx array to ensure the
*     constancy of magnetic flux along flux tubes.
*     ------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agbs_st')
      call ipgetr ('b2agmt_1d_width', width)

*   ..preparatory quantities
      hz(:,:) = vol(:,:)/gs(:,:,2)
      if (isymm.eq.0) then
        hzl(:,:) = width
      else if (isymm.eq.1 .or. isymm.eq.2) then
        hzl(:,:) = pi * ( crx(:,:,0) + crx(:,:,2) )
      else if (isymm.eq.3 .or. isymm.eq.4) then
        hzl(:,:) = pi * ( cry(:,:,0) + cry(:,:,2) )
      end if

*   ..toroidal magnetic field
      if (isymm.eq.1 .or. isymm.eq.3) then
        if (maxval(abs(ffbz(:,:,:))).gt.0.0_R8) then
          bb(:,:,2) =
     .     (ffbz(:,:,0) + ffbz(:,:,1) + ffbz(:,:,2) + ffbz(:,:,3)) /
     .     (4.0_R8 * hz(:,:))
          wbbl(:,:,2) = bb(:,:,2)
          where (hzl(:,:).ne.0.0_R8)
            wbbl(:,:,2) = (ffbz(:,:,0) + ffbz(:,:,2)) /
     .       (2.0_R8 * hzl(:,:))
          end where
        end if
      end if
      wbbl(:,:,1) = 0.0_R8

*   ..identify flux tubes
      iseg = 0
      isnum = 0
      ntubes = 0
      ixleft = -2
      do iy = -1, ny
        do ix = -1, nx
          if (region(ix,iy,0).eq.0) cycle
          if (iseg(ix,iy).ne.0) cycle
          if (rightiy(ix,iy).ne.iy) cycle
          if (leftiy(ix,iy).ne.iy) cycle
          if (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          ntubes(iy) = ntubes(iy) + 1
          call xertst (ntubes(iy).le.nncutmax,
     .     'Too many flux tube segments. Increase DEF_NCUT!')
          iseg(ix,iy) = ntubes(iy)
          isnum(iseg(ix,iy),iy) = 1
          ixl = ix
          ixleft(iseg(ix,iy),iy) = ix
          ixr = ix
   20     continue
          if (leftix(ixl,iy).ne.-2.and.leftiy(ixl,iy).eq.iy) then
            ixl = leftix(ixl,iy)
            iyl = leftiy(ixl,iy)
            if (iseg(ixl,iy).ne.0.and.iyl.eq.iy) then
              call xertst(iseg(ixl,iy).eq.iseg(ix,iy),
     .         'Mismatched segments (left)!')
              goto 30
            else if (iyl.eq.iy) then
              iseg(ixl,iy) = iseg(ix,iy)
              isnum(iseg(ix,iy),iy) = isnum(iseg(ix,iy),iy) + 1
              if (leftix(ixl,iy).eq.-2) ixleft(iseg(ix,iy),iy) = ixl
            else if (iyl.ne.iy) then
              goto 30
            end if
            goto 20
          end if
   30     continue
          if (rightix(ixr,iy).ne.nx+1.and.rightiy(ixr,iy).eq.iy) then
            ixr = rightix(ixr,iy)
            iyr = rightiy(ixr,iy)
            if (iseg(ixr,iy).ne.0.and.iyr.eq.iy) then
              call xertst(iseg(ixr,iy).eq.iseg(ix,iy),
     .         'Mismatched segments (right)!')
              goto 40
            else if (iyr.eq.iy) then
              iseg(ixr,iyr) = iseg(ix,iy)
              isnum(iseg(ix,iy),iy) = isnum(iseg(ix,iy),iy) + 1
            else
              goto 40
            end if
            goto 30
          end if
   40     continue
        end do
      end do

*   ..select reference location for each flux tube
*   ..we use the outer midplane whenever possible, then the inner midplane,
*   ..otherwise the middle of the flux tube.
*   ..The reference cell must belong to the flux tube and its western face
*   ..must be an interior face of that same tube. The reference poloidal
*   ..flux is taken from this face; a wall face at the end of the tube has
*   ..a grazing field and would corrupt the poloidal flux of the whole tube.
      ixref = -2
      call get_jsep(nx,ny,jxi,jxa,jsep)
      do iy = -1, ny
        do i = 1, ntubes(iy)
          ! preferred poloidal position of the reference cell
          ixpref = -2
          if (i.eq.iseg(jxa,iy)) ixpref = jxa
          if (ixpref.eq.-2 .and. i.eq.iseg(jxi,iy)) ixpref = jxi
          ncel = isnum(i,iy)
          ! position of this cell along the tube (the tube is a chain of
          ! cells which may not be consecutive in the ix index)
          lpref = -1
          if (ixpref.gt.-2) then
            ix = ixleft(i,iy)
            do l = 1, ncel
              if (ix.eq.ixpref) then
                lpref = l
                exit
              endif
              if ((l.lt.ncel) .and. (rightiy(ix,iy).eq.iy)) then
                ix = rightix(ix,iy)
              else
                exit
              endif
            end do
          endif
          ! look for the best reference among the cells of this tube
          ixsel = -2
          dist_best = 0.0_R8
          ix = ixleft(i,iy)
          do l = 1, ncel
            if (ix.lt.0 .or. ix.gt.nx) exit
            ! the western face must be an interior face of this tube
            ixw = leftix(ix,iy)
            iyw = leftiy(ix,iy)
            okw = (ixw.gt.-2) .and. (iyw.ge.-1) .and. (iyw.le.ny) .and.
     &            (iyw.eq.iy) .and. (iseg(ixw,iy).eq.i) .and.
     &            (.not.isGhostCell(cflags(ixw,iy,CELLFLAG_TYPE))) .and.
     &            (.not.isBoundaryCell(cflags(ixw,iy,CELLFLAG_TYPE)))
            if (okw) then
              if (lpref.gt.0) then
                dist = abs(real(l-lpref))
              else
                dist = abs(real(2*l-ncel-1))
              endif
              if (ixsel.eq.-2 .or. dist.lt.dist_best) then
                dist_best = dist
                ixsel = ix
              endif
            endif
            if (l.lt.ncel .and. rightiy(ix,iy).eq.iy) then
              ix = rightix(ix,iy)
            else
              exit
            endif
          end do
          ixref(i,iy) = ixsel
        end do
      end do

*   ..correct pbsx, hx, and hy to enforce:
*   ..B_{pol}*hy1*hz constant (on western cell faces and at cell centers)
*   ..vol = hy1*hx*hz
      do iy = -1, ny
        do ix = -1, nx
          if (region(ix,iy,0).eq.0) cycle
          if (iseg(ix,iy).eq.0) cycle
          if (ixref(iseg(ix,iy),iy).eq.-2) cycle

          ! poloidal flux through the western face of the reference cell
          ! (qac is the field angle, as in b2xgbs, not the mesh angle qc)
          if (leftix(ixref(iseg(ix,iy),iy),iy).ne.-2) then
            pol_flux = wbbl(ixref(iseg(ix,iy),iy),iy,0) *
     &                   gs(ixref(iseg(ix,iy),iy),iy,0) *
     &                   qac(ixref(iseg(ix,iy),iy),iy,0)
          else
            pol_flux = bb(ixref(iseg(ix,iy),iy),iy,0) *
     &                 hy(ixref(iseg(ix,iy),iy),iy)
          end if
          wrk = wbbl(ix,iy,0)*gs(ix,iy,0)*qac(ix,iy,0)
          if ((abs(wrk).lt.abs(pol_flux)) .and.
     &        isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) ) cycle
          ! do not force a uniform flux where the poloidal flux of the
          ! cell (its psi width) differs from the reference value, as in
          ! the unstructured b2agbs
          if (abs(abs(wrk)-abs(pol_flux))
     &        .gt. cutoff_width*abs(pol_flux)) cycle

          wrk = pol_flux /
     .         (hz(ix,iy)*qz(ix,iy,1)*bb(ix,iy,0))
          if (wrk .lt. 0.0_R8) then
            write(*,'(a,i2,a,i2,a)') "For cell (ix,iy) (", ix,",", iy,
     &        ") redefined hy < 0!"
            cycle
          endif
          if (leftix(ix,iy).ne.-2) then ! do not adjust pbs if leftix=-2
                                        ! since there is no western face
            wbbl(ix,iy,0) = pol_flux / (gs(ix,iy,0) * qac(ix,iy,0))
            wbbl(ix,iy,3) =
     &       sqrt(wbbl(ix,iy,0)**2+wbbl(ix,iy,1)**2+wbbl(ix,iy,2)**2)
            pbs(ix,iy,0) = pol_flux / wbbl(ix,iy,3)
          else
            wbbl(ix,iy,0) = pol_flux /
     &       (gs(rightix(ix,iy),rightiy(ix,iy),0) *
     &        qac(rightix(ix,iy),rightiy(ix,iy),0))
            wbbl(ix,iy,3) =
     &       sqrt(wbbl(ix,iy,0)**2+wbbl(ix,iy,1)**2+wbbl(ix,iy,2)**2)
          end if
          hy(ix,iy) = pol_flux /
     .               (hz(ix,iy)*qz(ix,iy,1)*bb(ix,iy,0))
          hx(ix,iy) = vol(ix,iy)/
     .               (hz(ix,iy)*hy(ix,iy)*qz(ix,iy,1))
        end do
      end do

*   ..copying over modified data
      cell_width = hx
      cell_height = hy

      bb(:,:,3)=sqrt(bb(:,:,0)**2+bb(:,:,1)**2+bb(:,:,2)**2)

*   ..return
      call subend ()
      return
      end subroutine b2agbs_st

