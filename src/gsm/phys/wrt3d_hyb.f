      subroutine wrt3d_hyb(ioproc,kdt,global_lats_r,lonsperlar)
      use resol_def
      use layout1
      use namelist_physics_def
      use d3d_def
      use machine
      implicit none
!! do you need zhour and fhour ?
      integer, intent(in) :: ioproc, kdt
      integer, intent(in), dimension(latr) :: global_lats_r, lonsperlar
! Local variables
      integer :: nv,k,lan,lon,i,iblk,lons_lat,njeff,il,lat
      real (kind=kind_io4) :: wrkga(lonr,latr)
      real (kind=kind_io8) :: rtime
      real (kind=kind_io8) :: glolal(lonr,lats_node_r)
      real (kind=kind_io8) :: buffo(lonr,lats_node_r)
!!
      integer:: kmsk(lonr*latr)
      integer, dimension(nodes) :: lats_nodes_r
!
!
!      if(fhour.gt.zhour) then
!        rtime=1./(3600.*(fhour-zhour))
!      else
!        rtime=0.
!      endif
!..........................................................
!     temperature tendencies
!
      kmsk = 0
      do nv=1,1 ! do not output other variables
        do k=1,levs
          glolal=0.
          do lan=1,lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+lan)
            lons_lat = lonsperlar(lat)
            iblk=0
            il=1
            do lon=1,lons_lat,ngptc
              njeff=min(ngptc,lons_lat-lon+1)
              iblk=iblk+1
              do i=1,njeff
                glolal(il,lan)=dt3dt(i,k,nv,iblk,lan) ! *rtime
                il=il+1
              enddo
            enddo
          enddo
          call uninterpred(1,kmsk,buffo,glolal,global_lats_r,lonsperlar)
          call unsplit2d_phys(ioproc,wrkga,buffo,global_lats_r)
          if(me.eq.ioproc) call ncwrite(nv, k, kdt, wrkga)
        enddo
      enddo

      return
      end subroutine wrt3d_hyb

!!

      subroutine ncwrite(nv, k, kdt, workga)
      use netcdf
      use layout1
      use gg_def
      use machine
      use resol_def
      implicit none
!
      integer, intent(in) :: nv, k, kdt
      real(kind=kind_io4),dimension(lonr,latr), intent(in) :: workga
! Local variables
      real(kind=kind_io4),dimension(lonr) :: lons
      real :: pi

      integer :: id, x_dimid, x_varid, y_dimid, y_varid, t_dimid,
     &           t_varid, varid, t_len, i, adj, z_dimid, z_varid
      integer, parameter :: numv = 1
      integer, parameter :: num_variables = 6
      integer, dimension(4) :: dimids
      integer, dimension(levs) :: nlevs

      character(len=*), dimension(num_variables), parameter :: vars  =
     &           (/"qno","un1","euv","uv","un2","un3"/)
      character(len=*), dimension(num_variables), parameter :: lvars =
     &           (/"NO Cooling","Unknown","EUV Heating","UV Heating",
     &             "Unknown","Unknown"/)
      character(len=8), parameter :: fname = "dt3dt.nc"

      logical :: exists
!
      inquire( file = fname, exist=exists )

      if (.not. exists) then
        call check(nf90_create(fname, nf90_clobber, id))

        call check(nf90_def_dim(id,"lon", lonr, x_dimid))
        call check(nf90_def_dim(id,"lat", latr, y_dimid))
        call check(nf90_def_dim(id,"levs",levs, z_dimid))
        call check(nf90_def_dim(id,"time",NF90_UNLIMITED,t_dimid))

        call check(nf90_def_var(id,"lon",NF90_REAL,x_dimid,x_varid))
        call check(nf90_put_att(id,x_varid,"axis","X"))
        call check(nf90_put_att(id,x_varid,"long_name","longitude"))
        call check(nf90_put_att(id,x_varid,"units","degrees_east"))
        call check(nf90_def_var(id,"lat",NF90_REAL,y_dimid,y_varid))
        call check(nf90_put_att(id,y_varid,"axis","Y"))
        call check(nf90_put_att(id,y_varid,"long_name","latitude"))
        call check(nf90_put_att(id,y_varid,"units","degrees_north"))
        call check(nf90_def_var(id,"levs",NF90_INT,z_dimid,z_varid))
        call check(nf90_put_att(id,y_varid,"axis","Z"))
        call check(nf90_put_att(id,y_varid,"long_name","model levels"))
        call check(nf90_put_att(id,y_varid,"units","1"))
        call check(nf90_def_var(id,"time",NF90_REAL,t_dimid,t_varid))
        call check(nf90_put_att(id,t_varid,"axis","T"))
        call check(nf90_put_att(id,t_varid,"long_name","time"))
        call check(nf90_put_att(id,t_varid,"units",
     &                          "hours since 1970-01-01"))

        dimids = (/ x_dimid, y_dimid, z_dimid, t_dimid /)

        do i=1,numv
          call check(nf90_def_var(id,vars(i),NF90_REAL,dimids,varid))
          call check(nf90_put_att(id,varid,"units","1"))
          call check(nf90_put_att(id,varid,"long_name", lvars(nv)))
        end do

        call check(nf90_enddef(id))

        do i=1,lonr
          lons(i) = (i-1.) / lonr * 360
        end do

        do i=1,levs
          nlevs(i) = i
        end do

        pi = acos(-1.0)

        call check(nf90_put_var(id, x_varid, lons))
        call check(nf90_put_var(id, y_varid, asin(sinlat_r)*180/pi))
        call check(nf90_put_var(id, z_varid, nlevs))

        call check(nf90_close(id))
      end if

      call check(nf90_open(fname, NF90_WRITE, id))

      call check(nf90_inq_varid(id, vars(nv), varid))
      call check(nf90_put_var(id,varid,workga,start=(/1,1,k,kdt/),
     &                              count=(/lonr,latr,1,1/)))

      call check(nf90_close(id))

      end subroutine ncwrite

!!

      subroutine check(istatus)

      use netcdf

      implicit none

      integer, intent(in) :: istatus

      if (istatus /= nf90_noerr) then
        write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
      end if

      end subroutine check
