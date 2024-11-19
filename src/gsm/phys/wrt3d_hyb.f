      subroutine wrt3d_hyb(dt6dt,kdt,global_lats_r,lonsperlar,nblck,
     &                     idate)
      use resol_def
      use layout1
      use namelist_physics_def
      use machine
      implicit none
!!
      integer, parameter :: numv = 7
      integer, parameter :: l0 = 90
      integer, parameter :: ioproc=0
!!
      integer, intent(in) :: kdt,nblck
      real (kind=kind_rad),dimension(ngptc,levs,numv,nblck,lats_node_r),
     &                     intent(in) :: dt6dt
      integer, intent(in), dimension(latr) :: global_lats_r, lonsperlar
      integer, dimension(4), intent(in) :: idate
! Local variables
      integer :: nv,k,lan,lon,i,iblk,lons_lat,njeff,il,lat
      real (kind=kind_io4) :: wrkga(lonr,latr)
      real (kind=kind_io8) :: rtime
      real (kind=kind_io8) :: glolal(lonr,lats_node_r)
      real (kind=kind_io8) :: buffo(lonr,lats_node_r)
      integer:: kmsk(lonr*latr)
      integer, dimension(nodes) :: lats_nodes_r
!..........................................................
!     temperature tendencies
!
      kmsk = 0
      do nv=1,numv
        if (nv .ge. 2 .and. nv .le. 4) cycle
        do k=l0,levs
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
                glolal(il,lan)=dt6dt(i,k,nv,iblk,lan) ! *rtime
                il=il+1
              enddo
            enddo
          enddo
          call uninterpred(1,kmsk,buffo,glolal,global_lats_r,lonsperlar)
          call unsplit2d_phys(ioproc,wrkga,buffo,global_lats_r)
          if(me.eq.ioproc) call ncwrite(nv, k, kdt, wrkga, idate)
        enddo
      enddo

      return
      end subroutine wrt3d_hyb

!!

      subroutine ncwrite(nv, k, kdt, workga, idate)
      use netcdf
      use layout1
      use gg_def
      use machine
      use resol_def
      implicit none
!
      integer, intent(in) :: nv, k, kdt
      integer, intent(in), dimension(4) :: idate
      real(kind=kind_io4),dimension(lonr,latr), intent(in) :: workga
! Local variables
      real(kind=kind_io4),dimension(lonr) :: lons
      real :: pi

      integer :: id, x_dimid, x_varid, y_dimid, y_varid, t_dimid,
     &           t_varid, varid, t_len, i, adj, z_dimid, z_varid
      integer, parameter :: numv = 7
      integer, parameter :: num_variables = 7
      integer, parameter :: l0 = 90
      integer, dimension(4) :: dimids
      integer, dimension(levs) :: nlevs

      character(len=*), dimension(num_variables), parameter :: vars  =
     &           (/"qno","wtot","euv","uv","cp","rho","z"/)
      character(len=*), dimension(num_variables), parameter :: lvars =
     &           (/"NO Cooling","Total Heating","EUV Heating",
     &             "UV Heating","Specific Heat","Total Density",
     &             "Geometric Height"/)
      character(len=*), dimension(num_variables), parameter :: units =
     &           (/"K/s","K/s","K/s","K/s","J/kg/K","kg/m^3","m"/)
      character(len=8), parameter :: fname = "dt6dt.nc"

      character(len=31) :: time_units

      logical :: exists
!
      inquire( file = fname, exist=exists )

      if (.not. exists) then
        write(time_units, "(A12,I0.4,A1,I0.2,A1,I0.2,A1,I0.2,A6)"),
     &   "hours since ", idate(4), "-", idate(2), "-", idate(3),
     &   " ", idate(1),":00:00"
        call check(nf90_create(fname, nf90_clobber, id))

        call check(nf90_def_dim(id,"lon", lonr, x_dimid))
        call check(nf90_def_dim(id,"lat", latr, y_dimid))
        call check(nf90_def_dim(id,"levs",levs-l0+1, z_dimid))
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
        call check(nf90_put_att(id,t_varid,"units",time_units))

        dimids = (/ x_dimid, y_dimid, z_dimid, t_dimid /)

        do i=1,numv
          if (i .ge. 2 .and. i .le. 4) cycle
          call check(nf90_def_var(id,vars(i),NF90_REAL,dimids,varid))
          call check(nf90_put_att(id,varid,"units",units(i)))
          call check(nf90_put_att(id,varid,"long_name", lvars(i)))
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
        call check(nf90_put_var(id, z_varid, nlevs(l0:levs)))

        call check(nf90_close(id))
      end if

      call check(nf90_open(fname, NF90_WRITE, id))

      call check(nf90_inq_varid(id, vars(nv), varid))
      call check(nf90_put_var(id,varid,workga,start=(/1,1,k-l0+1,kdt/),
     &                              count=(/lonr,latr,1,1/)))

      if (nv .eq. 1) then
        call check(nf90_inq_varid(id, "time", varid))
        call check(nf90_put_var(id,varid,real(kdt),start=(/kdt/)))
      end if

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
