      subroutine wrt3d_hyb(ioproc,kdt,global_lats_r,lonsperlar)
      use resol_def
      use layout1
      use namelist_physics_def
      use d3d_def
      implicit none
!! do you need zhour and fhour ?
      integer, intent(in) :: ioproc, kdt
      integer, intent(in), dimension(latr) :: global_lats_r, lonsperlar
! Local variables
      integer :: nv,k,lan,lon,i,iblk,lons_lat,njeff
      real (kind=kind_io4) :: wrkga(lonr*latr)
      real (kind=kind_io8) :: rtime
      real (kind=kind_io8) :: glolal(lonr,lats_node_r)
      real (kind=kind_io8) :: buffo(lonr,lats_node_r)
!!
      integer, parameter :: kmsk0(lonr*latr) = 0
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
      do nv=1,6
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
          call uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!          call unsplit2d_phys(ioproc,wrkga,buffo,global_lats_r) ! think I don't need this?
          if(me.eq.ioproc) call ncwrite(nv, buffo)
        enddo
      enddo

      return
      end subroutine wrt3d_hyb

!!

      subroutine ncwrite(nv, buffo)
      use netcdf
      use layout1
      use gg_def
      implicit none
!
      integer, intent(in) :: nv
      real(kind=kind_io8),dimension(lonr,lats_node_r,1),intent(in) ::   
     &                                                            buffo
! Local variables
      real(kind=kind_io4),dimension(lonr) :: lons

      integer :: id, x_dimid, x_varid, y_dimid, y_varid, t_dimid,       
     &           t_varid, varid, t_len, i, adj
      integer, parameter :: num_variables = 6
      integer, dimension(3) :: dimids

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
        call check(nf90_create(fname, nf90_clobber, ncid))

        call check(nf90_def_dim(id,"lon", lonr,          x_dimid))
        call check(nf90_def_dim(id,"lat", lats_node_r,   y_dimid))
        call check(nf90_def_dim(id,"time",NF90_UNLIMITED,t_dimid))

        call check(nf90_def_var(id,"lon",NF90_REAL,x_dimid,x_varid))
        call check(nf90_put_att(id,x_varid,"axis","X"))
        call check(nf90_put_att(id,x_varid,"long_name","longitude"))
        call check(nf90_put_att(id,x_varid,"units","degrees_east"))
        call check(nf90_def_var(id,"lat",NF90_REAL,y_dimid,y_varid))
        call check(nf90_put_att(id,y_varid,"axis","Y"))
        call check(nf90_put_att(id,y_varid,"long_name","latitude"))
        call check(nf90_put_att(id,y_varid,"units","degrees_north"))
        call check(nf90_def_var(id,"time",NF90_DOUBLE,t_dimid,t_varid))
        call check(nf90_put_att(id,t_varid,"axis","T"))
        call check(nf90_put_att(id,t_varid,"long_name","time"))
        call check(nf90_put_att(id,t_varid,"units",                     
       &                        "days since 1970-01-01"))

        dimids = (/ x_dimid, y_dimid, t_dimid /)

        do i=1,nv
          call check(nf90_def_var(id,vars(i),NF90_DOUBLE,dimids,varid)
          call check(nf90_put_att(id,varid,"units","1"))
          call check(nf90_put_att(id,varid,"long_name", lvars(nv)))
        end do

        call check(nf90_enddef(ncid))

        do i=1,lonr
          lon(i) = (i-1) / lonr * 360
        end do

        call check(nf90_put_var(ncid, x_varid, lon))
        call check(nf90_put_var(ncid, y_varid, asin(sinlat_r)))

        call check(nf90_close(id))
      end if

      call check(nf90_inq_dimid(id, "time", t_dimid))
      call check(nf90_inquire_dimension(id, t_dimid, len=t_len))

      adj = 0
      if ( nv = 1 ) then
        adj = 1
        call check(nf90_inq_varid(id, "time", t_varid))
        call check(nf90_put_var(id, t_varid, (/0./), start=t_len+1)
      end if

      call check(nf90_inq_varid(id, vars(nv), varid))
      call check(nf90_put_var(id, varid, buffo, start=t_len+adj)

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
