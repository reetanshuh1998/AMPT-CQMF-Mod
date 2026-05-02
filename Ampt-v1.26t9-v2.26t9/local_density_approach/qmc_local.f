c=============================================================
c  LOCAL DENSITY CQMF SUBROUTINES (iqmc=2 mode)
c  Nuclear saturation density rho0 = 0.16 fm^-3
c  CSV table density column is in units of rho/rho0
c=============================================================

      subroutine pos_to_cell(x,y,z,ix,iy,iz,inside)
      implicit double precision (a-h,o-z)
      integer ix,iy,iz,inside
      common /ilist3/ size1, size2, size3, v1, v2, v3, size
      save

      xmin = -5d0*size1
      xmax =  5d0*size1
      ymin = -5d0*size2
      ymax =  5d0*size2
      zmin = -5d0*size3
      zmax =  5d0*size3

      inside = 1
      if (x.lt.xmin .or. x.gt.xmax .or.
     &    y.lt.ymin .or. y.gt.ymax .or.
     &    z.lt.zmin .or. z.gt.zmax) then
         inside = 0
         ix=1
         iy=1
         iz=1
         return
      endif

      dx = (xmax-xmin)/10d0
      dy = (ymax-ymin)/10d0
      dz = (zmax-zmin)/10d0

      ix = int((x - xmin)/dx) + 1
      iy = int((y - ymin)/dy) + 1
      iz = int((z - zmin)/dz) + 1

      if(ix.lt.1) ix=1
      if(ix.gt.10) ix=10
      if(iy.lt.1) iy=1
      if(iy.gt.10) iy=10
      if(iz.lt.1) iz=1
      if(iz.gt.10) iz=10

      return
      end

c-------------------------------------------------------------
      subroutine build_rhob_grid(tcur)
      implicit double precision (a-h,o-z)
      parameter (MAXPTN=400001)
      common /para1/ mul
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN),py(MAXPTN),pz(MAXPTN),e(MAXPTN),
     &     xmass(MAXPTN),ityp(MAXPTN)
      common /ilist3/ size1, size2, size3, v1, v2, v3, size
      common /qmcgrid/ rhob(10,10,10),
     &     mu_c(10,10,10),md_c(10,10,10),ms_c(10,10,10),
     &     vu_c(10,10,10),vd_c(10,10,10),vs_c(10,10,10)
      save /qmcgrid/
      double precision Bcell(10,10,10)
      integer ix,iy,iz,inside
      save

      do iz=1,10
      do iy=1,10
      do ix=1,10
         Bcell(ix,iy,iz)=0d0
      enddo
      enddo
      enddo

      do i=1,mul
         if (tcur .lt. ft(i)) goto 10
         call pos_to_cell(gx(i),gy(i),gz(i),ix,iy,iz,inside)
         if (inside .eq. 0) goto 10
         if (abs(ityp(i)).ge.1 .and. abs(ityp(i)).le.3) then
            if (ityp(i) .gt. 0) then
               Bcell(ix,iy,iz) = Bcell(ix,iy,iz) + 1d0/3d0
            else
               Bcell(ix,iy,iz) = Bcell(ix,iy,iz) - 1d0/3d0
            endif
         endif
 10      continue
      enddo

      dx = (10d0*size1)/10d0
      dy = (10d0*size2)/10d0
      dz = (10d0*size3)/10d0
      Vcell = dx*dy*dz
      if (Vcell .lt. 1d-20) Vcell = 1d-20

      do iz=1,10
      do iy=1,10
      do ix=1,10
         rhob(ix,iy,iz) = Bcell(ix,iy,iz)/Vcell
      enddo
      enddo
      enddo

      return
      end

c-------------------------------------------------------------
      subroutine smooth_rhob()
      implicit double precision (a-h,o-z)
      common /qmcgrid/ rhob(10,10,10),
     &     mu_c(10,10,10),md_c(10,10,10),ms_c(10,10,10),
     &     vu_c(10,10,10),vd_c(10,10,10),vs_c(10,10,10)
      save /qmcgrid/
      double precision tmp(10,10,10)
      integer ix,iy,iz,jx,jy,jz,dx,dy,dz,n
      save

      do iz=1,10
      do iy=1,10
      do ix=1,10
         sum1=0d0
         n=0
         do dz=-1,1
         do dy=-1,1
         do dx=-1,1
            jx=ix+dx
            jy=iy+dy
            jz=iz+dz
            if(jx.ge.1.and.jx.le.10.and.
     &         jy.ge.1.and.jy.le.10.and.
     &         jz.ge.1.and.jz.le.10) then
               sum1 = sum1 + rhob(jx,jy,jz)
               n = n + 1
            endif
         enddo
         enddo
         enddo
         if(n.gt.0) then
            tmp(ix,iy,iz) = sum1/dble(n)
         else
            tmp(ix,iy,iz) = 0d0
         endif
      enddo
      enddo
      enddo

      do iz=1,10
      do iy=1,10
      do ix=1,10
         rhob(ix,iy,iz) = tmp(ix,iy,iz)
      enddo
      enddo
      enddo

      return
      end

c-------------------------------------------------------------
      subroutine interp_qmc_table(rho_ratio,
     &     omu,omd,oms,ovu,ovd,ovs)
      implicit double precision (a-h,o-z)
      common /qmctbl/ den_t(200),mu_t(200),md_t(200),
     &     ms_t(200),vu_t(200),vd_t(200),vs_t(200),ntbl
      save /qmctbl/
      double precision rho_ratio,omu,omd,oms,ovu,ovd,ovs,frac
      save

      if(ntbl.le.0) then
         omu=316.8d0
         omd=316.8d0
         oms=497.4d0
         ovu=0d0
         ovd=0d0
         ovs=0d0
         return
      endif

      rr = rho_ratio
      if(rr.lt.den_t(1)) rr=den_t(1)
      if(rr.gt.den_t(ntbl)) rr=den_t(ntbl)

      omu=mu_t(1)
      omd=md_t(1)
      oms=ms_t(1)
      ovu=vu_t(1)
      ovd=vd_t(1)
      ovs=vs_t(1)

      do i=1,ntbl-1
         if(rr.ge.den_t(i) .and. rr.le.den_t(i+1)) then
            frac=(rr-den_t(i))/(den_t(i+1)-den_t(i))
            omu=mu_t(i)+frac*(mu_t(i+1)-mu_t(i))
            omd=md_t(i)+frac*(md_t(i+1)-md_t(i))
            oms=ms_t(i)+frac*(ms_t(i+1)-ms_t(i))
            ovu=vu_t(i)+frac*(vu_t(i+1)-vu_t(i))
            ovd=vd_t(i)+frac*(vd_t(i+1)-vd_t(i))
            ovs=vs_t(i)+frac*(vs_t(i+1)-vs_t(i))
         endif
      enddo

      return
      end

c-------------------------------------------------------------
      subroutine load_qmc_table()
      implicit double precision (a-h,o-z)
      common /qmctbl/ den_t(200),mu_t(200),md_t(200),
     &     ms_t(200),vu_t(200),vd_t(200),vs_t(200),ntbl
      save /qmctbl/
      character*200 line
      double precision cd,cu,cdq,cs,cb,cvu,cvd,cvs
      save

      ntbl=0
      open(89, file='model_data.csv', status='old', err=999)
      read(89, '(A)') line
 100  continue
      read(89, *, end=200) cd,cu,cdq,cs,cb,cvu,cvd,cvs
      ntbl=ntbl+1
      den_t(ntbl)=cd
      mu_t(ntbl)=cu
      md_t(ntbl)=cdq
      ms_t(ntbl)=cs
      vu_t(ntbl)=cvu
      vd_t(ntbl)=cvd
      vs_t(ntbl)=cvs
      goto 100
 200  continue
      close(89)
      write(6,*) 'CQMF table loaded, npts=',ntbl
      return
 999  write(6,*) 'WARNING: model_data.csv missing for local mode'
      return
      end

c-------------------------------------------------------------
      subroutine update_cell_fields()
      implicit double precision (a-h,o-z)
      parameter (rho0 = 0.16d0)
      common /qmcgrid/ rhob(10,10,10),
     &     mu_c(10,10,10),md_c(10,10,10),ms_c(10,10,10),
     &     vu_c(10,10,10),vd_c(10,10,10),vs_c(10,10,10)
      save /qmcgrid/
      double precision omu,omd,oms,ovu,ovd,ovs,rr
      save

      do iz=1,10
      do iy=1,10
      do ix=1,10
         rr = dabs(rhob(ix,iy,iz)) / rho0
         call interp_qmc_table(rr,omu,omd,oms,ovu,ovd,ovs)
         mu_c(ix,iy,iz) = omu/1000d0
         md_c(ix,iy,iz) = omd/1000d0
         ms_c(ix,iy,iz) = oms/1000d0
         vu_c(ix,iy,iz) = ovu/1000d0
         vd_c(ix,iy,iz) = ovd/1000d0
         vs_c(ix,iy,iz) = ovs/1000d0
      enddo
      enddo
      enddo

      return
      end

c-------------------------------------------------------------
      subroutine assign_parton_local_medium(tcur)
      implicit double precision (a-h,o-z)
      parameter (MAXPTN=400001)
      common /para1/ mul
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN),py(MAXPTN),pz(MAXPTN),e(MAXPTN),
     &     xmass(MAXPTN),ityp(MAXPTN)
      common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
      common /qmcgrid/ rhob(10,10,10),
     &     mu_c(10,10,10),md_c(10,10,10),ms_c(10,10,10),
     &     vu_c(10,10,10),vd_c(10,10,10),vs_c(10,10,10)
      save /qmcgrid/
      common /qmcv0p/ v0p(MAXPTN)
      save /qmcv0p/
      integer ix,iy,iz,inside,ifl
      double precision xm_loc,v_loc
      save

      do i=1,mul
         if (tcur .lt. ft(i)) goto 20

         call pos_to_cell(gx(i),gy(i),gz(i),ix,iy,iz,inside)
         if (inside .eq. 0) goto 20

         ifl = ityp(i)
         if (abs(ifl).eq.21 .or. abs(ifl).eq.9) then
            xm_loc = 0d0
            v_loc  = 0d0
         elseif (abs(ifl).eq.3) then
            xm_loc = ms_c(ix,iy,iz)
            v_loc  = vs_c(ix,iy,iz)
         elseif (abs(ifl).eq.1) then
            xm_loc = md_c(ix,iy,iz)
            v_loc  = vd_c(ix,iy,iz)
         else
            xm_loc = mu_c(ix,iy,iz)
            v_loc  = vu_c(ix,iy,iz)
         endif

         if (ifl .lt. 0) v_loc = -v_loc

         xmass(i) = xm_loc
         e(i) = dsqrt(px(i)**2+py(i)**2+pz(i)**2
     &              +xmass(i)**2)
         v0p(i) = v_loc

         if(e(i).gt.1d-10) then
            vx(i) = px(i)/e(i)
            vy(i) = py(i)/e(i)
            vz(i) = pz(i)/e(i)
         endif

 20      continue
      enddo

      return
      end

c-------------------------------------------------------------
      subroutine qmc_medium_step(tcur)
      implicit double precision (a-h,o-z)
      common /qmcpar/ xmu_q, xmd_q, xms_q, iqmc
      common /qmcupd/ dt_med, t_last_med
      save /qmcupd/
      save

      if (iqmc .ne. 2) return
      if (tcur - t_last_med .lt. dt_med) return

      call build_rhob_grid(tcur)
      call smooth_rhob()
      call update_cell_fields()
      call assign_parton_local_medium(tcur)
      t_last_med = tcur

      return
      end

c-------------------------------------------------------------
c  PHASE 2: Mean-field force subroutines
c-------------------------------------------------------------
      subroutine compute_gradients()
      implicit double precision (a-h,o-z)
      common /ilist3/ size1, size2, size3, v1, v2, v3, size
      common /qmcgrid/ rhob(10,10,10),
     &     mu_c(10,10,10),md_c(10,10,10),ms_c(10,10,10),
     &     vu_c(10,10,10),vd_c(10,10,10),vs_c(10,10,10)
      save /qmcgrid/
      common /qmcgrad/
     &     dvu_x(10,10,10),dvu_y(10,10,10),dvu_z(10,10,10),
     &     dvd_x(10,10,10),dvd_y(10,10,10),dvd_z(10,10,10),
     &     dvs_x(10,10,10),dvs_y(10,10,10),dvs_z(10,10,10),
     &     dmu_x(10,10,10),dmu_y(10,10,10),dmu_z(10,10,10),
     &     dmd_x(10,10,10),dmd_y(10,10,10),dmd_z(10,10,10),
     &     dms_x(10,10,10),dms_y(10,10,10),dms_z(10,10,10)
      save /qmcgrad/

      dx = (10d0*size1)/10d0
      dy = (10d0*size2)/10d0
      dz = (10d0*size3)/10d0

      do iz=1,10
      do iy=1,10
      do ix=1,10
         ixp=min(ix+1,10)
         ixm=max(ix-1,1)
         iyp=min(iy+1,10)
         iym=max(iy-1,1)
         izp=min(iz+1,10)
         izm=max(iz-1,1)
         ddx=dble(ixp-ixm)*dx
         ddy=dble(iyp-iym)*dy
         ddz=dble(izp-izm)*dz
         if(ddx.lt.1d-20) ddx=1d0
         if(ddy.lt.1d-20) ddy=1d0
         if(ddz.lt.1d-20) ddz=1d0

         dvu_x(ix,iy,iz)=(vu_c(ixp,iy,iz)-vu_c(ixm,iy,iz))/ddx
         dvu_y(ix,iy,iz)=(vu_c(ix,iyp,iz)-vu_c(ix,iym,iz))/ddy
         dvu_z(ix,iy,iz)=(vu_c(ix,iy,izp)-vu_c(ix,iy,izm))/ddz
         dvd_x(ix,iy,iz)=(vd_c(ixp,iy,iz)-vd_c(ixm,iy,iz))/ddx
         dvd_y(ix,iy,iz)=(vd_c(ix,iyp,iz)-vd_c(ix,iym,iz))/ddy
         dvd_z(ix,iy,iz)=(vd_c(ix,iy,izp)-vd_c(ix,iy,izm))/ddz
         dvs_x(ix,iy,iz)=(vs_c(ixp,iy,iz)-vs_c(ixm,iy,iz))/ddx
         dvs_y(ix,iy,iz)=(vs_c(ix,iyp,iz)-vs_c(ix,iym,iz))/ddy
         dvs_z(ix,iy,iz)=(vs_c(ix,iy,izp)-vs_c(ix,iy,izm))/ddz

         dmu_x(ix,iy,iz)=(mu_c(ixp,iy,iz)-mu_c(ixm,iy,iz))/ddx
         dmu_y(ix,iy,iz)=(mu_c(ix,iyp,iz)-mu_c(ix,iym,iz))/ddy
         dmu_z(ix,iy,iz)=(mu_c(ix,iy,izp)-mu_c(ix,iy,izm))/ddz
         dmd_x(ix,iy,iz)=(md_c(ixp,iy,iz)-md_c(ixm,iy,iz))/ddx
         dmd_y(ix,iy,iz)=(md_c(ix,iyp,iz)-md_c(ix,iym,iz))/ddy
         dmd_z(ix,iy,iz)=(md_c(ix,iy,izp)-md_c(ix,iy,izm))/ddz
         dms_x(ix,iy,iz)=(ms_c(ixp,iy,iz)-ms_c(ixm,iy,iz))/ddx
         dms_y(ix,iy,iz)=(ms_c(ix,iyp,iz)-ms_c(ix,iym,iz))/ddy
         dms_z(ix,iy,iz)=(ms_c(ix,iy,izp)-ms_c(ix,iy,izm))/ddz
      enddo
      enddo
      enddo

      return
      end

c-------------------------------------------------------------
      subroutine apply_force_kick(ii, dtf)
      implicit double precision (a-h,o-z)
      parameter (MAXPTN=400001)
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN),py(MAXPTN),pz(MAXPTN),e(MAXPTN),
     &     xmass(MAXPTN),ityp(MAXPTN)
      common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
      common /qmcgrad/
     &     dvu_x(10,10,10),dvu_y(10,10,10),dvu_z(10,10,10),
     &     dvd_x(10,10,10),dvd_y(10,10,10),dvd_z(10,10,10),
     &     dvs_x(10,10,10),dvs_y(10,10,10),dvs_z(10,10,10),
     &     dmu_x(10,10,10),dmu_y(10,10,10),dmu_z(10,10,10),
     &     dmd_x(10,10,10),dmd_y(10,10,10),dmd_z(10,10,10),
     &     dms_x(10,10,10),dms_y(10,10,10),dms_z(10,10,10)
      save /qmcgrad/
      integer ix,iy,iz,inside,ifl
      double precision fx,fy,fz,moe,sgn
      save

      ifl = ityp(ii)
      if(abs(ifl).eq.21.or.abs(ifl).eq.9) return

      call pos_to_cell(gx(ii),gy(ii),gz(ii),ix,iy,iz,inside)
      if(inside.eq.0) return

      sgn = 1d0
      if(ifl.lt.0) sgn = -1d0

      if(e(ii).gt.1d-10) then
         moe = xmass(ii)/e(ii)
      else
         moe = 0d0
      endif

      if(abs(ifl).eq.3) then
         fx = -(sgn*dvs_x(ix,iy,iz) + moe*dms_x(ix,iy,iz))
         fy = -(sgn*dvs_y(ix,iy,iz) + moe*dms_y(ix,iy,iz))
         fz = -(sgn*dvs_z(ix,iy,iz) + moe*dms_z(ix,iy,iz))
      elseif(abs(ifl).eq.1) then
         fx = -(sgn*dvd_x(ix,iy,iz) + moe*dmd_x(ix,iy,iz))
         fy = -(sgn*dvd_y(ix,iy,iz) + moe*dmd_y(ix,iy,iz))
         fz = -(sgn*dvd_z(ix,iy,iz) + moe*dmd_z(ix,iy,iz))
      else
         fx = -(sgn*dvu_x(ix,iy,iz) + moe*dmu_x(ix,iy,iz))
         fy = -(sgn*dvu_y(ix,iy,iz) + moe*dmu_y(ix,iy,iz))
         fz = -(sgn*dvu_z(ix,iy,iz) + moe*dmu_z(ix,iy,iz))
      endif

      px(ii) = px(ii) + fx*dtf
      py(ii) = py(ii) + fy*dtf
      pz(ii) = pz(ii) + fz*dtf
      e(ii)  = dsqrt(px(ii)**2+py(ii)**2+pz(ii)**2
     &              +xmass(ii)**2)
      if(e(ii).gt.1d-10) then
         vx(ii) = px(ii)/e(ii)
         vy(ii) = py(ii)/e(ii)
         vz(ii) = pz(ii)/e(ii)
      endif

      return
      end

c-------------------------------------------------------------
      subroutine qmc_evolve_all_partons(t_old, t_new)
      implicit double precision (a-h,o-z)
      parameter (MAXPTN=400001)
      common /para1/ mul
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN),py(MAXPTN),pz(MAXPTN),e(MAXPTN),
     &     xmass(MAXPTN),ityp(MAXPTN)
      common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
      common /qmcpar/ xmu_q, xmd_q, xms_q, iqmc
      common /qmcupd/ dt_med, t_last_med
      save /qmcupd/
      double precision dtf, tcur, dt_force
      save

      if(iqmc.ne.2) return
      if(t_new.le.t_old) return

      dt_force = 0.05d0
      tcur = t_old

 50   continue
      if(tcur.ge.t_new) return
      dtf = min(dt_force, t_new - tcur)

      if(tcur - t_last_med .ge. dt_med) then
         call build_rhob_grid(tcur)
         call smooth_rhob()
         call update_cell_fields()
         call compute_gradients()
         call assign_parton_local_medium(tcur)
         t_last_med = tcur
      endif

      do i=1,mul
         if(tcur.lt.ft(i)) goto 60
         call apply_force_kick(i, dtf)
         gx(i) = gx(i) + vx(i)*dtf
         gy(i) = gy(i) + vy(i)*dtf
         gz(i) = gz(i) + vz(i)*dtf
 60      continue
      enddo

      tcur = tcur + dtf
      goto 50

      end
