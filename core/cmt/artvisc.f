      subroutine compute_entropy(s)
! computes entropy at istep and pushes the stack down for previous
! steps needed to compute ds/dt via finite difference (for now).
! hardcoded for Burgers equation. More later when folded into CMT-nek
! for Burgers, s=energy=1/2 U^2
      include 'SIZE'
      include 'TOTAL'  ! tlag is lurking. be careful
      include 'CMTDATA'
! I've always seen lorder=3, but I still need three steps
!          s(:,               1,       1)  ! entropy at current step
!          s(:,               2,       1)  ! entropy at step n-1
!          s(:,               1,       2)  ! entropy at step n-2
      real s(lx1*ly1*lz1*lelt,lorder-1,*)
      real ntol
      integer e
      data icalld /0/
      save icalld

      n=nx1*ny1*nz1
      ntot=n*nelt
      ntol=1.0e-10

      if (icalld .eq. 0) then
         if(nid.eq.0)write(6,*) 'zeroing out entropy stack',istep
         icalld=1
         call rzero(s,ntot)
         call rzero(s(1,1,2),ntot) ! s_{n-1}
         call rzero(s(1,2,1),ntot) ! s_n
      endif

! compute the current entropy. This actually needs to go back in the
! usr file because it's EOS-dependent
      rgam=rgasref/(gmaref-1.0)
      do i=1,ntot
         rho=max(vtrans(i,1,1,1,irho),ntol)
         s(i,1,1)=rgam*rho*log(pr(i,1,1,1)/(rho**gmaref))
      enddo

      if (stage .eq. 1) then
! push the stack
         call copy(s(1,1,2),s(1,2,1),ntot) ! s_{n-1}=s_n
         call copy(s(1,2,1),s(1,1,1),ntot) ! fill s_n
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine entropy_viscosity
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      parameter (lxyz=lx1*ly1*lz1)
      common /scrns/ scrent(lxyz,lelt)
      integer e
      character*132 deathmessage

      pi=4.0*atan(1.0)

      n=nx1*ny1*nz1
      ntot=n*nelt

! entropy at this and lorder prior steps
      call compute_entropy(tlag)
! compute maxval(|S-<S>|)
      savg    =    glsc2(tlag,bm1,ntot)
      savg    = -savg/volvm1
      call cadd2(scrent,tlag,savg,ntot)
      maxdiff =     glamax(scrent,ntot)
      if (maxdiff.le.0.0) then
         write(deathmessage,*) 'zero maxdiff usually means NAN$'
         call exittr(deathmessage,maxdiff,istep)
!     else
!        if (nio .eq. 0) write (6,*) 'max(s-<s>)=',maxdiff, meshh(1)
      endif
      call entropy_residual(tlag) ! fill res2
      call copy(res2(1,1,1,1,2),res2,ntot) ! raw residual in res2
      call wavevisc(t(1,1,1,1,3))
      call resvisc(res2) ! overwrite res2
      call evmsmooth(res2,t(1,1,1,1,3),.true.) ! endpoints=.false.
                                               ! is intended to
                                               ! preserve face states,
                                               ! but this is easier to
                                               ! test 1D
!     call evmsmooth(res2,t(1,1,1,1,3),.true.) ! And again.
      call dsavg(res2) ! you DEFINITELY don't want a min here


      return
      end

c-----------------------------------------------------------------------

      subroutine entropy_residual(s)
! COMPUTE R=ds/dt + div F(s) for entropy pair s and (hardcoded) F
! and store its norm functional thing |R| in res2 where it will
! provide artificial viscosity according to the code in
! entropy_viscosity and the method of Guermond
      include 'SIZE'
      include 'TOTAL' ! tlag lurks
      include 'CMTDATA'
      integer e
      real dsdtcoef(3) ! put this in /TIMESTEPCOEF/ someday
      data dsdtcoef /1.0,1.0,0.5/
      real s(lx1*ly1*lz1*lelt,lorder-1,ldimt) ! because it's really tlag

      n=nx1*ny1*nz1
      ntot=n*nelt

      if (istep .eq. 1) return
      rdt=1.0/(dsdtcoef(stage)*DT_cmt)
      if (stage .eq. 1) then ! THE MOST SABOTAGEABLE PART OF THE CODE
         call sub3(res2,s(1,1,1),s(1,1,2),ntot) ! EVALUATE s_n-s_{n-1}
      else
         call sub3(res2,s(1,1,1),s(1,2,1),ntot) ! EVALUATE s^{(stage)}-s_n
      endif
      call cmult(res2,rdt,ntot)
! res2=ds/dt. now,

!-----------------------------------------------------------------------
! cons approach: strong-conservation form of flux divergence in entropy residual
!-----------------------------------------------------------------------
! get around to expanding totalh to store fluxes for whole fields and
! properly vectorize evaluate_*_h and flux_div_integral
      do e=1,nelt
         call evaluate_entropy_flux(e) ! diffh. zero it before diffusion 
         call flux_div_mini(e) ! into res2 it goes.
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine evaluate_entropy_flux(e)
! entropy flux function for entropy residual.
! just vel*s for now
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'CMTDATA'
      integer e

      call rzero(totalh,3*nxd*nyd*nzd)
      n=nx1*ny1*nz1

      call col3(totalh(1,1),vx(1,1,1,e),tlag(1,1,1,e,1,1),n)
      call col3(totalh(1,2),vy(1,1,1,e),tlag(1,1,1,e,1,1),n)
      if (if3d) call col3(totalh(1,3),vz(1,1,1,e),tlag(1,1,1,e,1,1),n)

      return
      end

c-----------------------------------------------------------------------

      subroutine flux_div_mini(e)
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'DXYZ'
      include 'SOLN'
      include 'CMTDATA'
      parameter (ldd=lx1*ly1*lz1)
      parameter (ldg=lx1**3,lwkd=2*ldg)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju

      integer e

      nrstd=ldd
      nxyz=nx1*ny1*nz1
      mdm1=nx1-1

      call rzero(ur,nrstd)
      call rzero(us,nrstd)
      call rzero(ut,nrstd)
      call rzero(ud,nrstd)
      call rzero(tu,nrstd)

      if (if3d) then
         call local_grad3(ur,us,ut,totalh(1,1),mdm1,1,dxm1,dxtm1)
         do i=1,nxyz
            ud(i) = jacmi(i,e)*(rxm1(i,1,1,e)*ur(i)
     $                        + sxm1(i,1,1,e)*us(i)
     $                        + txm1(i,1,1,e)*ut(i))
         enddo
         call local_grad3(ur,us,ut,totalh(1,2),mdm1,1,dxm1,dxtm1)
         do i=1,nxyz ! confirmed to have no effect in 1D
            ud(i)=ud(i)+jacmi(i,e)*(rym1(i,1,1,e)*ur(i)
     $                            + sym1(i,1,1,e)*us(i)
     $                            + tym1(i,1,1,e)*ut(i))
         enddo
         call local_grad3(ur,us,ut,totalh(1,3),mdm1,1,dxm1,dxtm1)
         do i=1,nxyz ! confirmed to have no effect in 1D
            ud(i)=ud(i)+jacmi(i,e)*(rzm1(i,1,1,e)*ur(i)
     $                            + szm1(i,1,1,e)*us(i)
     $                            + tzm1(i,1,1,e)*ut(i))
         enddo
      else
         call local_grad2(ur,us,totalh(1,1),mdm1,1,dxm1,dxtm1)
         do i=1,nxyz
            ud(i) = jacmi(i,e)*(rxm1(i,1,1,e)*ur(i)
     $                        + sxm1(i,1,1,e)*us(i))
         enddo
         call local_grad2(ur,us,totalh(1,2),mdm1,1,dxm1,dxtm1)
         do i=1,nxyz
            ud(i)=ud(i)+jacmi(i,e)*(rym1(i,1,1,e)*ur(i)
     $                            + sym1(i,1,1,e)*us(i))
         enddo
      endif
      call add2(res2(1,1,1,e,1),ud,nxyz)

      return
      end

!-----------------------------------------------------------------------

      subroutine resvisc(residual)
! smooth residual-based entropy visc, defined by Guermond, Popov, whoever
! DSAVG assumes IFFLOW to work correctly. IFHEAT still doesn't work with CMT-nek
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'DG'
      integer lfq,heresize,hdsize
      parameter (lxyz=lx1*ly1*lz1)
      real residual(lxyz,nelt)
      integer e,f

      nxyz =nx1*ny1*nz1
      nxz  =nx1*nz1
      nface=2*ndim
      nxzf =nxz*nface
      nfq  =nxzf*nelt

! ensure continuity at faces. doing this before |abs| causes some cancellation that,
! so far, appears to beneficially reduce spikiness at faces.

      do e=1,nelt
         do i=1,nxyz
            residual(i,e)=residual(i,e)*meshh(e)**2
         enddo
      enddo

      call dsavg(residual) ! signed, can cancel at faces. Hope it does

      do e=1,nelt
         do i=1,nxyz
            residual(i,e)=abs(residual(i,e))
         enddo
      enddo

      call cmult(residual,c_sub_e,nxyz*nelt)

      if (maxdiff .ne. 0) then
         const=1.0/maxdiff
         call cmult(residual,const,nxyz*nelt)
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine evmsmooth(resvisc,wavevisc,endpoints)
      include 'SIZE'
      include 'INPUT'
      real resvisc(nx1,ny1,nz1,nelt),wavevisc(nx1,ny1,nz1,nelt)
      real rtmp
      common /ctmp1/ rtmp(lx1,ly1,lz1)
! are faces included in smoothing? if not (say they're fixed by dsavg) then
! endpoints should be .false.
      logical endpoints
      integer e
! clip residual viscosity and smooth it.
! actually just smooth first and then clip it. Just average-with-my-neighbors for now
! weight the neighbors according to GLL-spacing instead of uniform-grid
! at a later date.

      nxyz=nx1*ny1*nz1
      kstart=1
      kend=nz1
      rndim=1.0/ndim

      if (endpoints) then
         istart=1
         jstart=1
         iend=nx1
         jend=ny1
      else
         istart=2
         jstart=2
         iend=nx1-1
         jend=ny1-1
         if (if3d) then
            kstart=2
            kend=nz1-1
         endif
      endif

! yes I know this loop sucks and I need to write a matrix or something
      do e=1,nelt
         do i=1,nxyz
            if (wavevisc(i,1,1,e).le. resvisc(i,1,1,e))
     >          resvisc(i,1,1,e)= wavevisc(i,1,1,e)
         enddo
         call copy(rtmp,resvisc(1,1,1,e),nxyz) ! really only for .false.
         do iz=kstart,kend
            if (if3d) then
               km1=iz-1
               kp1=iz+1
               izm=km1
!              if (km1 .lt. 1) izm=iz ! bias towards {{face point}}
               if (km1 .lt. 1) izm=kp1 ! Guermond symmetry
               izp=kp1
!              if (kp1 .gt. nz1) izp=iz ! bias towards {{face point}}
               if (kp1 .gt. nz1) izp=km1 ! Guermond symmetry
            else
               izm=iz
               izp=iz
            endif
            do iy=jstart,jend
               jm1=iy-1
               jp1=iy+1
               iym=jm1
!              if (jm1 .lt. 1) iym=iy ! bias towards {{face point}}
               if (jm1 .lt. 1) iym=jp1 ! Guermond symmetry
               iyp=jp1
!              if (jp1 .gt. ny1) iyp=iy ! bias toward {{face point}}
               if (jp1 .gt. ny1) iyp=jm1 ! Guermond symmetry
               do ix=istart,iend
                  im1=ix-1
                  ip1=ix+1
                  ixm=im1
!                 if (im1 .lt. 1) ixm=ix ! bias towards {{face point}}
                  if (im1 .lt. 1) ixm=ip1 ! Guermond symmetry
                  ixp=ip1
!                 if (ip1 .gt. nx1) ixp=ix ! bias towards {{face point}}
                  if (ip1 .gt. nx1) ixp=im1 ! Guermond symmetry
                  x0 = resvisc(ix ,iy ,iz ,e)
                  x1 = resvisc(ixm,iy ,iz ,e)
                  x2 = resvisc(ixp,iy ,iz ,e)
                  x3 = resvisc(ix ,iym,iz ,e)
                  x4 = resvisc(ix ,iyp,iz ,e)
                  if (if3d) then
                     x5 = resvisc(ix ,iy ,izm,e)
                     x6 = resvisc(ix ,iy ,izp,e)
                  else
                     x5=0.0
                     x6=0.0
                  endif
                  rtmp(ix,iy,iz)=0.25*(2.0*ndim*x0+x1+x2+x3+x4+x5+x6)
     >                               *rndim
               enddo
            enddo
         enddo
         call copy(resvisc(1,1,1,e),rtmp,nxyz)
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine wavevisc(numax)
! max entropy visc, defined by Guermond, Popov, whoever (chapter-verse)
! as 
! numax = c_max*h*max(dH/dU)
! in a given element
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      parameter (lxyz=lx1*ly1*lz1)
      common /scrns/ wavespeed(lxyz)
      real wavespeed
      real maxeig
      real numax(lxyz,nelt)
      integer e

      nxyz=nx1*ny1*nz1

      do e=1,nelt
         do i=1,nxyz
            wavespeed(i)=csound(i,1,1,e)+
     >      sqrt(vx(i,1,1,e)**2+vy(i,1,1,e)**2+vz(i,1,1,e)**2)
         enddo
         maxeig=vlamax(wavespeed,nxyz)
         rhomax(e)=vlamax(vtrans(1,1,1,e,irho),nxyz)
         do i=1,nxyz
            numax(i,e)=c_max*maxeig*meshh(e)
         enddo
      enddo

      call max_to_trilin(numax)

      return
      end

!-----------------------------------------------------------------------

      subroutine max_to_trilin(field)
! stupid subroutine to take a stupid uniform field and compute a trilinear
! tent between maximum shared values at the vertices.
      include 'SIZE'
      include 'TOTAL'
      real field(lx1,ly1,lz1,nelt)
      integer e

      character*32 fname ! diagnostic

      nxyz=nx1*ny1*nz1

! get maxima on faces
      call dsop(field,'MAX',nx1,ny1,nz1)

! trilinear interpolation. you should adapt xyzlin to your needs instead
      do e=1,nelt
         p000=field(1,  1,  1,  e)
         p100=field(nx1,1,  1,  e)
         p010=field(1,  ny1,1,  e)
         p110=field(nx1,ny1,1,  e)
         p001=field(1,  1,  nz1,e)
         p101=field(nx1,1,  nz1,e)
         p011=field(1,  ny1,nz1,e)
         p111=field(nx1,ny1,nz1,e)
         c1=p100-p000
         c2=p010-p000
         c3=p001-p000
         c4=p110-p010-p100+p000
         c5=p011-p001-p010+p000
         c6=p101-p001-p100+p000
         c7=p111-p011-p101-p110+p100+p001+p010-p000
         rdx=1.0/(xm1(nx1,1,1,e)-xm1(1,1,1,e)) ! cubes only!!!
         rdy=1.0/(ym1(1,ny1,1,e)-ym1(1,1,1,e))
         rdz=0.0
         if(if3d) rdz=1.0/(zm1(1,1,nz1,e)-zm1(1,1,1,e))
         do i=1,nxyz
            deltax=rdx*(xm1(i,1,1,e)-xm1(1,1,1,e)) ! cubes only!!!
            deltay=rdy*(ym1(i,1,1,e)-ym1(1,1,1,e))
            deltaz=0.0
            if (if3d) deltaz=rdz*(zm1(i,1,1,e)-zm1(1,1,1,e))
            field(i,1,1,e)=p000+c1*deltax+c2*deltay+c3*deltaz+
     >                          c4*deltax*deltay+c5*deltay*deltaz+
     >                       c6*deltaz*deltax+c7*deltay*deltaz*deltax
         enddo
      enddo

      return
      end
