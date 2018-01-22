c-----------------------------------------------------------------------
      subroutine readat_lb
C
C     Read in data from preprocessor input file (.rea)
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'ZPER'

      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload

      logical ifbswap,ifre2,parfound
      character*132 string
      integer idum(3*numsts+3)
      ierr = 0
      call flush_io

      ! check if new rea file version exists
      if(nid.eq.0) inquire(file=parfle, exist=parfound)
      call bcast(parfound,lsize)
      if (parfound) then
         if(nio.eq.0) write(6,'(A,A)') ' Reading ', parfle
         call readat_new
         return
      endif

      call nekgsync()
      starttime = dnekclock()
      if (.not.ifgtp) call mapelpr  ! read .map file, est. gllnid, etc.
      endtime = dnekclock()
      if(nid .eq. 0) print *, "mapelpr time ", endtime-starttime

      return
      END
c------------------------------------------------------------------------
      subroutine usrdat2_lb
      include 'SIZE'
      include 'TOTAL'
      include 'TORO'
      include 'CMTBCDATA'
      include 'CMTDATA'
      include 'PERFECTGAS'

      outflsub=.true.
      IFCNTFILT=.false.
      ifrestart=.false.
      ifsip=.false.
      gasmodel = 1
! JH080714 Now with parameter space to sweep through
      starttime1 = dnekclock_sync()
      open(unit=81,file="riemann.inp",form="formatted")
      read (81,*) domlen
      read (81,*) xdiaph
      read (81,*) gmaref
      read (81,*) dleft
      read (81,*) uleft
      read (81,*) pleft
      read (81,*) dright
      read (81,*) uright
      read (81,*) pright
      read (81,*) zerotime
      close(81)

!     molmass=8314.3
!     muref=0.0
!     coeflambda=-2.0/3.0
!     suthcoef=1.0
!     prlam = 0.72
!      rgasref    = MixtPerf_R_M(molmass,dum)
!      cvgref     = rgasref/(gmaref-1.0)
!      cpgref     = MixtPerf_Cp_CvR(cvgref,rgasref)
!      gmaref     = MixtPerf_G_CpR(cpgref,rgasref)

      c_max=0.5     ! should be 0.5, really
      c_sub_e=1.0e36
      endtime = dnekclock_sync()
      if( mod(nid, np/2) .eq. np/2-2) then
         print *, 'usrdat2_check ', endtime-starttime1
      endif

!      call e1rpex(domlen,xdiaph,gmaref,dleft,uleft,pleft,dright,uright,
!     >            pright,1.0)
c     CALL SAMPLE(PMstar, UM, 0.0, rhohere, uhere, pinfty)
!      reftemp=pleft/dleft/rgasref
!      aleft=sqrt(gmaref*pleft/dleft)
c     call domain_size(xmin,xmax,ymin,ymax,zmin,zmax)
c     if(nio.eq.0)then
c        write(6,*) 'domlen',domlen
c        write(6,*) 'xdiaph',xdiaph
c        write(6,*) 'gamma ',gmaref
c        write(6,*) 'rhol  ',dleft
c        write(6,*) 'ul    ',uleft
c        write(6,*) 'pl    ',pleft
c        write(6,*) 'rhor  ',dright
c        write(6,*) 'ur    ',uright
c        write(6,*) 'pr    ',pright
c        write(6,*) 'sound ',aleft
c        write(6,*) 'ustar ',um
c        write(6,*) 'dt    ',dt
c        write(6,*) 'nsteps',nsteps
c        write(6,*) 'final time ',(xdiaph-xmin)/aleft
c     endif
      return
      end
!-----------------------------------------------------------------------
      subroutine nek_cmt_init_lb
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'CMTDATA'
      if (nio.eq.0) write(6,*)'Set up CMT-Nek'
      if (toteq.ne.5) then
         if (nio.eq.0) write(6,*)'toteq is low ! toteq = ',toteq
         if (nio.eq.0) write(6,*) 'Reset toteq in SIZE to 5'
         call exitt
      endif
      if (ifrestart) then
         ifheat = .true. ! almost certainly incorrect
      endif
      call setup_cmt_commo_lb

      return
      end
!-----------------------------------------------------------------------

      subroutine setup_cmt_commo_lb
      include 'SIZE'
      include 'TOTAL'
      include 'DG'

      parameter(lf=lx1*lz1*2*ldim*lelt)
      common /c_is1/ glo_num_face(lf)
     $             , glo_num_vol((lx1+2)*(ly1+2)*(lz1+2)*lelt)
      integer*8 glo_num_face,glo_num_vol,ngv

      call setup_cmt_gs(dg_hndl,nx1,ny1,nz1,nelt,nelgt,vertex,
     >                  glo_num_vol,glo_num_face)
      call cmt_set_fc_ptr(nelt,nx1,ny1,nz1,ndg_face,iface_flux)

      return
      end

!-----------------------------------------------------------------------

      subroutine setup_cmt_gs_lb(dg_hndl,nx,ny,nz,nel,melg,vertex,gnv
     >     ,gnf)

!     Global-to-local mapping for gs

      include 'SIZE'
      include 'TOTAL'

      integer   dg_hndl
      integer   vertex(*)

      integer*8 gnv(*),gnf(*),ngv

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      mx = nx+2
      call set_vert(gnv,ngv,mx,nel,vertex,.false.) ! lives in navier8.f

      mz0 = 1
      mz1 = 1
      if (if3d) mz0 = 0
      if (if3d) mz1 = nz+1
      call iface_vert_int8cmt(nx,ny,nz,gnf,gnv,mz0,mz1,nelt)

      nf = nx*nz*2*ndim*nelt !total number of points on faces BETTA BE 4-byte!
      call gs_setup(dg_hndl,gnf,nf,nekcomm,np)

      return
      end

!-----------------------------------------------------------------------

      subroutine cmt_set_fc_ptr_lb(nel,nx,ny,nz,nface,iface)

!     Set up pointer to restrict u to faces ! NOTE: compact
! JH062314 Now 2D so we can strip faces by element and not necessarily
!          from the whole field

      include 'SIZE'
      include 'TOTAL'

      integer nx, ny, nz, nel
      integer nface,iface(nx*nz*2*ldim,*)
      integer e,f,ef

      call dsset(nx,ny,nz) ! set skpdat. lives in connect1.f

      nxyz = nx*ny*nz
      nxz  = nx*nz
      nfpe = 2*ndim
      nxzf = nx*nz*nfpe ! red'd mod to area, unx, etc.

      do e=1,nel
      do f=1,nfpe

         ef     = eface(f)
         js1    = skpdat(1,f)
         jf1    = skpdat(2,f)
         jskip1 = skpdat(3,f)
         js2    = skpdat(4,f)
         jf2    = skpdat(5,f)
         jskip2 = skpdat(6,f)

         i = 0
         do j2=js2,jf2,jskip2
         do j1=js1,jf1,jskip1

            i = i+1
            k = i+nxz*(ef-1)           ! face numbering
            iface(k,e) = j1+nx*(j2-1)  ! cell numbering

         enddo
         enddo

      enddo
      enddo
      nface = nxzf*nel

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_topo_lb
C
C     Parallel compatible routine to find 
C     connectivity of element structure.
C
C     On Processor 0:
C
C     .Verify right-handedness of elements.
C     .Verify element-to-element reciprocity of BC's
C     .Verify correlation between E-E BC's and physical coincidence
C     .Set rotations
C     .Determine multiplicity
C     .Set up direct stiffness summation arrays.
C
C     All Processors:
C
C     .Disperse/Receive BC and MULT temporary data read from
C     preprocessor.
C
C
      include 'SIZE'
      include 'TOTAL'
      include 'NONCON'
      include 'ZPER'
      include 'SCRCT'
c
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C
      common /c_is1/ glo_num(1*lx1*ly1*lz1*lelv)
      integer*8 glo_num
      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer vertex

      if(nio.eq.0) write(6,*) 'setup mesh topology'
C
C     Initialize key arrays for Direct Stiffness SUM.
C
      NXL=3
      NYL=3
      NZL=1+2*(NDIM-2)

      !call initds
      !call dsset (nxl,nyl,nzl)
      !call setedge
C
C=================================================
C     Establish (global) domain topology
C=================================================
C
C     .Generate topologically correct mesh data.
C     .Set up element centers, face centers, etc. 
C     .Check  right handedness of elements.
C     .Check  element boundary conditions.
C     .Establish Element-Element rotations
C     .Construct the element to processor map and

      !call genxyzl
      call setside
      !call verify

      !CALL SETCDOF
      !IF (IFAXIS            ) CALL SETRZER
      !IF (IFMVBD            ) CALL CBCMESH
      !IF (IFMODEL.AND.IFKEPS) CALL CBCTURB
      !CALL CHKAXCB
C
C========================================================================
C     Set up element-processor mapping and establish global numbering
C========================================================================
C
      mfield=2
      if (ifflow) mfield=1
      if (ifmvbd) mfield=0

      ncrnr = 2**ndim

      if (nelgv.eq.nelgt) then
         !if (ifgtp) then
         !   call gen_gtp_vertex    (vertex, ncrnr)
         !else
         !   call get_vert
         !endif
         call setupds(gsh_fld(1),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)
         gsh_fld(2)=gsh_fld(1)
      else

c
c        For conjugate heat transfer, it is assumed that fluid
c        elements are listed both globally and locally with lower
c        element numbers than the solid elements.
c
c        We currently assume that there is at least one fluid elem.
c        per processor.
c

         call get_vert
c        call outmati(vertex,4,nelv,'vrtx V')
         call setupds(gsh_fld(1),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)

c        call get_vert  (vertex, ncrnr, nelgt, '.mp2')  !  LATER !
c        call outmati(vertex,4,nelt,'vrtx T')
         call setupds(gsh_fld(2),nx1,ny1,nz1,nelt,nelgt,vertex,glo_num)
      endif
C========================================================================
C     Set up multiplicity and direct stiffness arrays for each IFIELD
C========================================================================

      ntotv = nx1*ny1*nz1*nelv
      ntott = nx1*ny1*nz1*nelt


c     print *, 'ifflow', ifflow, ifheat, ifmvbd !T F F 
      if (ifflow) then
         ifield = 1
         call rone    (vmult,ntotv)
         call dssum   (vmult,nx1,ny1,nz1)
         vmltmax=glmax(vmult,ntotv)
         ivmltmax=vmltmax
         if (nio.eq.0) write(6,*) ivmltmax,' max multiplicity'
         call invcol1 (vmult,ntotv)
      endif
      if (ifheat) then
         ifield = 2
         call rone    (tmult,ntott)
         call dssum   (tmult,nx1,ny1,nz1)
         call invcol1 (tmult,ntott)
      endif
      if (.not.ifflow) call copy(vmult,tmult,ntott)
      if (ifmvbd)  call copy (wmult,vmult,ntott)
c     do ifield=3,nfield                  ! Additional pass. scalrs.
c        if (nelg(ifield).eq.nelgv) then
c           gsh_fld(ifield) = gsh_fld(1)
c           call copy (tmult(1,1,1,1,ifield-1),vmult,ntotv)
c        else
c           gsh_fld(ifield) = gsh_fld(2)
c           call copy (tmult(1,1,1,1,ifield-1),tmult,ntott)
c        endif
c     enddo

      ifgsh_fld_same = .true.
      do ifield=2,nfield
         if (gsh_fld(ifield).ne.gsh_fld(1)) then
            ifgsh_fld_same = .false.
            goto 100
         endif
      enddo
 100  continue
      if(nio.eq.0) then
        write(6,*) 'done :: setup mesh topology'
        write(6,*) ' '
      endif

      return
      end


