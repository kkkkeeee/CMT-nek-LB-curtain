c-----------------------------------------------------------------------
      subroutine nek_init(intracomm)
c

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'ZPER'
c
      include 'OPCTR'
      include 'CTIMER'

C     used scratch arrays
C     NOTE: no initial declaration needed. Linker will take 
c           care about the size of the CBs automatically
c
c      COMMON /CTMP1/ DUMMY1(LCTMP1)
c      COMMON /CTMP0/ DUMMY0(LCTMP0)
c
c      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
c      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
c      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)
  
      common /rdump/ ntdump

      real kwave2
      real*8 t0, tpp

      logical ifemati,ifsync_
      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload

      gfirst=1
      inoassignd=0

      call get_session_info(intracomm)

      etimes = dnekclock()
      istep  = 0
      tpp    = 0.0

      call opcount(1)

      call initdim         ! Initialize / set default values.
      call initdat
      call files

      etime = dnekclock()
      call readat          ! Read .rea +map file
      etims0 = dnekclock_sync()
      if (nio.eq.0) then
         write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,12) 'lx1  /lx2  /lx3 :',lx1,lx2,lx3
         write(6,'(A,g13.5,A,/)')  ' done :: read .rea file ',
     &                             etims0-etime,' sec'
 12      format(1X,A,4I12,/,/)
      endif 

      ifsync_ = ifsync
      ifsync = .true.

      call setvar          ! Initialize most variables

      instep=1             ! Check for zero steps
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

      igeom = 2
      call setup_topo      ! Setup domain topology  

      call genwz           ! Compute GLL points, weights, etc.

      call io_init         ! Initalize io unit

      if(nio.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat' 

      call gengeom(igeom)  ! Generate geometry, after usrdat 

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)

      if(nio.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat2' 

      call geom_reset(1)    ! recompute Jacobians, etc.
      call vrdsmsh          ! verify mesh topology

      call setlog  ! Initalize logical flags

      call bcmask  ! Set BC masks for Dirichlet boundaries.

      if (fintim.ne.0.0.or.nsteps.ne.0) 
     $   call geneig(igeom) ! eigvals for tolerances

      call vrdsmsh     !     Verify mesh topology

      call dg_setup    !     Setup DG, if dg flag is set.

c     commented by keke, set_overlap takes too much time
c     if (ifflow.and.(fintim.ne.0.or.nsteps.ne.0)) then    ! Pressure solver 
c        call estrat                                       ! initialization.
c        if (iftran.and.solver_type.eq.'itr') then         ! Uses SOLN space 
c           call set_overlap                               ! as scratch!
c        elseif (solver_type.eq.'fdm'.or.solver_type.eq.'pdm')then
c           ifemati = .true.
c           kwave2  = 0.0
c           if (ifsplit) ifemati = .false.
c           call gfdm_init(nx2,ny2,nz2,ifemati,kwave2)
c        elseif (solver_type.eq.'25D') then
c           call g25d_init
c        endif
c     endif

      if(ifcvode) call cv_setsize

      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

#ifdef CMTNEK
        call nek_cmt_init
        if (nio.eq.0) write(6,*)'Initialized DG machinery'
#endif

      gfirst = 0       !     for correct executaion of reinitialize call in particles code (userchk)

      call setics      !     Set initial conditions 
      call setprop     !     Compute field properties

      if (instep.ne.0) then !USRCHK
        if(nio.eq.0) write(6,*) 'call userchk'
         if (ifneknek) call userchk_set_xfer
         if (ifneknek) call bcopy
         if (ifneknek) call chk_outflow
         call userchk
         if(nio.eq.0) write(6,'(A,/)') ' done :: userchk' 
      endif

      if (ifcvode .and. nsteps.gt.0) call cv_init

      call comment
      call sstest (isss) 

      call dofcnt

      jp = 0            ! Set perturbation field count to 0 for baseline flow

      call in_situ_init()

      call time00       !     Initalize timers to ZERO
      call opcount(2)

      ntdump=0
      if (timeio.ne.0.0) ntdump = int( time/timeio )

      etims0 = dnekclock_sync()
      if (nio.eq.0) then
        write (6,*) ' '
        if (time.ne.0.0) write (6,'(a,e14.7)') ' Initial time:',time
        write (6,'(a,g13.5,a)') 
     &              ' Initialization successfully completed ',
     &              etims0-etimes, ' sec'
      endif

      ifsync = ifsync_ ! restore initial value

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'PARALLEL'


      real*4 papi_mflops
      integer*8 papi_flops
      integer modstep
      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload
c     integer reinit_step  !added by keke
      integer counter !added by keke
      integer last_kstep !added by keke
      real diff_time, diff_time2, reinit_interval
      real timet 
      integer adaptivelb, stepvalue, rebal
c     real starttime


      call nekgsync()
      reinit_step=0
      diff_time = 0.0
      diff_time2 = 0.0
      counter = 0
      last_kstep = 0


      if (instep.eq.0) then
        if(nid.eq.0) write(6,'(/,A,/,A,/)') 
     &     ' nsteps=0 -> skip time loop',
     &     ' running solver in post processing mode'
      else
        if(nio.eq.0) write(6,'(/,A,/)') 'Starting time loop ...'
      endif

      isyc  = 0
      if(ifsync) isyc=1
      itime = 0
#ifdef TIMER
      itime = 1
#endif
      call nek_comm_settings(isyc,itime)
      call nek_comm_startstat()

      istep  = 0
      msteps = 1

c     Upon redistribution of elements, the following stages happen in sequence:
c     1) reinitialize is called
c     2) elements are moved to different processors to achieve load balancing based on current distribution of particles (particles have not been moved yet)
c     3) Fluid forces are updated in nek__multi_advance
c     4) CMT-bone: Particles are moved based on predefined forces
c     5) Particles are moved to new processors based on which element they are in.
c     To add support of particle movement in CMT-nek, we *must* move particles along
c     with the elements to different processors as part of the load balancing in step (2)
c     above. Then the particle locations are updated based on the new fluid force calculations done in step(3).
      do kstep=1,nsteps,msteps
         timet = DNEKCLOCK()
c        starttime = DNEKCLOCK()
         call nek__multi_advance(kstep,msteps)
c        if(nid.eq. 0) print *, 'nek__multi_advance',
c    $      DNEKCLOCK()-starttime 
c        starttime = DNEKCLOCK()
         call check_ioinfo  
c        if(nid.eq. 0) print *, 'check_ioinfo',
c    $      DNEKCLOCK()-starttime 
c        starttime = DNEKCLOCK()

         call set_outfld
c        if(nid.eq. 0) print *, 'set_outfld',
c    $      DNEKCLOCK()-starttime 
c        starttime = DNEKCLOCK()

         call userchk
c        if(nid.eq. 0) print *, 'userchk',
c    $      DNEKCLOCK()-starttime 
c        starttime = DNEKCLOCK()

         call prepost (ifoutfld,'his')
c        if(nid.eq. 0) print *, 'prepost',
c    $      DNEKCLOCK()-starttime 
c        starttime = DNEKCLOCK()

         call in_situ_check()
c        if(nid.eq. 0) print *, 'in_situ_check()',
c    $      DNEKCLOCK()-starttime 
c        starttime = DNEKCLOCK()

         resetFindpts = 0
         if (lastep .eq. 1) goto 1001

         adaptivelb = param(77)
         if(nid .eq. 0) print *, 'adaptivelb', adaptivelb
         if (adaptivelb .eq. 0) then
             stepvalue = param(78)
             modstep = mod(kstep, stepvalue)
             if (modstep .eq. 0) then
                resetFindpts = 1
                call reinitialize
                !call printVerify
             endif
         else if(adaptivelb .eq. 1) then 
c           auto load balancing
            if(nid .eq. 0) then
               if(kstep .le. reinit_step+10) then !for the first 10 step after
                                            !rebalance, pick the minimum
                                            !one as the init_time
                  if((INIT_TIME .gt. TTIME_STP) 
     $                           .and. (TTIME_STP .ne. 0)) then
                      INIT_TIME = TTIME_STP
                  endif
               else if(kstep .gt. reinit_step+100) then
                  diff_time = (TTIME_STP-INIT_TIME)/INIT_TIME
                  if(nid .eq. 0) then
                     print *, "nid:", nid, "ttime_stp:", TTIME_STP, 
     $                                INIT_TIME, diff_time
                  endif
               endif
            endif

            call bcast(diff_time, 8)
            if (diff_time .gt. 0.3) then
               if (last_kstep .eq. 0) then
                   counter = counter + 1
               else if((counter .le. 2) .and.
     $                     (last_kstep .eq. kstep-1))then
                   counter = counter + 1
               else
                   counter = 0
               endif
               last_kstep = kstep
               if (counter .gt. 2) then
                   !print *, "into the reinit, nid:", nid, "diff_time:",
     $            !diff_time
                   resetFindpts = 1
                   call reinitialize
                   !call printVerify
                   reinit_step = kstep
                   if(nid .eq. 0) then
                      print *, "reintilize, reinit_step:", reinit_step
                   endif
                   diff_time = 0
                   INIT_TIME = 100
                   counter = 0
               endif
            endif
         else if (adaptivelb .eq. 2) then !for the new adaptive lb algorithm
            if(nid .eq. 0 .and. reinit_step .eq. 0) then
               if(kstep .le. reinit_step+100) then !for the first 10 step after
                                            !rebalance, pick the minimum
                                            !one as the init_time
                  if((INIT_TIME .gt. TTIME_STP)
     $                           .and. (TTIME_STP .ne. 0)) then
                      INIT_TIME = TTIME_STP
                  endif
               else if(kstep .gt. reinit_step+100) then
                  diff_time = (TTIME_STP-INIT_TIME)/INIT_TIME
                  diff_time2 = TTIME_STP-INIT_TIME
                  if(nid .eq. 0) then
                     print *, "nid:", nid, "ttime_stp:", TTIME_STP,
     $                                INIT_TIME, diff_time
                  endif
               endif
            endif

            call bcast(diff_time, 8)
            call bcast(diff_time2, 8)
            if(reinit_step .ne. 0) then
               if(kstep-reinit_step .eq. 1) then
                   lb_time=ttime_stp !added by keke for adaptive lb
                   if(nid .eq. 0) print *,'istep',istep,'reinit_step',
     :                     reinit_step,'lb_time', lb_time
                  call bcast(lb_time,8)
                  rebal=int(sqrt(2*reinit_interval*lb_time/diff_time2))
                  if(nid .eq. 0) then
                     print *, 'rebal:', rebal, 'lb_time:',lb_time
     $                ,'reinit_interval', reinit_interval
     $                ,'diff_time2', diff_time2,kstep, nid
                  endif
               else if(kstep .le. reinit_step+10) then !for the first 10 step after
                                            !rebalance, pick the minimum
                                            !one as the init_time
                  if((INIT_TIME .gt. TTIME_STP)
     $                           .and. (TTIME_STP .ne. 0)) then
                      INIT_TIME = TTIME_STP
                  endif
               else if(kstep-reinit_step.ge.rebal)then
                  diff_time2 = TTIME_STP-INIT_TIME
                  call bcast(diff_time2, 8)
                  if(diff_time2 .gt. 0) then 
                     resetFindpts = 1
                     call reinitialize
                     reinit_interval = kstep - reinit_step
                     reinit_step = kstep
                     INIT_TIME = 100
                     if(nid .eq. 0) then
                        print *, "reintilize, reinit_step:", reinit_step
                     endif
                  endif
               endif
            endif
            if (diff_time .gt. 0.2 .and. reinit_step .eq. 0) then
               if (last_kstep .eq. 0) then
                   counter = counter + 1
               else if((counter .le. 2) .and.
     $                     (last_kstep .eq. kstep-1))then
                   counter = counter + 1
               else
                   counter = 0
               endif
               last_kstep = kstep
               if (counter .gt. 2) then 
                   !print *, "into the reinit, nid:", nid, "diff_time:",
     $            !diff_time
                   resetFindpts = 1
                   call reinitialize
                   !call printVerify
                   reinit_interval = kstep - reinit_step
                   reinit_step = kstep
                   diff_time = 0
                   INIT_TIME = 100
                   counter = 0
                   if(nid .eq. 0) then
                      print *, "first reintilize, reinit_step:"
     $                   , reinit_step
                   endif
               endif
            endif
         else if (adaptivelb .eq. 3) then
            call adaptive_loadbalance(kstep) 
         endif 
      enddo
 1001 lastep=1


      call nek_comm_settings(isyc,0)

      call comment


c     check for post-processing mode
      if (instep.eq.0) then
         nsteps=0
         istep=0
         if(nio.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nio.eq.0) write(6,*) 'done :: userchk'
         call prepost (.true.,'his')
      else
         if (nio.eq.0) write(6,'(/,A,/)') 
     $      'end of time-step loop' 
      endif


      RETURN
      END

c-----------------------------------------------------------------------
      subroutine adaptive_loadbalance(kstep)
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'PARALLEL'

      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload

      integer kstep
      integer c1_step
      integer c2_step
      integer tt, timesLoadBalance, nmax
      real lb_overhead, m, runtime
      real start_time, lb_time2, t2_time
      real t1_time
      save t1_time
      data t1_time / 0.0/
      integer lb_interval
      save lb_interval
      data lb_interval /1000000/ 

      if (lb_interval .eq. 1000000 .and. nid .eq. 0) then
          if (kstep .gt. 3 .and. kstep .lt. 14) then
              t1_time = t1_time + TTIME_STP
          endif  
          if (kstep .eq. 14) then
              t1_time = t1_time/10.0
              c1_step = 7
          endif  
          if(kstep .gt. 14) then
              diff_time = (TTIME_STP-t1_time)/t1_time
              print *, "nid:", nid, "ttime_stp:", TTIME_STP,
     $                                t1_time, diff_time
          endif
      endif 

      if (lb_interval .eq. 1000000) then
          call bcast(diff_time, 8) 
          if (diff_time .gt. 0.2) then
              c2_step = kstep 
              t2_time = TTIME_STP
              m = (t2_time - t1_time) / (c2_step - c1_step)
              start_time = dnekclock()
              resetFindpts = 1
              call reinitialize
              lb_time2 = dnekclock() - start_time
              if (nid .eq. 0) then
                   lb_overhead = lb_time2/t1_time
                   tt = nsteps - c2_step
                   nmax = int(ceiling(tt/lb_overhead))
                   timesLoadBalance = 0
                   runtime = 1.0/2.0 * tt * (tt * m)
                   do i=1, nmax
                       r=lb_time2*i + (i+1)*0.5*(tt/(i+1))*(tt/(i+1))*m
                       if (r < runtime) then
                           runtime = r
                           timesLoadBalance = i
                       else
                           exit
                       endif
                   enddo 
                   lb_interval = tt/timesLoadBalance
              endif
              call bcast(lb_interval, 4) 
              if (nid .eq. 0) then
                  print * , "lb_interval", lb_interval, c2_step
              endif
          endif
      endif
      if (kstep - c2_step .eq. lb_interval) then
          resetFindpts = 1
          call reinitialize
          c2_step = kstep
      endif
      end
c-----------------------------------------------------------------------
      subroutine nek_advance

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      real timet

      common /cgeom/ igeom

      ntot = nx1*ny1*nz1*nelv

      call nekgsync

      call setup_convect(2) ! Save conv vel

      if (iftran) call settime
      if (ifmhd ) call cfl_check
      call setsolv
      call comment

#ifdef CMTNEK
      if (nio.eq.0.and.istep.le.1) write(6,*) 'CMT branch active'
      call cmt_nek_advance
      return
#endif


      if (ifsplit) then   ! PN/PN formulation


         do igeom=1,ngeom


         ! within cvode we use the lagged wx for 
         ! extrapolation, that's why we have to call it before gengeom 
         if (ifheat .and. ifcvode) call heat_cvode (igeom)   

         if (ifgeom) then
            call gengeom (igeom)
            call geneig  (igeom)
         endif


         if (ifheat) call heat (igeom)

         if (igeom.eq.2) then  
            call setprop
            call rzero(qtl,ntot)
            if (iflomach) call qthermal
         endif


         if (ifflow)          call fluid    (igeom)
         if (ifmvbd)          call meshv    (igeom)
         if (param(103).gt.0) call q_filter (param(103))

         enddo

      else                ! PN-2/PN-2 formulation

         call setprop
         do igeom=1,ngeom

            if (igeom.gt.2) call userchk_set_xfer

            if (ifgeom) then
               call gengeom (igeom)
               call geneig  (igeom)
            endif

            if (ifneknekm.and.igeom.eq.2) call multimesh_create

            if (ifmhd) then
               if (ifheat)      call heat     (igeom)
                                call induct   (igeom)
            elseif (ifpert) then
               if (ifbase.and.ifheat)  call heat          (igeom)
               if (ifbase.and.ifflow)  call fluid         (igeom)
               if (ifflow)             call fluidp        (igeom)
               if (ifheat)             call heatp         (igeom)
            else  ! std. nek case
               if (ifheat)             call heat          (igeom)
               if (ifflow)             call fluid         (igeom)
               if (ifmvbd)             call meshv         (igeom)
            endif

            if (igeom.eq.ngeom.and.param(103).gt.0) 
     $          call q_filter(param(103))
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine nek_end

      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'OPCTR'

      if(instep.ne.0)  call runstat
      if(xxth(1).gt.0) call crs_stats(xxth(1))

   
      call in_situ_end()
      return
      end
c-----------------------------------------------------------------------
      subroutine nek__multi_advance(kstep,msteps)

      include 'SIZE'
      include 'TOTAL'
      real timet

      do i=1,msteps
         istep = istep+i
         timet = DNEKCLOCK()
         call nek_advance
c        if(nid .eq . 0) print *, "nek_advance",
c    $          DNEKCLOCK()-timet
c        timet = DNEKCLOCK()

         if (ifneknek) call userchk_set_xfer
c        if(nid .eq . 0) print *, "userchk_set_xfer",
c    $          DNEKCLOCK()-timet
c        timet = DNEKCLOCK()

         if (ifneknek) call bcopy
c        if(nid .eq . 0) print *, "bcopy",
c    $          DNEKCLOCK()-timet
c        timet = DNEKCLOCK()

         if (ifneknek) call chk_outflow
c        if(nid .eq . 0) print *, "chk_outflow",
c    $          DNEKCLOCK()-timet
c        timet = DNEKCLOCK()


      enddo

      return
      end
c-----------------------------------------------------------------------
