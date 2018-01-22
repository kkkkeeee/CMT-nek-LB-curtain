c-----------------------------------------------------------------------
      subroutine mapelpr()
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SCRCT'
      include 'SOLN'
      include 'TSTEP'
      real starttime, starttime1, endtime
c
      logical ifverbm
c
      if(nio.eq.0) write(6,'(/,A)') ' mapping elements to processors'

      MFIELD=2
      IF (IFFLOW) MFIELD=1
      IF (IFMVBD) MFIELD=0

c     Set up TEMPORARY value for NFIELD - NFLDT
      NFLDT = 1
      IF (IFHEAT) NFLDT = 2 + NPSCAL
c
c     Distributed memory processor mapping
      IF (NP.GT.NELGT) THEN
         IF(NID.EQ.0) THEN
           WRITE(6,1000) NP,NELGT
 1000      FORMAT(2X,'ABORT: Too many processors (',I8
     $          ,') for to few elements (',I8,').'
     $          ,/,2X,'ABORTING IN MAPELPR.')
         ENDIF
         call exitt
      ENDIF
c     endtime = dnekclock_sync()
c     if(nid.eq.0) print *, "mapelpr before set_proc", endtime-starttime 
      starttime1 = dnekclock_sync()
      call set_proc_map()
      endtime = dnekclock_sync()
      if( mod(nid, np/2) .eq. np/2-2) then
         print *, 'set_proc_map ', endtime-starttime1
      endif
c     if(nid.eq.0) print *, "mapelpr in set_proc", endtime-starttime 
c     starttime = dnekclock_sync()
c
      DO 1200 IFIELD=MFIELD,NFLDT
         IF (IFTMSH(IFIELD)) THEN
            NELG(IFIELD)      = NELGT
         ELSE
            NELG(IFIELD)      = NELGV
         ENDIF
 1200 CONTINUE

C     Output the processor-element map:
      starttime1 = dnekclock_sync()
      ifverbm=.true.
      if (np.gt.2050.or.nelgt.gt.40000) ifverbm=.false.

      if(ifverbm) then
        idum = 1
        if(nid.eq.0) then
           N8 = min(8,nelt)
           write(6 ,1310) node-1,(lglel(ie),ie=1,n8)
           if (NELT.GT.8) write(6 ,1315) (lglel(ie),ie=9,NELT)
           DO inid=1,NP-1
              mtype = inid
              call csend(mtype,idum,4,inid,0)            ! handshake
              call crecv(mtype,inelt,4)               ! nelt of other cpus
              N8 = min(8,inelt)
c             write(6 ,1310) inid+1,(lglel(ie,inid+1),ie=1,n8)
c             IF (inelt.gt.8) 
c    &           write(6 ,1315) (lglel(ie,inid+1),ie=9,inelt)
           ENDDO
 1310      FORMAT(' RANK',I6,' IEG',8I8)
 1315      FORMAT('     ',6X,'    ',8I8)
        else
           mtype = nid
           call crecv(mtype,idum,4)                ! hand-shake
           call csend(mtype,nelt,4,0,0)            ! nelt
        endif
      endif
      endtime = dnekclock_sync()
      if( mod(nid, np/2) .eq. np/2-2) then
         print *, 'mapelpr_output ', endtime-starttime1
      endif

C     Check elemental distribution
C
C      IF (IPASS.EQ.2.AND.PARAM(156).eq.9) THEN
C         NXYZ=NX1*NY1*NZ1
C         DO 1400 IE=1,NELT
C            VTMP1=NODE
c            VTMP2=IE
C            CALL CFILL(VX(1,1,1,IE) ,VTMP1,NXYZ)
C            CALL CFILL(VY(1,1,1,IE) ,VTMP2,NXYZ)
C            CALL CFILL(T(1,1,1,IE,1),VTMP1,NXYZ)
C 1400    CONTINUE
C         call prepost(.true.,'   ')
C      ENDIF

      nn = iglmin(nelt,1)
      nm = iglmax(nelt,1)
      if(nio.eq.0) write(6,*) 'element load imbalance: ',nm-nn,nn,nm

      if(nio.eq.0) then
        write(6,*) 'done :: mapping elements to processors'
        write(6,*) ' '
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_proc_map()
C
C     Compute element to processor distribution according to (weighted) 
C     physical distribution in an attempt to minimize exposed number of
C     element interfaces.
C
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'SCRCT'
      include 'TSTEP'
      include 'ZPER'
      include 'CMTPART' 
      include 'mpif.h'
      common /ctmp0/ iwork(lelt)
      REAL*8 dnekclock,t0
 
c     added by keke
      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload
      !common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c      parameter (lr=16*ldim,li=5+6)
c     parameter (lr=76,li=10)
c     common  /cpartr/ rpart(lr,llpart) ! Minimal value of lr = 14*ndim+1
c     common  /cparti/ ipart(li,llpart) ! Minimal value of lr = 14*ndim+1
c     common  /iparti/ n,nr,ni
c     common /ptpointers/ jrc,jpt,je0,jps,jpid1,jpid2,jpid3,jpnn,jai
c    >                ,nai,jr,jd,jx,jy,jz,jx1,jx2,jx3,jv0,jv1,jv2,jv3
c    >                ,ju0,ju1,ju2,ju3,jf0,jar,jaa,jab,jac,jad,nar,jpid
c     common /ptpointers/ jrc,jpt,je0,jps,jpid1,jpid2,jpid3,jpnn,jpid
c    >                   ,jai,nai,    jr,jd,jx,jy,jz,jv0,ju0,jf0,jfusr
c    >                   ,jfqs,jfun,jfiu,jtaup,jcd,jdrhodt,jre,jDuDt
c    >                   ,jtemp,jrho,jrhop,ja,jvol,jdp,jar,jx1,jx2,jx3
c    >                   ,jv1,jv2,jv3,ju1,ju2,ju3,nar,jvol1,jgam
      real starttime1, endtime
      integer iwork2(lelg)
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl
      integer ip, e
      logical partl         ! This is a dummy placeholder, used in cr()
      nl = 0                ! No logicals exchanged
c     end added by keke

      t0 = dnekclock()
c     if (.not.(ifgtp.or.ifgfdm)) then
c     starttime = dnekclock_sync()
      if (.not.ifgtp) then
c
c        rsb element to processor mapping 
c
         if (ifgfdm)       call gfdm_elm_to_proc(gllnid,np) ! gfdm w/ .map
         starttime1 = dnekclock_sync() 
         call get_map
         endtime = dnekclock_sync()
         if( mod(nid, np/2) .eq. np/2-2) then
         print *, 'get_map ', endtime-starttime1
         endif

      endif

      if(ifzper.or.ifgtp) call gfdm_elm_to_proc(gllnid,np) ! special processor map

c     compute global to local map (no processor info)
c
      starttime1 = dnekclock_sync() 
      IEL=0
      CALL IZERO(GLLEL,NELGT)
      DO IEG=1,NELGT
         IF (GLLNID(IEG).EQ.NID) THEN
            IEL = IEL + 1
            GLLEL(IEG)=IEL
            NELT = IEL
            if (ieg.le.nelgv) NELV = IEL
         ENDIF
c        write(6,*) 'map2 ieg:',ieg,nelv,nelt,nelgv,nelgt
      ENDDO
c
c     dist. global to local map to all processors
c
c     endtime = dnekclock_sync()
c     if(nid.eq.0) print *, "set_proc 1", endtime-starttime 

      ! commented out by keke
c      npass = 1 + nelgt/lelt
c      k=1
c      do ipass = 1,npass
c         m = nelgt - k + 1
c         m = min(m,lelt)
c         if (m.gt.0) call igop(gllel(k),iwork,'+  ',m)
c         !if (m.gt.0) then
c         !   call mpi_allreduce (gllel(k),iwork,m,mpi_integer,
c     >   !          mpi_sum ,nekcomm,ierr)
c         !   call icopy(gllel(k),iwork,m)
c         !endif 
cc        endtime = dnekclock_sync()
cc        if(nid.eq.0) print *, "set_proc 2", endtime-starttime, ipass, 
cc    >      istep  
c         k = k+m
c      enddo

            call mpi_allreduce (gllel(1),iwork2,nelgt,mpi_integer,
     >             mpi_sum ,nekcomm,ierr)
            call icopy(gllel(1),iwork2,nelgt)
      endtime = dnekclock_sync()
      if( mod(nid, np/2) .eq. np/2-2) then
      print *, 'gllel_comp ', endtime-starttime1
      endif
c
c     compute local to global map
c     (i.e. returns global element number given local index and proc id)
c
c     starttime = dnekclock_sync()
      starttime1 = dnekclock_sync() 
      do ieg=1,nelgt
         mid  =gllnid(ieg)
         ie   =gllel (ieg)
         if (mid.eq.nid) lglel(ie)=ieg
      enddo
      endtime = dnekclock_sync()
      if( mod(nid, np/2) .eq. np/2-2) then
      print *, 'lglel_comp ', endtime-starttime1
      endif

c     added by keke to convert ipart(je0, i) to local element id
      starttime1 = dnekclock_sync() 
      if ( gfirst .eq. 0) then
         ip=0
         do ip = 1, n
            e = ipart(je0, ip)
            ipart(je0, ip) = gllel(e) - 1  ! je0 start from 0
         enddo
c     Sort by element number 
          call crystal_tuple_sort(i_cr_hndl,n
     $              , ipart,ni,partl,nl,rpart,nr,je0,1)
      endif
      endtime = dnekclock_sync()
      if( mod(nid, np/2) .eq. np/2-2) then
      print *, 'part_sort ', endtime-starttime1
      endif

c     endtime = dnekclock_sync()
c     if(nid.eq.0) print *, "set_proc 3", endtime-starttime 
c     starttime = dnekclock_sync()
c
c     All Done.
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_elm_to_proc(gllnid,np)
c
c
      include 'SIZE'
      include 'ZPER'
c
      integer gllnid(1)
c
      common /ctmp1/  map_st(lelg_sm)
      common /vptsol/ iwork(0:lp)
      integer nelbox(3),nstride_box(3)
c
      call gfdm_set_pst(ip,is,it,nelbox,nstride_box,nx2,ny2,nz2)
c
      nep = nelbox(ip)
      nes = nelbox(is)
     
      if(nelbox(it).eq.0) nelbox(it)=1
      net = nelbox(it)
c
      nst = nes*net
      if (nst.lt.np) then
         if (nid.eq.0) 
     $   write(6,*) 'ERROR, number of elements in plane must be > np'
     $   ,nst,np,nep,nes,net
         call exitt
      endif
c
c
      call gfdm_map_2d(map_st,nes,net,iwork,np)
      call gfdm_build_global_el_map (gllnid,map_st,nes,net
     $                                     ,nelbox,nstride_box,ip,is,it)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_map_2d(map_st,nes,net,num_el,np)
c
c     Set up a 2D processor decomposition of an NES x NET array.
c
      integer map_st(nes,net),num_el(0:np)
c
c     First a stupid dealing algorithm to determine how many
c     elements on each processor

      do i=0,np
         num_el(i) = 0
      enddo
c
      k = np-1
      do i=1,nes*net
         num_el(k) = num_el(k)+1
         k=k-1
         if (k.lt.0) k = np-1
      enddo
c
      jnid = 0
      nel_cnt = 0
      nel_cur = num_el(jnid)
      do j=1,net,2
         do i=1,nes                 ! Count down
            nel_cnt = nel_cnt + 1
            if (nel_cnt.gt.nel_cur) then
               jnid=jnid+1
               nel_cur = num_el(jnid)
               nel_cnt = 1
            endif
            map_st(i,j) = jnid
         enddo
c
         j1 = j+1
         if (j1.le.net) then
            do i=nes,1,-1                ! Count up
               nel_cnt = nel_cnt + 1
               if (nel_cnt.gt.nel_cur) then
                  jnid=jnid+1
                  nel_cur = num_el(jnid)
                  nel_cnt = 1
               endif
               map_st(i,j1) = jnid
            enddo
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_set_pst(ip,is,it,nelbox,nstride_box,nxp,nyp,nzp)
c
      include 'SIZE'
      include 'INPUT'
      include 'ZPER'
c
      integer nelbox(3),nstride_box(3)
c
      if (if3d) then
         if (param(118).lt.0) then
            ip = 3
            is = 1
            it = 2
         elseif (param(117).lt.0) then
            ip = 2
            is = 3
            it = 1
         else
            ip = 1
            is = 2
            it = 3
         endif
      else
         if (param(117).lt.0) then
            ip = 2
            is = 1
            it = 3
         else
            ip = 1
            is = 2
            it = 3
         endif
      endif
c
      pst2lex(1)=ip       ! identify x-, y- or z with primary direction
      pst2lex(2)=is
      pst2lex(3)=it
c
      lex2pst(ip)=1
      lex2pst(is)=2
      lex2pst(it)=3
c
      nelbox(1) = nelx
      nelbox(2) = nely
      nelbox(3) = nelz
c
      nstride_box(1) = 1
      nstride_box(2) = nelx
      nstride_box(3) = nelx*nely
c
      ngfdm_p(1) = nelx*nxp
      ngfdm_p(2) = nely*nyp
      ngfdm_p(3) = nelz*nzp
      write(6,*) 'ngfdm:',(ngfdm_p(k),k=1,3)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_build_global_el_map (gllnid,map_st,nes,net
     $                                     ,nelbox,nstride_box,ip,is,it)
c
      include 'SIZE'
      integer gllnid(1)
      integer map_st(nes,net)
      integer nelbox(3),nstride_box(3)
c
      integer proc
c
      do jt=1,nelbox(it)
      do js=1,nelbox(is)
         proc = map_st(js,jt)
         do jp=1,nelbox(ip)
            ieg = 1 + nstride_box(ip)*(jp-1)  ! nstride_p=nes*net
     $              + nstride_box(is)*(js-1)
     $              + nstride_box(it)*(jt-1)
            gllnid(ieg) = proc
         enddo
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outmati(u,m,n,name6)
      integer u(m,n)
      character*6 name6
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
c     Print out copies of a global matrix
c
      do mid=0,np-1
        call nekgsync
        if (mid.eq.nid) then
         n20 = min(n,20)
         write(6,1) nid,m,n,name6
   1     format(//,3i6,'  Matrix:',2x,a6,/)
         do i=1,m
            write(6,2) nid,name6,(u(i,j),j=1,n20)
         enddo
   2     format(i3,1x,a6,20i6)
        endif
        call nekgsync
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine outmati8(u,m,n,name6)
      integer*8 u(m,n)
      character*6 name6
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
c     Print out copies of a global matrix
c
      do mid=0,np-1
        call nekgsync
        if (mid.eq.nid) then
         n20 = min(n,20)
         write(6,1) nid,m,n,name6
   1     format(//,3i6,'  Matrix:',2x,a6,/)
         do i=1,m
            write(6,2) nid,name6,(u(i,j),j=1,n20)
         enddo
   2     format(i3,1x,a6,20i6)
        endif
        call nekgsync
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine get_map

      call get_vert

      return
      end
c-----------------------------------------------------------------------
