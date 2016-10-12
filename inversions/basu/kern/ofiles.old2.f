      subroutine ofiles
c
c  reads unit numbers and file names from standard input, in
c  format
c
c  <unit number>   <file name>
c
c  if file name is given as 0, /dev/null is used for the file.
c  input ends with EOF or a line containing -1.
c  returns number of files in nfiles, unit numbers in ids(.),
c  and file names in file(.).
c
c  19/8/87: modified for HP9000 by taking out action option
c     in s/r openf.
c
c  21/8/87: modified for HP9000 by replacing multiple occurences
c     of /dev/null by scratch files, since the HP9000, unfortunately,
c     does not allow several unit numbers to be associated with
c     /dev/null
c  ..............................................
c
      character*80 file, filein
      common/cofile/ nfiles, ids(20), file(20)
c
c  the following line required on CR32
c
c..      save nfiles, ids, file
c
      write(6,*) 'Input format: <unit number>   <file name>'
      write(6,*) 'input ends with EOF or a line containing -1.'
      write(6,*)
     *  'if file name is given as 0, /dev/null is used for the file.'
c
      nfiles=0
   10 read(5,*, end=30) idsin, filein
      if(idsin.lt.0) go to 30
c
c  test for /dev/null
c
      if(filein.eq.'0') then
        filein='/dev/null'
      end if
c
      nfiles=nfiles+1
      ids(nfiles)=idsin
      file(nfiles)=filein
      go to 10
c
   30 continue
c
c  output file information
c
      write(6,100)
      do 40 n=1,nfiles
   40 write(6,110) ids(n),file(n)
      write(6,120)
      return
  100 format(/' files set in s/r ofiles:'/)
  110 format(i3,2x,'''',a,'''      @')
  120 format('-1      ''''          @')
      end
      subroutine stfile(idsst, nfst)
c
c  find number of file nfst corresponding to unit number idsst.
c  list of unit numbers and file names in ids and file must have been
c  set up in common/cofile/ by call of ofiles.
c
      character*80 file
      common/cofile/ nfiles, ids(20), file(20)
c
c  the following line required on CR32
c
c..      save nfiles, ids, file
c
      do 10 i=1,nfiles
      if(idsst.eq.ids(i)) then
        nfst=i
        go to 20
      end if
   10 continue
c
c  idsst not found, print diagnostics
c
      write(6,*) idsst,' not found'
      write(6,*) 'List of files available:'
      do 15 i=1,nfiles
   15 write(6,*) ids(i),'  ',file(i)
c
      nfst=-1
c
   20 continue
      return
      end
      subroutine openf(id,status,form)
c
c  open file with unit number id, status as in string status,
c  and format as in string form.
c  status and form may be abbreviated to a single character,
c  as i.e.'n' for 'new', 'u' for 'unformatted'.
c
c  s/r ofiles must have been called previously to set up
c  nfiles, ids and file in common /cofile/.
c
c  original version 30/9/86
c
c            ....................................
c
      character*(*) status, form
      character*80 stat1, form1, file, files
      character*24 timest,fdate
      character ss1*1, ss2*2
      common/cofile/ nfiles, ids(20), file(20)
c
      save 
      data nnul /0/
c
c  find file name
c
      call stfile(id,nfin)
      if(nfin.lt.0) go to 90
c
c  set full status
c
      if(status(1:1).eq.'o') then
        stat1='old'
      else if(status(1:1).eq.'n') then
        stat1='new'
      else if(status(1:1).eq.'s') then
        stat1='scratch'
      else 
        stat1='unknown'
      end if
c
c  for /dev/null, set status to old
c
      if(file(nfin).eq.'/dev/null') then
        nnul = nnul+1
	write(6,*) 'Now nnul =',nnul
        if(nnul.eq.1) then
          stat1='old'
        else
          stat1='scratch'
        end if
      end if
c
c  set full format
c
      if(form(1:1).eq.'u') then
        form1='unformatted'
      else
        form1='formatted'
      end if 
c
c  open file
c
      if(stat1.eq.'scratch') then
	write(ss1,'(i1)') nnul
        files='/scratch/jcd/scratch/ttt.'//fdate()//'.'//ss1
	write(6,'(2a)') 'Scratch file: ',files
        open(id,status='unknown',disp='delete',file=files,form=form1)
c
c  diagnostic output
c
        write(6,110) id,stat1,form1
c
      else
        open(id,file=file(nfin),status=stat1,form=form1)
c
c  diagnostic output
c
        write(6,100) id,file(nfin),stat1,form1
c
      end if
c
      return
c
c  error in locating file name. exit
c
   90 stop
  100 format(' open(',i3,',file=',a30,',status=',a10,',form=',a10,')')
  110 format(' open(',i3,',status=',a,',form=',a,')')
      end
      subroutine openfc(id,idp,status,form)
c
c  open file with unit number id, status as in string status,
c  format as in string form. For details, see s/r openf.
c
c  open only takes place if id .ne. idp. If idp .gt. 0,
c  also closes unit idp. idp is returned as id.
c
      character*(*) status, form
      if(id.ne.idp) then
        if(idp.gt.0) close(idp)
        call openf(id,status,form)
        idp=id
      end if
c
      return
      end
