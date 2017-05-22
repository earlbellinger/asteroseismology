      subroutine ofiles(nfiles, ids, file)
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
c  ..............................................
c
      character*80 file, filein
      dimension ids(1),file(1)
      write(6,*) 'Input format: <unit number>   <file name>'
      write(6,*) 'input ends with EOF or a line containing -1.'
      write(6,*)
     *  'if file name is given as 0, /dev/null is used for the file.'
c
      nfiles=0
   10 read(9,*, end=30) idsin, filein
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
   40 write(6,110) ids(n),trim(file(n))
      return
  100 format(/' files set in s/r ofiles:'/)
  110 format(i3,2x,a)
      end
      subroutine stfile(nfiles, ids, file, idsst, nfst)
c
c  given list of unit numbers and file names in ids and file,
c  find number of file nfst corresponding to unit number idsst.
c
      character*80 file
      dimension ids(1),file(1)
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
   15 write(6,*) ids(i),'  ',trim(file(i))
c
      nfst=-1
c
   20 continue
      return
      end
      subroutine openf(id,status,form,nfiles,ids,file)
c
c  open file with unit number id, status as in string status,
c  and format as in string form.
c  status and form may be abbreviated to a single character,
c  as i.e.'n' for 'new', 'u' for 'unformatted'.
c
c  s/r ofiles must have been called previously to set up
c  nfiles, ids and file.
c
c  original version 30/9/86
c
c            ....................................
c
      character*(*) status, form, file
      character*80 stat1, form1
      dimension ids(1),file(1)
c
c  find file name
c
      call stfile(nfiles,ids,file,id,nfin)
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
      open(id,file=file(nfin),status=stat1,form=form1)
c
c  diagnostic output
c
      write(6,100) id,trim(file(nfin)),trim(stat1),trim(form1)
      return
c
c  error in locating file name. exit
c
   90 stop
  100 format(' open(',i3,',file=',a,',status=',a,',form=',a,')')
      end
