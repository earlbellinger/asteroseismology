      subroutine setusr(user)
c
c  returns user name of user running the programme.
c  Note: this is based on C38 user IDs. List must be updated
c  when new relevant users are added
c
c  Original version: 3/12/92
c
      character*(*) user
      integer getuid
c
      id=getuid()
c
      if(id.eq.534) then
        user='aake'
      else if(id.eq.535) then
        user='brandenb'
      else if(id.eq.503) then
        user='bt'
      else if(id.eq.519) then
        user='fgj'
      else if(id.eq.509) then
        user='hans'
      else if(id.eq.501) then
        user='jcd'
      else if(id.eq.530) then
        user='jens'
      else if(id.eq.515) then
        user='jesm'
      else if(id.eq.554) then
        user='jones'
      else if(id.eq.542) then
        user='kj'
      else if(id.eq.544) then
        user='michel'
      else if(id.eq.505) then
        user='mjt'
      else if(id.eq.514) then
        user='mlo'
      else if(id.eq.537) then
        user='ms'
      else if(id.eq.547) then
        user='mv'
      else if(id.eq.541) then
        user='nha'
      else if(id.eq.556) then
        user='nuspl'
      else if(id.eq.513) then
        user='pen'
      else if(id.eq.512) then
        user='pg'
      else if(id.eq.543) then
        user='pot'
      else if(id.eq.510) then
        user='schou'
      else if(id.eq.500) then
        user='srf'
      else if(id.eq.529) then
        user='tang'
      else if(id.eq.536) then
        user='th'
      else if(id.eq.581) then
        user='sushant'
      else if(id.eq.502) then
        user='oz'
      else
	write(6,100) id
	stop
      endif
      return
  100 format(//' ***** Error in setusr. ID = ',i5,' not in table')
      end
