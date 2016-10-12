      subroutine setscr(scrusr)
c
c  returns path for scratch files for user running program
c  in the form /scratch/user. This is specific to the Convexes
c
c  Original version: 3/12/92
c
      character*(*) scrusr
      character*80 user
c
      call setusr(user)
      scrusr='/scratch/'//user
      return
      end
