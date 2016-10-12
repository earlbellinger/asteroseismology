      block data blstio
c
c  define unit numbers for standard input and output
c
c  istdin: standard input (typically terminal)
c  istdou: standard output for diagnostics, etc (typically terminal)
c  istdpr: printed output (terminal or printer)
c
c
c  hp9000 version
c  ************
c
      common/cstdio/ istdin, istdou, istdpr
c
      data istdin, istdou, istdpr
     *  /    5,       6,     6    /
      end
