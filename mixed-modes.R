DF <- read.table('profile170-freqs.dat', col.names=c('l', 'n', 'nu', 'inertia'))
#for ( l_deg in 0:4 ) {
#    if ( any( diff( DF[DF[['l']] == l_deg, ][['inertia']] ) < 0 ) ) {
#        invisible(cat(1))
#        quit()
#    }
#}
#invisible(cat(0))

attach(DF)
invisible(cat(ifelse( # check that all of the nonradial modes aren't bumped 
    length(n==0) == 4
 && any( inertia[ l==1 & n==min( n[ l==1 & n>0 ] ) ] > inertia[ l==1 & n==0 ] )
 && any( inertia[ l==2 & n==min( n[ l==2 & n>0 ] ) ] > inertia[ l==2 & n==0 ] )
 && any( inertia[ l==3 & n==min( n[ l==3 & n>0 ] ) ] > inertia[ l==3 & n==0 ] ),
    1, 0 )))

