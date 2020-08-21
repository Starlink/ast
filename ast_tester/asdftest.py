#!/usr/bin/python3
import math
import asdf
from gwcs import wcs

asdf_file = asdf.open("tanSipWcs.asdf")
wcsobj = asdf_file.tree['wcs']

ab_targets = ( (0.74087499500083731, 0.76582017043755046),
               (0.73971051729252979, 0.76580407766915781),
               (0.74040832140984847, 0.76461426513798736) ,
               (0.74111042497286139, 0.76572737656088874),
               (0.74200902465461160, 0.76703478086374388) ,
               (0.74167589835205705, 0.76703043529783410) )

xylist = ((-25.0,0.0), (-200.0,0.0), (-100.0,-250.0), (10.0,-20.0),
          (150.0,250.0), (100.0,250.0) )

ok = True
i = 0
for (x,y) in xylist:
   (a,b) = wcsobj( x, y )

   a_target = math.degrees(ab_targets[ i ][ 0 ])
   b_target = math.degrees(ab_targets[ i ][ 1 ])

   if abs( a - a_target ) > 1E-12:
      ok = False
      print("Point {0}: a={1} (expected {2})".format(i,a,a_target))
      break
   elif abs( b - b_target ) > 1E-12:
      ok = False
      print("Point {0}: b={1} (expected {2})".format(i,b,b_target))
      break

   (xi,yi) = wcsobj.invert( a, b )

   if abs( xi - x ) > 1E-8:
      ok = False
      print("Point {0}: xi={1} (expected {2})".format(i,xi,x))
      break
   elif abs( yi - y ) > 1E-8:
      ok = False
      print("Point {0}: yi={1} (expected {2})".format(i,yi,y))
      break

   i += 1

if ok:
   print( "asdftest: all test passed" )
else:
   print( "asdftest: tests failed" )
