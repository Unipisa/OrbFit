#!	/bin/sh

libdir=`( cd ../../lib ; pwd )`
echo "default library directory is $libdir"
cat > parlib.h << END
* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: September 2, 1996
* ---------------------------------------------------------------------
* Default library directory
      CHARACTER*100 dlibd
      PARAMETER (dlibd=
     +   '${libdir}')
END
