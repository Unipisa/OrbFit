#!	/bin/sh

docdir=`( cd ../../doc ; pwd )`
echo "default documentation directory is $docdir"
cat > doclib.h << END
* Copyright (C) 1998 by Andrea Milani (milani@dm.unipi.it)
* Version: August 4, 1998
* ---------------------------------------------------------------------
* Default library directory
      CHARACTER*100 ddocd
      PARAMETER (ddocd=
     +   '${docdir}')
END
