* Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: January 13, 1998
* ---------------------------------------------------------------------
* LIST OF RECOGNIZED KEYWORDS
*
* nkls        -  Number of entries in the list
* ns2it       -  Total number of string to integer translations 
* ns2i(i)     -  Number of string to integer translations 
*                for each keyword (i=1,nkls)
* ipos2i(i)   -  Position of the zero-th integer translation descriptors
*                for each keyword (i=1,nkls)
* intlst(k)   -  Integer value of the keyword (k=1,ns2it)
* keytyp      -  Keyword type:
*                      1 - integer
*                      2 - real
*                      3 - character string
*                      4 - logical
*                      5 - MJD
*                      6 - REF (reference system: see rdnref.f)
* iickls      -  Initialization check
* keylst      -  List of keywords
* vallst(k)   -  String value of the keyword (k=1,ns2it)
*
      INTEGER nkls,ns2it,ns2i,ipos2i,intlst,keytyp,iickls
      CHARACTER*(lckx) keylst
      CHARACTER*(lcvx) vallst

      COMMON/ckls1/nkls,ns2it,ns2i(nklsx),ipos2i(nklsx),
     +             intlst(ns2itx),keytyp(nklsx),iickls
      COMMON/ckls2/keylst(nklsx)
      COMMON/ckls3/vallst(ns2itx)
