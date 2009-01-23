all: 1996RO13_gauss.oel 1996RO13_diffcor.oel 1996RO13_ephem1.oep \
     1996RO13_ephem2.oep 1996RO13_ephem3.oep 1996RO13_ident1.oel
#     1996RO13_ident2.oel 4067.oel 1998PB1.oel

OK:  1996RO13_gauss.oel 1996RO13_diffcor.oel 1996RO13_ephem1.oep \
     1996RO13_ephem2.oep 1996RO13_ephem3.oep 

distclean: clean

clean:
	rm -f *.pro *.clo *.odc *.oel *.olg *.oep *.rwo *~ reject.rep *.err

4067.oel: 4067.obs ../../src/orbfit/orbfit.x
	echo 4067 | ../../src/orbfit/orbfit.x

1998PB1.oel: 1998PB1.obs ../../src/orbfit/orbfit.x
	echo 1998PB1 | ../../src/orbfit/orbfit.x

1996RO13_gauss.oel: 1996RO13_gauss.oop ../../src/orbfit/orbfit.x
	echo 1996RO13_gauss | ../../src/orbfit/orbfit.x

1996RO13_diffcor.oel: 1996RO13_diffcor.oop 1996RO13_gauss.oel \
                      ../../src/orbfit/orbfit.x
	echo 1996RO13_diffcor | ../../src/orbfit/orbfit.x

1996RO13_ephem1.oep: 1996RO13_diffcor.oel ../../src/orbfit/orbfit.x
	echo 1996RO13_ephem1 | ../../src/orbfit/orbfit.x

1996RO13_ephem2.oep: astorb.sub ../../src/orbfit/orbfit.x
	echo 1996RO13_ephem2 | ../../src/orbfit/orbfit.x

1996RO13_ephem3.oep: astorb.sub ../../src/orbfit/orbfit.x
	echo 1996RO13_ephem3 | ../../src/orbfit/orbfit.x

1996RO13_ident1.oel: 1996RO13_ident1.oop astorb.sub \
                    ../../src/orbfit/orbfit.x
#	rm -f 1996RO13=1220T-2.rwo
	cp 1996RO13=1220T-2.rwo.sav 1996RO13=1220T-2.rwo
	echo 1996RO13_ident1 | ../../src/orbfit/orbfit.x
	mv 1996RO13=1220T-2.rwo 1996RO13=1220T-2.ident1.rwo

1996RO13_ident2.oel: 1996RO13_ident2.oop 1996RO13_ident1.oel \
                     ../../src/orbfit/orbfit.x 
#	discard 1996RO13=1220T-2.ident1.rwo 2.5 > 1996RO13=1220T-2.rwo
# 	Next line is new.
	cp 1996RO13=1220T-2.rwo.sav 1996RO13=1220T-2.rwo
	cp 1996RO13=1220T-2.rwo 1996RO13=1220T-2.in.rwo
	echo 1996RO13_ident2 | ../../src/orbfit/orbfit.x
	mv 1996RO13=1220T-2.rwo 1996RO13=1220T-2.ident2.rwo

../../src/orbfit/orbfit.x:
	(cd ../..; make)
