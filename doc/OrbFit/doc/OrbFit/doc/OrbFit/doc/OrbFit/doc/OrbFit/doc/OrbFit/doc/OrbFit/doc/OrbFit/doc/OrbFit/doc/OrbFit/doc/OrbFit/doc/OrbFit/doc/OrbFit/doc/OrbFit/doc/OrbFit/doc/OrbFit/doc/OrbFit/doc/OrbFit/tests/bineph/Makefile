all: CPV.bep

distclean: 
	rm -f *.bep *.bai core gmon.out *~ *.*~

CPV.bep: CPV.bop CPV.inc ../../src/bineph/bineph.x
	echo CPV | bineph.x

../../src/bineph/bineph.x:
	(cd ../..; make)
