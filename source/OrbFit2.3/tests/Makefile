all: fitobs orbfit bineph
	(cd fitobs; make)
	(cd bineph; make)
#	(cd orbfit; make)

clean:
	(cd bineph; make distclean)
	(cd orbfit; make distclean)
	(cd fitobs; make distclean)

distclean:
	(cd bineph; make distclean)
	(cd orbfit; make distclean)
	(cd fitobs; make distclean)
