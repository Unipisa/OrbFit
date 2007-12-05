all: fitobs orbfit bineph
	(cd bineph; make)
	(cd fitobs; make)
	(cd orbfit; make)
	(cd moid; make)

clean: distclean

distclean:
	(cd bineph; make distclean)
	(cd orbfit; make distclean)
	(cd fitobs; make distclean)
	(cd moid; make distclean)  