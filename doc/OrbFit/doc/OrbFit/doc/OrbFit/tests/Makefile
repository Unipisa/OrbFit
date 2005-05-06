all: fitobs orbfit bineph
	(cd fitobs; make)
	(cd bineph; make)
	(cd orbfit; make)
	(cd moid; make)

clean: distclean

distclean:
	(cd bineph; make distclean)
	(cd orbfit; make distclean)
	(cd fitobs; make distclean)
	(cd moid; make distclean)  