all: fitobs orbfit bineph
	(cd orbfit; make)
	(cd fitobs; make)
	(cd bineph; make)

clean:

distclean:
	(cd bineph; make distclean)
	(cd orbfit; make distclean)
	(cd fitobs; make distclean)
