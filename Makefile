all: .conf.status
	(cd src; make)

tests: test

test:
	(cd tests; make)

clean:
	(cd src; make clean)
	(cd lib; make clean)

del_obsolete_files:
	rm -f `cat dellist`

distclean:
	rm -f *~
	(cd src; make distclean)
	(cd lib; make distclean)
	(cd tests; make distclean)
	-(cd src/lib; rm *.a)
	(rm -f .conf.status)

.conf.status:
	@ echo "Please run \"config\" before doing \"make\""
	@ exit 1

distribution: distclean notar
	tar -cvf ../OrbFit3.4.2.tar -X notar . ; \
	gzip ../OrbFit3.4.2.tar

doctar:
	tar -cvf ../doc.tar --exclude-from notar ./doc; gzip ../doc.tar

additional_doc:
	tar -cvf ../additional_doc.tar ./doc/additional_doc; gzip ../additional_doc.tar 

nondistribute: 
	cd src; make nondistclean; cd .. 
	tar -T notar -cvf ../OrbFitwork342.tar ; gzip ../OrbFitwork342.tar

panst: 
	cd src/panst; make clean; cd ../../tests/panst; make distclean
	tar -czf ../panst3.4.2.tgz src/panst tests/panst lib/orbsrv.key

starcat:
	cd ../starcat; make clean; cd ../skymap; make clean; cd ../..;
	tar -cvzf ../starcat.tgz src/starcat src/skymap 

nondistclean:
	cd src; make nondistclean

patch: 
	tar -T patchlist -cvf ../patch3.4.2.tar ; gzip ../patch3.4.2.tar

