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
	tar -cvf ../OrbFit2.3.tar -X notar . ; \
	gzip ../OrbFit2.3.tar

doctar:
	tar -cvf ../doc.tar --exclude-from notar ./doc; gzip ../doc.tar

additional_doc:
	tar -cvf ../additional_doc.tar ./doc/additional_doc; gzip ../additional_doc.tar 

nondistribute: 
	cd src; make nondistclean; cd .. 
	tar -T notar -cvf ../OrbFitwork23.tar ; gzip ../OrbFitwork23.tar

starcat:
	cd ../starcat; make clean; cd ../skymap; make clean; cd ../..;
	tar -cvzf ../starcat.tgz src/starcat src/skymap 

nondistclean:
	cd src; make nondistclean

orb9dist:
	cd src/orb9; make distclean; cd ../..; cd tests/orbit9; make clean; \
	cd ../..; tar -cf ../orbit9.tar src/orb9 tests/orbit9 doc/ORBIT9 \
	doc/ORBHELP doc/READMEorb ; gzip ../orbit9.tar

skymapdist:
	cd src/starcat;make clean; cd ../skymap; make clean; cd ../..;\
	tar -czf ../skymap.tgz --exclude src/starcat/cfitsio src/skymap src/starcat 

patch: 
	tar -T patchlist -cvf ../patch2.3.1.tar ; gzip ../patch2.3.1.tar

#Windows Targets
winstall:
	copy conf\make.flags.win src\make.flags
	copy conf\sysdep.win src\include\sysdep.h
	copy conf\doclib.win src\fitobs\doclib.h
	copy conf\parlib.win src\suit\parlib.h
	copy conf\fszer2.win src\suit\fszer2.f
	cd src
	nmake /nologo win

windist:
	@del lib\jpleph
	@cd tests\fitobs
	@call cleanfit.bat
	@cd ..\orbfit
	@call cleanorb.bat
	@cd ..\bineph
	@call cleanbep.bat
	@cd ..\..\src
	@nmake /nologo windist
