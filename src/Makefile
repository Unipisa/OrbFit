all: depend
	(cd suit; make)
	(cd prelim; make)
	(cd propag; make)
	(cd moid; make)
	(cd bineph; make)
	(cd fitobs; make)
	(cd orbfit; make)

nondistribution:
#	-(cd twoobs; make depend; make)
	-(cd allobs; make depend; make)
	-(cd idconf2; make depend; make)
	-(cd stat; make depend; make)
	-(cd orb9; make depend; make)
	-(cd astdys; make depend; make)
	-(cd neodys; make depend; make)
	-(cd newret; make depend; make)
#	(cd skymap; make depend; make)
#	(cd starcat; make)

nondistclean:
#	-(cd twoobs; make clean)
	-(cd allobs; make clean)
	-(cd idconf2; make clean)
	-(cd stat; make clean)
	-(cd orb9; make clean)
	-(cd astdys; make clean)
	-(cd neodys; make clean)
	-(cd newret; make clean)
#	-(cd skymap; make clean)
#	-(cd starcat; make clean)


clean:
	(cd suit; make clean)
	(cd prelim; make clean)
	(cd propag; make clean)
	(cd moid; make clean)
	(cd bineph; make clean)
	(cd fitobs; make clean)
	(cd orbfit; make clean)
	(cd jpleph; make clean)
	rm -f *~

distclean:
	(cd suit; make distclean)
	(cd prelim; make distclean)
	(cd moid; make clean)
	(cd propag; make distclean)
	(cd bineph; make distclean)
	(cd fitobs; make distclean)
	(cd orbfit; make distclean)
	(cd jpleph; make distclean)
	rm -f lib/lib*.a
	cp make.flags.dist make.flags
	rm -f *~

depend:
	(cd suit; make depend)
	(cd prelim; make depend)
	(cd propag; make depend)
	(cd bineph; make depend)
	(cd fitobs; make depend)
	(cd orbfit; make depend)

win:
	@copy dummy.dep suit\make.dep > junk
	@copy dummy.dep prelim\make.dep > junk
	@copy dummy.dep propag\make.dep > junk
	@copy dummy.dep moid\make.dep > junk
	@copy dummy.dep bineph\make.dep > junk
	@copy dummy.dep fitobs\make.dep > junk
	@copy dummy.dep orbfit\make.dep > junk
	cd suit
	nmake /nologo win
	cd ..\prelim
	nmake /nologo win
	cd ..\propag
	nmake /nologo win
	cd ..\moid
	nmake /nologo win
	cd ..\bineph
	nmake /nologo win
	cd ..\fitobs
	nmake /nologo win
	cd ..\orbfit
	nmake /nologo win
	cd ..\jpleph
	nmake /nologo win

winclean:
	@cd suit
	@nmake /nologo winclean
	@cd ..\prelim
	@nmake /nologo winclean
	@cd ..\propag
	@nmake /nologo winclean
	@cd ..\moid
	@nmake /nologo winclean
	@cd ..\bineph
	@nmake /nologo winclean
	@cd ..\fitobs
	@nmake /nologo winclean
	@cd ..\orbfit
	@nmake /nologo winclean
	@cd ..\jpleph
	@nmake /nologo winclean

windist:
	@cd suit
	@nmake /nologo winclean
	@cd ..\prelim
	@nmake /nologo winclean
	@cd ..\propag
	@nmake /nologo winclean
	@cd ..\moid
	@nmake /nologo winclean
	@cd ..\bineph
	@nmake /nologo windist
	@cd ..\fitobs
	@nmake /nologo windist
	@cd ..\orbfit
	@nmake /nologo windist
	@cd ..\orbfit
	@nmake /nologo windist







