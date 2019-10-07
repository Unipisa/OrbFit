all: test_fitobs

test_fitobs : fitobs.x 
	time fitobs.x < 1220T-2.inp

distclean:
	-(rm gmon.out core *.pro *.clo *.fga *.rwo *.fou *.fel *.cbd \
	giffv.* *.frw *.*~ *~ *.tmp; cd obsdata; rm *.frw *.rwo mult.ctc; \
	cd xephem; rm *.edb)


