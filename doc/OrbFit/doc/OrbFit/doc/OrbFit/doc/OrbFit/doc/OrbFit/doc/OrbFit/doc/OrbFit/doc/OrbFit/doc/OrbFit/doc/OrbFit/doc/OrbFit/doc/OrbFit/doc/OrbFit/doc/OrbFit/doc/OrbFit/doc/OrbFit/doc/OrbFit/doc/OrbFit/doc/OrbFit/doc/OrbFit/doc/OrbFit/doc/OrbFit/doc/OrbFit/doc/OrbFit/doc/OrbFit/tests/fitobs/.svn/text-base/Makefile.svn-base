all: test_fitobs test_catpro

test_fitobs : fitobs.x 1220T-2.fel 2003BH84.fel 2004XK3.fel

1220T-2.fel : fitobs.x
	time fitobs.x < 1220T-2.inp

2003BH84.fel: fitobs.x
	time fitobs.x < 2003BH84.inp

2004XK3.fel: fitobs.x
	time fitobs.x < 2004XK3.inp

test_catpro: 
	time catpro.x

distclean:
	-(rm gmon.out core *.pro *.clo *.fga *.rwo *.fou *.fel *.cbd \
	giffv.* *.frw *.*~ *~ *.tmp *.err *.done multprop.ctc *.mou \
	virimp.* *.vis mult.ctc)


