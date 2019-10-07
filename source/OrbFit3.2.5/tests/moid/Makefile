
all: test_moid

test_moid : min_comp.x 
	time min_comp.x < test_moid.inp

distclean: clean

clean:
	(rm -f *.minpts CHECK.* gmon.out mincomp.log mincomp.names.fail *~ \#*\#)