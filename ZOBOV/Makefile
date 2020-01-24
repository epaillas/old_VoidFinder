
all: qhull voboz

.PHONY : qhull voboz examples

qhull:
	make -C qhull/src

voboz: 
	make -C src install

examples:
	make -C examples
	make -C examples data

clean:
	make -C src clean
	make -C qhull/src cleanall
	make -C examples clean
	rm -f bin/vozinit bin/voz1b1* bin/voztie bin/jovoz bin/boz bin/jozov