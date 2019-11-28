all: SVF CVF

SVF:
	make -C src -f Makefile1 install

CVF:
	make -C src -f Makefile2 install

clean:
	make -C src -f Makefile1 clean
	make -C src -f Makefile2 clean
	rm -f bin/*.exe
