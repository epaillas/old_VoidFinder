all: SVF CVF

SVF:
	make -C src -f Makefile1 install

CVF:
	make -c src -f Makefile2 install

clean:
	make -C src clean
	rm -f bin/*.exe
