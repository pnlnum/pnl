GCC=gcc-4.8
GXX=g++-4.8

all: build build-gcc build-gcc-O2 build-O2

build:
	(cd build; \
	cmake -DCMAKE_BUILD_TYPE=Debug ..; \
	make -j 4; \
	make install)

build-O2:
	(cd build-O2; \
	cmake -DCMAKE_BUILD_TYPE=Release ..; \
	make -j 4; \
	make install)

build-gcc:
	(cd build-gcc; \
	cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=$(GCC) -DCMAKE_CXX_COMPILER=$(GXX) ..; \
	make -j 4; \
	make install)

build-gcc-O2:
	(cd build-gcc-O2; \
	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=$(GCC) -DCMAKE_CXX_COMPILER=$(GXX) ..; \
	make -j 4; \
	make install)


clean:
	rm -rf build/*
	rm -rf build-O2/*
	rm -rf build-gcc/*
	rm -rf build-gcc-O2/*

.PHONY: build build-gcc build-O2 build-gcc-O2
	
# vim: ft=Make:
