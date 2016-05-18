kentSrc = /mnt/lab_data/kundaje/projects/encodeenhancerpredict/igr/kentUtils/src
L += -lm -lz -lssl
CFLAGS += -I$(kentSrc)/inc

all:
	gcc -shared -fPIC bigWigFeaturize.c -o bigWigFeaturize.so \
	-I/users/nasa/nandi/include/python2.7 \
	-I ~/nandi/lib/python2.7/site-packages/numpy/core/include/  \
	-I/mnt/lab_data/kundaje/projects/encodeenhancerpredict/igr/kentUtils/src/inc \
	/mnt/lab_data/kundaje/projects/encodeenhancerpredict/igr/kentutils/src/lib/local/jkweb.a \
	/mnt/lab_data/kundaje/projects/encodeenhancerpredict/igr/kentutils/samtabix/lib*.a \
	-lz -lssl -lm -g

test:
	python test.py
