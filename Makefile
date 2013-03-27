#
# makefile for vdm.
#
# fec.S.980624a is an optimized version for use on 486 and old pentium
# machines. It is only for GF_BITS=8 and generally does not work
# fast on systems with multiple instruction pipelines (PentiumPro,
# Pentium2)
# same for fec.S16.980624a , for use with GF_BITS=16
#
# gcc does something strange, so check the various opt. levels for
# best performance (or write addmul1 in assembly code).
#
# Standard compilation with -O9 works well for PentiumPro and Pentium2
# machines.
#

COPT= -O3 -funroll-loops -DGF_BITS=8
CFLAGS=$(COPT) -Wall # -DTEST

rs-codec.a: rs-codec.o rs-codec.h
	echo $LD -o $< $@
