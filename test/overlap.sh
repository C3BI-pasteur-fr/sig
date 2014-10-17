#! /bin/sh

## Verbose mode
if test "x$VERBOSE" = xx; then
  set -x
fi

## Call program

NUM=`../src/sig -i -p '[RK]-X-V-X-[FW] (0,) F-X-X-[RK]-X-[RK]' $srcdir/seq-test | wc -l`

exit `expr $NUM - 3`
