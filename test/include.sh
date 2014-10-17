#! /bin/sh

## Verbose mode
if test "x$VERBOSE" = xx; then
  set -x
fi

## Call program

NUM=`../src/sig -p "[RK]-X(0,1)-V-X-[FW] (0,) F-x(0,1)-V-x-F" $srcdir/included | wc -l`

exit `expr $NUM - 4` 
