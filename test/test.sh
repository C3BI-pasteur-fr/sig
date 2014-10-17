#! /bin/sh

## Verbose mode
if test "x$VERBOSE" = xx; then
  set -x
fi

## Call program
../src/sig -p '[RK]-X-V-X-[FW] (0,) F-X-X-[RK]-X-[RK]' $srcdir/seq-test > /dev/null
