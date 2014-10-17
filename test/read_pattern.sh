#! /bin/sh

## Verbose mode
if test "x$VERBOSE" = xx; then
  set -x
fi

# read pattern from file

# reading syntax one
ret=`../src/sig -f $srcdir/pattern $srcdir/seq-test  | wc -l|| exit 1`
if [ $ret != 2 ]
	then 
		exit 1
fi



# reading syntax two
ret=`../src/sig -f $srcdir/pattern2 $srcdir/seq-test  | wc -l|| exit 1`
if [ $ret != 2 ]
        then
	                exit 1
fi


