bin_PROGRAMS = sig

sig_SOURCES = main.c error.c scan_motif.c seq-reader.c new_parse.c

noinst_HEADERS = params.h error.h scan_motif.h seq-reader.h new_parse.h

## Maintainer parano check
LINTDEFS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS)

parano:
	$(LINT) $(LINTFLAGS) $(LINTDEFS) $(sig_SOURCES)
