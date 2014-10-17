SUFFIXES = .pod .1

man_MANS = sig.1

EXTRA_DIST = $(man_MANS) sig.pod

PODARGS = --center='User Manuals' --release='Unix'

.pod.1:
	$(POD2MAN) $(PODARGS) $< > $@ && touch $@

parano:
	(for p in `ls *.pod`; do $(PODCHECKER) $$p; done)
