SUBDIRS = lib src
if ALPS_FULL
SUBDIRS += tool example
endif
SUBDIRS += test doc .
DIST_SUBDIRS = $(SUBDIRS)

distdir = alps-$(VERSION)
distdir_light = alps-light-$(VERSION)

dist_noinst_SCRIPTS = bootstrap config/make_snapshot config/run-test config/update_preamble
dist_noinst_DATA = ChangeLog ChangeLog.0 INSTALL TODO
if ALPS_FULL
dist_noinst_DATA += CITATIONS.txt LICENSE.txt README config/ac_alps.m4 config/preamble.in
endif
if ALPS_LIGHT
dist_noinst_DATA += LICENSE-light.txt README-light config/preamble-light.in
endif

nodist_data_DATA = config/config.site
nodist_pkgdata_DATA = config/include.mk

.PHONY : lib libs tool tools test tests example examples doc docs
lib libs :
	$(MAKE) -C lib all
tool tools : libs
	$(MAKE) -C tool all
test tests : libs
	$(MAKE) -C test check
example examples : libs
	$(MAKE) -C example check
doc docs : lib tool
	$(MAKE) -C doc doc

install-data-local : config/config.cache
	@$(mkinstalldirs) $(DESTDIR)$(pkgdatadir)
	sed "s/ac_cv_env/save_ac_cv_env/" config/config.cache | sed "s/\^\\\'\^/^/g" > $(pkgdatadir)/config.cache
uninstall-local :
	rm -f $(pkgdatadir)/config.cache

.PHONY : update-preamble
update-preamble:
	if test -x "$(top_srcdir)/config/update_preamble"; then \
	  find $(top_srcdir) -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs $(top_srcdir)/config/update_preamble; \
	fi

.PHONY : release snapshot boost-snapshot boost-sandbox-snapshot
release snapshot :
	rm -rf $(distdir) $(distdir_light)
	$(MAKE) dist
	tar zxf $(distdir).tar.gz
	mv $(distdir) $(distdir_light)
	cd $(distdir_light); \
	  rm -f CITATIONS.txt LICENSE.txt README config/preamble.in; \
	  rm -rf src/alps/lattice* src/alps/model* src/alps/plot* src/alps/scheduler* lib/xml tool example test/lattice test/model test/plot doc/index.html doc/lattice doc/model doc/scheduler doc/tool; \
	  find . -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs ./config/update_preamble -l > /dev/null; \
	  find . -name '*\.txt' -or -name '*\.html' | xargs chmod 644
	if test "$@" = "snapshot"; then \
	  date=`date +%Y%m%d`; \
          mv -f $(distdir_light) alps-light-$$date; \
	  tar zcf alps-light-$$date.tar.gz alps-light-$$date; \
	  rm -rf alps-light-$$date; \
	else \
	  tar zcf $(distdir_light).tar.gz $(distdir_light); \
	  rm -rf $(distdir_light); \
	fi
	tar zxf $(distdir).tar.gz
	cd $(distdir); \
	  rm -f LICENSE-light.txt README-light config/preamble-light.in doc/index-light.html; \
	  find . -name '*\.txt' -or -name '*\.html' | xargs chmod 644
	rm -f $(distdir).tar.gz
	if test "$@" = "snapshot"; then \
	  date=`date +%Y%m%d`; \
          mv -f $(distdir) alps-$$date; \
	  tar zcf alps-$$date.tar.gz alps-$$date; \
	  rm -rf alps-$$date; \
	else \
	  tar zcf $(distdir).tar.gz $(distdir); \
	  rm -rf $(distdir); \
	fi

boost-snapshot boost-sandbox-snapshot : 
	if test -f $(top_srcdir)/config/make_snapshot; then \
	  $(SHELL) $(top_srcdir)/config/make_snapshot $@; \
	fi
