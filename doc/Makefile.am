pkgdocdir = @docdir@/alps

doc_ALPS = bitops.html config.html expression.html fixed_capacity.html functional.html index_general.html math.html parameterlist.html parameters.html random.html setup.html stringvalue.html typetraits.html valuetype.html vectormath.html vectortraits.html
doc_ALEA = alea/binning.html alea/index.html alea/observable.html alea/observableset.html alea/obsvalue.html alea/simpleobservable.html alea/simpleobseval.html
doc_LATTICE = lattice/coordinatelattice.html lattice/graph.html lattice/graph_traits.html lattice/graphconcepts.html lattice/graphproperties.html lattice/helperconcepts.html lattice/hypercubic.html lattice/index.html lattice/lattice.html lattice/latticeconcepts.html lattice/latticedescriptor.html lattice/latticegraph.html lattice/latticegraphdescriptor.html lattice/latticelibrary.html lattice/parity.html lattice/simplecell.html lattice/simplelattice.html lattice/unitcell.html
doc_MODEL = model/index.html
doc_OSIRIS = osiris/comm.html osiris/dump.html osiris/index.html osiris/mpdump.html osiris/os.html osiris/parameterlist.html osiris/parameters.html osiris/process.html osiris/stdmap.html osiris/stdpair.html osiris/stdstring.html osiris/stdvalarray.html osiris/stdvector.html osiris/valuetype.html osiris/xdrdump.html
doc_PARSER = parser/index.html parser/parser.html parser/xmlattributes.html parser/xmlstream.html
doc_SCHEDULER = scheduler/implementing.html scheduler/index.html scheduler/montecarlo.html scheduler/running.html scheduler/scheduler.html scheduler/signal.html
doc_TOOL = tool/index.html

doc_CONFIG =
dist_noinst_DATA = $(doc_ALPS) $(doc_OSIRIS) $(doc_PARSER)
if ALPS_FULL
dist_noinst_DATA += index.html $(doc_ALEA) $(doc_LATTICE) $(doc_MODEL) $(doc_SCHEDULER) $(doc_TOOL)
endif
if ALPS_LIGHT
dist_noinst_DATA += index-light.html
endif

install-data-local:
	@for f in $(doc_CONFIG) $(dist_noinst_DATA); do \
	  d=`echo $$f | sed -e s,/[^/]*$$,,`; \
	  if test "$$f" = "$$d"; then \
	    $(mkinstalldirs) $(DESTDIR)$(pkgdocdir); \
	  else \
	    $(mkinstalldirs) $(DESTDIR)$(pkgdocdir)/$$d; \
	  fi; \
	  if test -f $(srcdir)/$$f; then p=$(srcdir)/$$f; else p=$$f; fi; \
	  echo "$(INSTALL_DATA) $$p $(DESTDIR)$(pkgdocdir)/$$f"; \
	  $(INSTALL_DATA) $$p $(DESTDIR)$(pkgdocdir)/$$f; \
	done
