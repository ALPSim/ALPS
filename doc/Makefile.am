dist_noinst_DATA =

#
#
#

BOOST_ROOT = @BOOST_INC_DIR@
ALPS_ROOT = $(top_srcdir)
ALPS_SRCDIR = $(ALPS_ROOT)/src/alps

DOCBOOK_DTD_DIR = @DOCBOOK_DTD_DIR@
DOCBOOK_XSL_DIR = @DOCBOOK_XSL_DIR@
FOP_DIR = @FOP_DIR@

QUICKBOOK = $(top_builddir)/tool/quickbook
XSLTPROC = @XSLTPROC@

DOC_ROOT = $(top_srcdir)/doc
DOC_SRCDIR = $(DOC_ROOT)/src
DOC_HTMLDIR = $(DOC_ROOT)/html
DOC_PDFDIR = $(DOC_ROOT)

DOCS =
if BUILD_DOC
DOCS += $(DOC_HTMLDIR)/index.html
if BUILD_PDF
DOCS += $(DOC_PDFDIR)/alps.pdf
endif
endif

docs : $(DOCS)

XMLDOCS = $(XMLDOCS_QBK:.qbk=.xml) $(XMLDOCS_GEN)

XMLDOCS_QBK = \
	$(DOC_SRCDIR)/alps.qbk \
	$(DOC_SRCDIR)/overview.qbk \
	$(DOC_SRCDIR)/getting_started.qbk \
	$(DOC_SRCDIR)/general.qbk \
	$(DOC_SRCDIR)/osiris.qbk \
	$(DOC_SRCDIR)/parser.qbk \
	$(DOC_SRCDIR)/alea.qbk \
	$(DOC_SRCDIR)/lattice.qbk \
	$(DOC_SRCDIR)/model.qbk \
	$(DOC_SRCDIR)/scheduler.qbk \
	$(DOC_SRCDIR)/tools.qbk

XMLDOCS_GEN =

AUTODOCS = $(DOC_SRCDIR)/general.autodoc
AUTODOC_GENERAL_SRC = \
	$(ALPS_SRCDIR)/bitops.h \
	$(ALPS_SRCDIR)/cctype.h \
	$(ALPS_SRCDIR)/config.h.in \
	$(ALPS_SRCDIR)/copyright.h \
	$(ALPS_SRCDIR)/expression.h \
	$(ALPS_SRCDIR)/factory.h \
	$(ALPS_SRCDIR)/fixed_capacity_deque.h \
	$(ALPS_SRCDIR)/fixed_capacity_fwd.h \
	$(ALPS_SRCDIR)/fixed_capacity_traits.h \
	$(ALPS_SRCDIR)/fixed_capacity_vector.h \
	$(ALPS_SRCDIR)/functional.h \
	$(ALPS_SRCDIR)/math.hpp \
	$(ALPS_SRCDIR)/matrix_as_vector.h \
	$(ALPS_SRCDIR)/multi_array.hpp \
	$(ALPS_SRCDIR)/parameterlist.h \
	$(ALPS_SRCDIR)/parameters.h \
	$(ALPS_SRCDIR)/random.h \
	$(ALPS_SRCDIR)/random/buffered_rng.h \
	$(ALPS_SRCDIR)/random/pseudo_des.h \
	$(ALPS_SRCDIR)/random/rngfactory.h \
	$(ALPS_SRCDIR)/random/seed.h \
	$(ALPS_SRCDIR)/random/uniform_on_sphere_n.h \
	$(ALPS_SRCDIR)/stringvalue.h \
	$(ALPS_SRCDIR)/typetraits.h \
	$(ALPS_SRCDIR)/vectorio.h \
	$(ALPS_SRCDIR)/vectormath.h \
	$(ALPS_SRCDIR)/vectortraits.h

AUTODOCS += $(DOC_SRCDIR)/osiris.autodoc
AUTODOC_OSIRIS_SRC = \
	$(ALPS_SRCDIR)/osiris.h \
	$(ALPS_SRCDIR)/osiris/buffer.h \
	$(ALPS_SRCDIR)/osiris/comm.h \
	$(ALPS_SRCDIR)/osiris/os.h \
	$(ALPS_SRCDIR)/osiris/process.h \
	$(ALPS_SRCDIR)/osiris/std/list.h \
	$(ALPS_SRCDIR)/osiris/std/map.h \
	$(ALPS_SRCDIR)/osiris/std/pair.h \
	$(ALPS_SRCDIR)/osiris/std/set.h \
	$(ALPS_SRCDIR)/osiris/std/string.h \
	$(ALPS_SRCDIR)/osiris/std/valarray.h \
	$(ALPS_SRCDIR)/osiris/std/vector.h

# 	$(ALPS_SRCDIR)/osiris/archivedump.h \
#	$(ALPS_SRCDIR)/osiris/dump.h \
#	$(ALPS_SRCDIR)/osiris/dumparchive.h \
#	$(ALPS_SRCDIR)/osiris/mpdump.h \
#	$(ALPS_SRCDIR)/osiris/xdrdump.h

AUTODOCS += $(DOC_SRCDIR)/parser.autodoc
AUTODOC_PARSER_SRC = \
	$(ALPS_SRCDIR)/xml.h \
	$(ALPS_SRCDIR)/parser/parser.h \
	$(ALPS_SRCDIR)/parser/xmlattributes.h \
	$(ALPS_SRCDIR)/parser/xmlhandler.h \
	$(ALPS_SRCDIR)/parser/xmlparser.h \
	$(ALPS_SRCDIR)/parser/xmlstream.h \
	$(ALPS_SRCDIR)/plot.h

AUTODOCS += $(DOC_SRCDIR)/alea.autodoc
AUTODOC_ALEA_SRC = \
	$(ALPS_SRCDIR)/alea.h \
	$(ALPS_SRCDIR)/alea/abstractbinning.h \
	$(ALPS_SRCDIR)/alea/abstractsimpleobservable.h \
	$(ALPS_SRCDIR)/alea/detailedbinning.h \
	$(ALPS_SRCDIR)/alea/hdf5.h \
	$(ALPS_SRCDIR)/alea/histogram.h \
	$(ALPS_SRCDIR)/alea/histogramdata.h \
	$(ALPS_SRCDIR)/alea/histogrameval.h \
	$(ALPS_SRCDIR)/alea/nan.h \
	$(ALPS_SRCDIR)/alea/nobinning.h \
	$(ALPS_SRCDIR)/alea/observable.h \
	$(ALPS_SRCDIR)/alea/observableset.h \
	$(ALPS_SRCDIR)/alea/obsvalue.h \
	$(ALPS_SRCDIR)/alea/output_helper.h \
	$(ALPS_SRCDIR)/alea/recordableobservable.h \
	$(ALPS_SRCDIR)/alea/signedobservable.h \
	$(ALPS_SRCDIR)/alea/simplebinning.h \
	$(ALPS_SRCDIR)/alea/simpleobsdata.h \
	$(ALPS_SRCDIR)/alea/simpleobservable.h \
	$(ALPS_SRCDIR)/alea/simpleobseval.h

AUTODOCS += $(DOC_SRCDIR)/lattice.autodoc
AUTODOC_LATTICE_SRC = \
	$(ALPS_SRCDIR)/lattice.h \
	$(ALPS_SRCDIR)/lattice/bond_compare.h \
	$(ALPS_SRCDIR)/lattice/boundary.h \
	$(ALPS_SRCDIR)/lattice/cell_traits.h \
	$(ALPS_SRCDIR)/lattice/coordinate_traits.h \
	$(ALPS_SRCDIR)/lattice/coordinatelattice.h \
	$(ALPS_SRCDIR)/lattice/dimensional_traits.h \
	$(ALPS_SRCDIR)/lattice/disorder.h \
	$(ALPS_SRCDIR)/lattice/graph.h \
	$(ALPS_SRCDIR)/lattice/graph_helper.h \
	$(ALPS_SRCDIR)/lattice/graph_traits.h \
	$(ALPS_SRCDIR)/lattice/graphproperties.h \
	$(ALPS_SRCDIR)/lattice/hypercubic.h \
	$(ALPS_SRCDIR)/lattice/lattice.h \
	$(ALPS_SRCDIR)/lattice/latticedescriptor.h \
	$(ALPS_SRCDIR)/lattice/latticegraph.h \
	$(ALPS_SRCDIR)/lattice/latticegraphdescriptor.h \
	$(ALPS_SRCDIR)/lattice/latticelibrary.h \
	$(ALPS_SRCDIR)/lattice/parity.h \
	$(ALPS_SRCDIR)/lattice/point_traits.h \
	$(ALPS_SRCDIR)/lattice/propertymap.h \
	$(ALPS_SRCDIR)/lattice/simplecell.h \
	$(ALPS_SRCDIR)/lattice/simplelattice.h \
	$(ALPS_SRCDIR)/lattice/unitcell.h

AUTODOCS += $(DOC_SRCDIR)/model.autodoc
AUTODOC_MODEL_SRC = \
	$(ALPS_SRCDIR)/model.h \
	$(ALPS_SRCDIR)/model/basisdescriptor.h \
	$(ALPS_SRCDIR)/model/basisstates.h \
	$(ALPS_SRCDIR)/model/blochbasisstates.h \
	$(ALPS_SRCDIR)/model/bondoperator.h \
	$(ALPS_SRCDIR)/model/bondterm.h \
	$(ALPS_SRCDIR)/model/globaloperator.h \
	$(ALPS_SRCDIR)/model/half_integer.h \
	$(ALPS_SRCDIR)/model/hamiltonian.h \
	$(ALPS_SRCDIR)/model/integer_state.h \
	$(ALPS_SRCDIR)/model/model_helper.h \
	$(ALPS_SRCDIR)/model/modellibrary.h \
	$(ALPS_SRCDIR)/model/operator.h \
	$(ALPS_SRCDIR)/model/operatordescriptor.h \
	$(ALPS_SRCDIR)/model/operatorsubstitution.h \
	$(ALPS_SRCDIR)/model/quantumnumber.h \
	$(ALPS_SRCDIR)/model/sign.h \
	$(ALPS_SRCDIR)/model/sitebasisdescriptor.h \
	$(ALPS_SRCDIR)/model/sitebasisstates.h \
	$(ALPS_SRCDIR)/model/siteoperator.h \
	$(ALPS_SRCDIR)/model/sitestate.h \
	$(ALPS_SRCDIR)/model/siteterm.h

AUTODOCS += $(DOC_SRCDIR)/scheduler.autodoc
AUTODOC_SCHEDULER_SRC = \
	$(ALPS_SRCDIR)/scheduler.h \
	$(ALPS_SRCDIR)/scheduler/factory.h \
	$(ALPS_SRCDIR)/scheduler/info.h \
	$(ALPS_SRCDIR)/scheduler/montecarlo.h \
	$(ALPS_SRCDIR)/scheduler/options.h \
	$(ALPS_SRCDIR)/scheduler/scheduler.h \
	$(ALPS_SRCDIR)/scheduler/task.h \
	$(ALPS_SRCDIR)/scheduler/types.h \
	$(ALPS_SRCDIR)/scheduler/worker.h

# quickbook -> boostbook

SUFFIXES = .xml .qbk
.qbk.xml :
	$(QUICKBOOK) $<

# *.h -> *.autodoc

$(DOC_SRCDIR)/general.autodoc : autodoc $(AUTODOC_GENERAL_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_GENERAL_SRC)
$(DOC_SRCDIR)/osiris.autodoc : autodoc $(AUTODOC_OSIRIS_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_OSIRIS_SRC)
$(DOC_SRCDIR)/parser.autodoc : autodoc $(AUTODOC_PARSER_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_PARSER_SRC)
$(DOC_SRCDIR)/alea.autodoc : autodoc $(AUTODOC_ALEAS_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_ALEA_SRC)
$(DOC_SRCDIR)/lattice.autodoc : autodoc $(AUTODOC_LATTICE_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_LATTICE_SRC)
$(DOC_SRCDIR)/model.autodoc : autodoc $(AUTODOC_MODEL_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_MODEL_SRC)
$(DOC_SRCDIR)/scheduler.autodoc : autodoc $(AUTODOC_SCHEDULER_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_SCHEDULER_SRC)

# boostbook -> docbook

$(DOC_SRCDIR)/alps.docbook : bb2db $(XMLDOCS) $(AUTODOCS)
	$(SHELL) bb2db $(DOC_SRCDIR)/alps.xml $@

# docbook -> html

$(DOC_HTMLDIR)/index.html : db2html $(DOC_SRCDIR)/alps.docbook html-setup
	$(SHELL) db2html $(DOC_SRCDIR)/alps.docbook $(DOC_HTMLDIR)/

html-setup :
	$(mkdir_p) $(DOC_HTMLDIR)
	cp -fp $(DOC_SRCDIR)/alps_logo.gif $(DOC_HTMLDIR)
	cp -fp $(BOOST_ROOT)/doc/html/boostbook.css $(DOC_HTMLDIR)
	$(mkdir_p) $(DOC_HTMLDIR)/images
	cp -fp $(BOOST_ROOT)/doc/html/images/*.png $(DOC_HTMLDIR)/images

# docbook -> pdf

$(DOC_PDFDIR)/alps.pdf : db2pdf $(DOC_SRCDIR)/alps.docbook
	$(SHELL) db2pdf $(DOC_SRCDIR)/alps.docbook $@

docclean :
	rm -f $(XMLDOCS_QBK:.qbk=.xml) $(AUTODOCS)
	rm -f $(DOC_SRCDIR)/*.docbook $(DOC_SRCDIR)/*.fo $(DOC_PDFDIR)/alps.pdf
	rm -rf $(DOC_HTMLDIR)

#
#
#

pkgdocdir = @docdir@/alps

doc_ALPS = bitops.html config.html expression.html fixed_capacity.html functional.html index_general.html math.html parameterlist.html parameters.html random.html setup.html stringvalue.html typetraits.html valuetype.html vectormath.html vectortraits.html
doc_ALEA = \
	alea/abstractsimpleobservable.html \
	alea/binning.html \
	alea/index.html \
	alea/observable.html \
	alea/observableset.html \
	alea/obsvalue.html \
	alea/recordableobservable.html \
	alea/signedobservable.html \
	alea/simpleobservable.html \
	alea/simpleobseval.html
doc_LATTICE = \
	lattice/boundary.html \
	lattice/coordinatelattice.html \
	lattice/disorder.html \
	lattice/graph.html \
	lattice/graph_helper.html \
	lattice/graph_traits.html \
	lattice/graphconcepts.html \
	lattice/graphproperties.html \
	lattice/helperconcepts.html \
	lattice/hypercubic.html \
	lattice/index.html \
	lattice/lattice.html \
	lattice/latticeconcepts.html \
	lattice/latticedescriptor.html \
	lattice/latticegraph.html \
	lattice/latticegraphdescriptor.html \
	lattice/latticelibrary.html \
	lattice/parity.html \
	lattice/propertymap.html \
	lattice/simplecell.html \
	lattice/simplelattice.html \
	lattice/unitcell.html
doc_MODEL = model/index.html
doc_OSIRIS = osiris/comm.html osiris/dump.html osiris/index.html osiris/mpdump.html osiris/os.html osiris/parameterlist.html osiris/parameters.html osiris/process.html osiris/stdmap.html osiris/stdpair.html osiris/stdstring.html osiris/stdvalarray.html osiris/stdvector.html osiris/valuetype.html osiris/xdrdump.html
doc_PARSER = parser/index.html parser/parser.html parser/xmlattributes.html parser/xmlstream.html
doc_SCHEDULER = scheduler/implementing.html scheduler/index.html scheduler/montecarlo.html scheduler/running.html scheduler/scheduler.html scheduler/signal.html
doc_TOOL = tool/index.html

dist_noinst_DATA += $(doc_ALPS) $(doc_OSIRIS) $(doc_PARSER)
if ALPS_FULL
dist_noinst_DATA += index.html $(doc_ALEA) $(doc_LATTICE) $(doc_MODEL) $(doc_SCHEDULER) $(doc_TOOL)
endif
if ALPS_LIGHT
dist_noinst_DATA += index-light.html
endif

install-data-local:
	@for f in $(dist_noinst_DATA); do \
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
