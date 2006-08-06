dist_noinst_DATA =

BOOST_DIR = @BOOST_DIR@
ALPS_SRCDIR = $(top_srcdir)/src/alps
DOC_SRCDIR = $(top_srcdir)/doc/src
DOC_XMLDIR = src
DOC_HTMLDIR = html
DOC_PDFDIR = .

DOCBOOK_DTD_DIR = @DOCBOOK_DTD_DIR@
DOCBOOK_XSL_DIR = @DOCBOOK_XSL_DIR@
FOP_DIR = @FOP_DIR@

QUICKBOOK = $(top_builddir)/tool/quickbook
XSLTPROC = @XSLTPROC@

DOCS =
if BUILD_DOC
DOCS += $(DOC_HTMLDIR)/index.html
if BUILD_PDF
DOCS += $(DOC_PDFDIR)/alps.pdf
endif
endif

.PHONY : doc docs
doc docs : $(DOCS)

dist_noinst_DATA += $(XMLDOCS_QBK) \
	src/alps_logo.gif \
	src/images/home.png src/images/next.png src/images/note.png src/images/prev.png src/images/up.png \
	xsl/html.xsl xsl/navbar.xsl

XMLDOCS_QBK = \
	src/alps.qbk \
	src/overview.qbk \
	src/getting_started.qbk \
	src/general.qbk \
	src/osiris.qbk \
	src/parser.qbk \
	src/alea.qbk \
	src/lattice.qbk \
	src/model.qbk \
	src/scheduler.qbk \
	src/tools.qbk

XMLDOCS = $(XMLDOCS_QBK:.qbk=.xml)

AUTODOCS = $(DOC_XMLDIR)/general.autodoc
AUTODOC_GENERAL_SRC = \
	$(ALPS_SRCDIR)/bitops.h \
	$(ALPS_SRCDIR)/cctype.h \
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
# $(ALPS_SRCDIR)/config.h.in

AUTODOCS += $(DOC_XMLDIR)/osiris.autodoc
AUTODOC_OSIRIS_SRC = \
	$(ALPS_SRCDIR)/osiris.h \
	$(ALPS_SRCDIR)/osiris/archivedump.h \
	$(ALPS_SRCDIR)/osiris/buffer.h \
	$(ALPS_SRCDIR)/osiris/comm.h \
	$(ALPS_SRCDIR)/osiris/dump.h \
	$(ALPS_SRCDIR)/osiris/dumparchive.h \
	$(ALPS_SRCDIR)/osiris/mpdump.h \
	$(ALPS_SRCDIR)/osiris/os.h \
	$(ALPS_SRCDIR)/osiris/process.h \
	$(ALPS_SRCDIR)/osiris/std/list.h \
	$(ALPS_SRCDIR)/osiris/std/map.h \
	$(ALPS_SRCDIR)/osiris/std/pair.h \
	$(ALPS_SRCDIR)/osiris/std/set.h \
	$(ALPS_SRCDIR)/osiris/std/string.h \
	$(ALPS_SRCDIR)/osiris/std/valarray.h \
	$(ALPS_SRCDIR)/osiris/std/vector.h \
	$(ALPS_SRCDIR)/osiris/xdrdump.h

AUTODOCS += $(DOC_XMLDIR)/parser.autodoc
AUTODOC_PARSER_SRC = \
	$(ALPS_SRCDIR)/xml.h \
	$(ALPS_SRCDIR)/parser/parser.h \
	$(ALPS_SRCDIR)/parser/xmlattributes.h \
	$(ALPS_SRCDIR)/parser/xmlhandler.h \
	$(ALPS_SRCDIR)/parser/xmlparser.h \
	$(ALPS_SRCDIR)/parser/xmlstream.h \
	$(ALPS_SRCDIR)/plot.h

AUTODOCS += $(DOC_XMLDIR)/alea.autodoc
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

AUTODOCS += $(DOC_XMLDIR)/lattice.autodoc
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

AUTODOCS += $(DOC_XMLDIR)/model.autodoc
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

AUTODOCS += $(DOC_XMLDIR)/scheduler.autodoc
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

src/alps.xml : $(DOC_SRCDIR)/alps.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/overview.xml : $(DOC_SRCDIR)/overview.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/getting_started.xml : $(DOC_SRCDIR)/getting_started.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/general.xml : $(DOC_SRCDIR)/general.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/osiris.xml : $(DOC_SRCDIR)/osiris.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/parser.xml : $(DOC_SRCDIR)/parser.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/alea.xml : $(DOC_SRCDIR)/alea.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/lattice.xml : $(DOC_SRCDIR)/lattice.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/model.xml : $(DOC_SRCDIR)/model.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/scheduler.xml : $(DOC_SRCDIR)/scheduler.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)
src/tools.xml : $(DOC_SRCDIR)/tools.qbk
	test -f $(@:.xml=.qbk) || $(LN_S) $< src
	$(QUICKBOOK) --input-file $(@:.xml=.qbk)

# *.h -> *.autodoc

$(DOC_XMLDIR)/general.autodoc : $(AUTODOC_GENERAL_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_GENERAL_SRC)
$(DOC_XMLDIR)/osiris.autodoc : $(AUTODOC_OSIRIS_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_OSIRIS_SRC)
$(DOC_XMLDIR)/parser.autodoc : $(AUTODOC_PARSER_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_PARSER_SRC)
$(DOC_XMLDIR)/alea.autodoc : $(AUTODOC_ALEAS_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_ALEA_SRC)
$(DOC_XMLDIR)/lattice.autodoc : $(AUTODOC_LATTICE_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_LATTICE_SRC)
$(DOC_XMLDIR)/model.autodoc : $(AUTODOC_MODEL_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_MODEL_SRC)
$(DOC_XMLDIR)/scheduler.autodoc : $(AUTODOC_SCHEDULER_SRC)
	$(SHELL) autodoc $@ $(AUTODOC_SCHEDULER_SRC)

# boostbook -> docbook

$(DOC_XMLDIR)/alps.docbook : $(AUTODOCS) $(XMLDOCS)
	$(SHELL) bb2db $(DOC_XMLDIR)/alps.xml $@

# docbook -> html

CSS = $(DOC_HTMLDIR)/boostbook.css
HTML_IMAGES = $(DOC_HTMLDIR)/alps_logo.gif \
	$(DOC_HTMLDIR)/images/home.png \
	$(DOC_HTMLDIR)/images/next.png \
	$(DOC_HTMLDIR)/images/note.png \
	$(DOC_HTMLDIR)/images/prev.png \
	$(DOC_HTMLDIR)/images/up.png

$(DOC_HTMLDIR)/index.html : $(DOC_XMLDIR)/alps.docbook $(CSS) $(HTML_IMAGES)
	$(SHELL) db2html $(DOC_XMLDIR)/alps.docbook $(DOC_HTMLDIR)/

$(DOC_HTMLDIR)/boostbook.css : $(BOOST_DIR)/doc/html/boostbook.css
	$(mkdir_p) $(DOC_HTMLDIR)
	cp -fp $< $@

$(DOC_HTMLDIR)/alps_logo.gif : $(DOC_SRCDIR)/alps_logo.gif
	$(mkdir_p) $(DOC_HTMLDIR)
	cp -fp $< $@

$(DOC_HTMLDIR)/images/%.png : $(DOC_SRCDIR)/images/%.png
	$(mkdir_p) $(DOC_HTMLDIR)/images
	cp -fp $< $@

# docbook -> pdf

$(DOC_PDFDIR)/alps.pdf : $(DOC_XMLDIR)/alps.docbook
	$(SHELL) db2pdf $(DOC_XMLDIR)/alps.docbook $@

# clean

doc-distclean :
	rm -f $(XMLDOCS_GEN) $(AUTODOCS)
	rm -f $(DOC_PDFDIR)/alps.pdf
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
