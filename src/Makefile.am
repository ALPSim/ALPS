headers_CONFIG = alps/config.h @BOOST_USER_CONFIG_H@
headers_BOOST = boost/function_objects.hpp
headers_ALPS = alps/bitops.h alps/cctype.h alps/copyright.h alps/expression.h alps/expression_impl.h alps/expression2.h alps/expression2_impl.h alps/evaluator.h alps/factory.h alps/functional.h alps/math.hpp alps/multi_array.hpp alps/parameterlist.h alps/parameters.h alps/random.h alps/stringvalue.h alps/typetraits.h alps/vectorio.h alps/vectormath.h alps/vectortraits.h
headers_ALEA = alps/alea.h alps/alea/detailedbinning.h alps/alea/histogram.h alps/alea/nan.h alps/alea/nobinning.h alps/alea/observable.h alps/alea/observableset.h alps/alea/obsvalue.h alps/alea/simplebinning.h alps/alea/simpleobsdata.h alps/alea/simpleobservable.h alps/alea/simpleobseval.h alps/alea/signedobservable.h alps/alea/hdf5.h alps/alea/recordableobservable.h alps/alea/abstractsimpleobservable.h alps/alea/abstractbinning.h alps/alea/output_helper.h
headers_FIXEDCAP = alps/fixed_capacity_deque.h alps/fixed_capacity_traits.h alps/fixed_capacity_fwd.h alps/fixed_capacity_vector.h alps/fixed_capacity/checking.h alps/fixed_capacity/deque_detail.h alps/fixed_capacity/uninitialized_array.h
headers_OSIRIS = alps/osiris.h alps/osiris/buffer.h alps/osiris/comm.h alps/osiris/dump.h alps/osiris/mpdump.h alps/osiris/os.h alps/osiris/process.h alps/osiris/xdrdump.h alps/osiris/std/impl.h alps/osiris/std/list.h alps/osiris/std/map.h alps/osiris/std/pair.h alps/osiris/std/set.h alps/osiris/std/string.h alps/osiris/std/valarray.h alps/osiris/std/vector.h
headers_LATTICE = alps/lattice.h alps/lattice/boundary.h alps/lattice/cell_traits.h alps/lattice/coordinate_traits.h alps/lattice/coordinatelattice.h alps/lattice/dimensional_traits.h alps/lattice/disorder.h alps/lattice/graph.h alps/lattice/graph_helper.h alps/lattice/graph_traits.h alps/lattice/graphproperties.h alps/lattice/hypercubic.h alps/lattice/lattice.h alps/lattice/latticedescriptor.h alps/lattice/latticegraph.h alps/lattice/latticegraphdescriptor.h alps/lattice/latticelibrary.h alps/lattice/parity.h alps/lattice/point_traits.h alps/lattice/propertymap.h alps/lattice/simplecell.h alps/lattice/simplelattice.h alps/lattice/unitcell.h
headers_MODEL = alps/model.h alps/model/basisdescriptor.h alps/model/basisstates.h alps/model/bondterm.h alps/model/half_integer.h alps/model/hamiltonian.h alps/model/integer_state.h alps/model/model_helper.h alps/model/modellibrary.h alps/model/operator.h alps/model/quantumnumber.h alps/model/sign.h alps/model/sitebasisdescriptor.h alps/model/sitebasisstates.h alps/model/sitestate.h alps/model/siteterm.h
headers_PARSER = alps/xml.h alps/parser/parser.h alps/parser/xmlattributes.h alps/parser/xmlhandler.h alps/parser/xmlparser.h alps/parser/xmlstream.h alps/parser/xslt.h
headers_RANDOM = alps/random.h alps/random/pseudo_des.h alps/random/random.h alps/random/rngfactory.h
headers_SCHEDULER = alps/scheduler.h alps/scheduler/factory.h alps/scheduler/info.h alps/scheduler/montecarlo.h alps/scheduler/options.h alps/scheduler/scheduler.h alps/scheduler/signal.hpp alps/scheduler/task.h alps/scheduler/types.h alps/scheduler/worker.h
headers_BINDINGS = alps/bindings/gels.hpp alps/bindings/heev.hpp alps/bindings/lapack.h alps/bindings/lapack_names.h alps/bindings/syev.hpp

headers = $(headers_BOOST) $(headers_ALPS) $(headers_FIXEDCAP) $(headers_OSIRIS) $(headers_PARSER) $(headers_RANDOM) $(headers_BINDINGS)
if ALPS_FULL
headers += $(headers_ALEA) $(headers_LATTICE) $(headers_MODEL) $(headers_SCHEDULER)
endif

noinst_HEADERS = $(headers)

install-data-local:
	@for f in $(headers_CONFIG) $(headers); do \
	  d=`echo $$f | sed -e s,/[^/]*$$,,`; \
	  if test "$$f" = "$$d"; then \
	    $(mkinstalldirs) $(DESTDIR)$(includedir); \
	  else \
	    $(mkinstalldirs) $(DESTDIR)$(includedir)/$$d; \
	  fi; \
	  if test -f $(srcdir)/$$f; then p=$(srcdir)/$$f; else p=$$f; fi; \
	  echo "$(INSTALL_DATA) $$p $(DESTDIR)$(includedir)/$$f"; \
	  $(INSTALL_DATA) $$p $(DESTDIR)$(includedir)/$$f; \
	done
uninstall-local:
	@for f in $(headers_CONFIG) $(headers); do \
	  echo "rm -f $(DESTDIR)$(includedir)/$$f"; \
	  rm -f $(DESTDIR)$(includedir)/$$f; \
	done
