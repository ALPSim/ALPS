SOURCES = \
	looper/alternating_tensor.h \
	looper/cluster.h \
	looper/correlation.h \
	looper/crop.h \
        looper/custom.h \
        looper/custom_impl.h \
	looper/divide_if_positive.h \
	looper/evaluator.h \
	looper/evaluator_impl.h \
	looper/graph.h \
	looper/graph_impl.h \
	looper/integer_range.h \
	looper/lapack.h \
	looper/lattice.h \
	looper/location.h \
	looper/location_impl.h \
	looper/matrix.h \
	looper/measurement.h \
	looper/model.h \
	looper/model_impl.h \
	looper/model_parameter.h \
	looper/montecarlo.h \
	looper/operator.h \
	looper/permutation.h \
	looper/power.h \
	looper/random_choice.h \
	looper/stiffness.h \
	looper/susceptibility.h \
	looper/temperature.h \
	looper/time.h \
	looper/type.h \
	looper/union_find.h \
	looper/version.h \
	looper/weight.h \
	looper/weight_impl.h

dist_noinst_DATA = doc/index.html

bin_PROGRAMS =
if BUILD_QMC_LOOPER
bin_PROGRAMS += loop
if HAVE_MPI
if HAVE_BOOST_MPI
bin_PROGRAMS += loop_mpi
endif
endif
endif

loop_SOURCES = loop.C loop_config.h loop_custom.C loop_model.C path_integral.C sse.C $(SOURCES)

loop_mpi_SOURCES = loop.C loop_config.h loop_custom.C loop_model.C path_integral.C sse.C $(SOURCES)
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@
