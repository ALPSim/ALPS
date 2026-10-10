CXXFLAGS ?= -O2

.PHONY: clean distclean

hello: hello.C
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

clean distclean:
	$(RM) hello
