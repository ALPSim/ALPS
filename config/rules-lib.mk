# library

$(LIB) : $(OBJ)
	$(ARXX) $(cxxflags) $(ldflags) -o $@ $(OBJ) $(libs)
