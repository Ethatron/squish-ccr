
include config

SRC = colourclusterfit.cpp colourblock.cpp colourfit.cpp colourset.cpp colourrangefit.cpp singlecolourfit.cpp alpha.cpp maths.cpp squish.cpp	\
	paletteclusterfit.cpp paletteblock.cpp palettefit.cpp paletteset.cpp paletterangefit.cpp singlepalettefit.cpp

OBJ = $(SRC:%.cpp=%.o)

LIB = libsquish.a

all : $(LIB)

install : $(LIB)
	install squish.h $(INSTALL_DIR)/include
	install libsquish.a $(INSTALL_DIR)/lib

uninstall:
	$(RM) $(INSTALL_DIR)/include/squish.h
	$(RM) $(INSTALL_DIR)/lib/libsquish.a

$(LIB) : $(OBJ)
	$(AR) cr $@ $?
	ranlib $@

%.o : %.cpp
	$(CXX) $(CPPFLAGS) -I. $(CXXFLAGS) -o$@ -c $<

clean :
	$(RM) $(OBJ) $(LIB)



