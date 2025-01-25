# which FFT?
FFTOPTION= -D FFTW3
# link editor otions to use FFT library:
FFTLIBOPTIONS= -lfftw3 -lm	
# directory containing FFT header file :
FFTINCDIR= -I /usr/include
# directory containing FFT library or librariies :
FFTLIBDIR= -L /usr/lib

INC_GL = -I /usr/include/i386_linux_gnu
LIB_GL = -lGL -lGLU -lglut

INC_X = -I /usr/include/X11
LIB_X = -L /usr/lib/X11

CXXFLAGS   = -Wall -O3 -g #-fast -Wefc++ 
CXX        = g++ -w $(FFTINCDIR) $(INC_GL) $(INC_X)

SRCS = $(wildcard SOURCE/*/*.cpp)
BIN  = $(wildcard SOURCE/*/*.o)
OBJS= $(SRCS:.cpp=.o)
EXEC= run_FFT

all: $(EXEC)

$(EXEC) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FFTLIBDIR) $(FFTLIBOPTIONS) $(LIB_GL) $(LIB_X)  -lrt	

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -o $@ -c  $< $(LIB_GL) $(LIB_X) -lrt
	  
SOURCE/MAIN/main.o: $(SRCS)

clean:
	-rm $(BIN)		
	-rm $(EXEC)


