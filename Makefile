TARGET	= Alignment
CXX	= g++
$(TARGET):
	$(CXX) -x c $@.cc -O2 -lm -o $@

clean:
	rm -f Alignment