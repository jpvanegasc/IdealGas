DEBUG = -g -Wall -fsanitize=address -fsanitize=leak -fsanitize=undefined
VAL = --track-origins=yes --leak-check=full
CACH = --tool=cachegrind
C_STD = -std=c++11
all : main.o  *.h

main.o : main.cpp *.h
	g++ $(C_STD) -O2 $< -o $@
	time ./$@

.PHONY : animation
animation : main.cpp *.h animation.hpp
	g++ $(C_STD) -O2 $<
	time ./a.out animate | gnuplot

.PHONY : save_gif
save_gif : main.cpp *.h animation.hpp
	g++ $(C_STD) -O2 $<
	time ./a.out animate save| gnuplot

.PHONY : debug
debug : main.cpp *.h
	g++ $(C_STD) $(DEBUG) $<
	./a.out

.PHONY : valgrind
valgrind : main.cpp *.h
	g++ $(C_STD) -g $<
	valgrind $(VAL) ./a.out

.PHONY : cachegrind
cachegrind : main.cpp *.h
	g++ $(C_STD) -g $<
	valgrind $(CACH) ./a.out

.PHONY : gprof
gprof : main.cpp *.h
	g++ $(C_STD) -Wall -pg $<
	./a.out
	gprof a.out gmon.out > analysis.txt
	most analysis.txt

.PHONY : perf
perf : main.cpp *.h
	g++ $(C_STD) -Wall -pg $<
	perf record ./a.out ; perf report

.PHONY : clean
clean :
	rm -f *.o *.txt *.out* debug main *.data

.PHONY : plot
plot : data.txt
	gnuplot plot.gp
