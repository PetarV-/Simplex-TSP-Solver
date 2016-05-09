CC = clang++ -Iinclude/
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic

tsp_solver : src/mst.cpp src/simplex.cpp src/tsp_solver.cpp
	$(CC) $(CFLAGS) src/mst.cpp src/simplex.cpp src/tsp_solver.cpp -o tsp_solver

.PHONY : clean
clean :
	-rm -f tsp_solver &> /dev/null
