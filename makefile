CC = clang++ -Iinclude/
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic $(DEBUG)

tsp_solver : src/simplex.cpp src/tsp_solver.cpp
	$(CC) $(CFLAGS) src/simplex.cpp src/tsp_solver.cpp -o tsp_solver

clean :
	-rm -f tsp_solver &> /dev/null
