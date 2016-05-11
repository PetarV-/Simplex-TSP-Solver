# Simplex-TSP-Solver
Iterative exact solver of the Travelling Salesman Problem, taking advantage of the Simplex Algorithm (by methods of Dantzig, Fulkerson and Johnson). The topmost layer of the solver also interfaces to a demo folder, which can be used for plotting the current solution as successive solutions are found.

The current form has been used for demonstrations at the [mgcsweek](http://www.csnedelja.mg.edu.rs) seminar at the [High School of Mathematics](http://www.mg.edu.rs) in Belgrade and the [Advanced Algorithms](http://www.cl.cam.ac.uk/teaching/1516/AdvAlgo/) lecture course for Part II of the Computer Science Tripos at the [University of Cambridge](http://www.cam.ac.uk).

## Building
Simply run

    $ make

while in the root directory of the repository, and the executable, `tsp_solver`, will be created therewith.

## Usage
In order to interface properly to a demo setting, four parameters need to be given right after executing

    $ ./tsp_solver

These are, in order:
- A path to the adjacency matrix specifying the input graph (`demo/adj_matrix` in this case);
- A path to the folder containing the demonstration .tex files (`demo/` in this case);
- File name of the .tex file containing the edge-drawing commands (`graph.tex` in this case);
- File name of the .tex file containing the map-drawing commands (`usa.tex` in this case);

Afterwards, an interactive window will open where one of the following commands may be performed repeatedly:
- `SOLVE`: launch the Simplex algorithm on the constraints obtained so far;
- `REM_LOOP`/`REM_LOOP_RNG`: add a cycle-eliminating constraint on a given set of nodes;
- `SET`: add a constraint that sets a particular edge's indicator variable to 0/1;
- `UNDO`: undo a number of previously added constraints (useful for branch&bound);
- `APPROX_MST`: perform a 2-approximation of the optimal tour using Prim's MST algorithm;
- `EXIT`: exit the program.

Upon each call to `SOLVE`, detailed information about all the useful iterations (that improve the objective value) will be printed on the terminal. Once a solution has been found, `pdflatex` will be called to regenerate the output PDF file (`demo/usa.pdf` in this case). Finally, if you are using Mac OS X, `qlmanage` will be called for automatically displaying the PDF.

## License
MIT
