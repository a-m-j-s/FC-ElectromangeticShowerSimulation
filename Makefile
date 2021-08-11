ROOTINC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
cc = $(shell root-config --cxx)

%.o: %.C 
	@echo  -n "Compiling" $@ "\n"
	$(cc) -std=c++11 -c $(ROOTINC) -I. $^ -o $@

projeto: projeto.C Particle.o Propagator.o Funcoes.o Integrator.o Derivator.o Func1D.o FCmatrix.o FCmatrixFull.o Vec1.o EqSolver.o cFCgraphics.C
	$(cc) -std=c++11 -o projeto.exe projeto.C Particle.o Propagator.o Funcoes.o Integrator.o Derivator.o Func1D.o FCmatrix.o FCmatrixFull.o Vec1.o EqSolver.o cFCgraphics.C -I $(ROOTLIB) $(ROOTINC)
	./projeto.exe

clean:
	@echo cleaning...
	rm -f *.o
	rm -f *.exe

wall:
	$(eval cc += -Wall)
