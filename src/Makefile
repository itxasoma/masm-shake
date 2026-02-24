FC = gfortran
FFLAGS = -O2 -g -Wall -Wextra -fcheck=all -fbacktrace

TARGET = mdinamics

all: $(TARGET)

$(TARGET): parameters.o force.o motion.o io_teacher.o rdf_module.o main.o
	$(FC) $(FFLAGS) -o $@ $^

parameters.o: parameters.f90
	$(FC) $(FFLAGS) -c $<

force.o: force.f90 parameters.o
	$(FC) $(FFLAGS) -c $<

motion.o: motion.f90 parameters.o force.o
	$(FC) $(FFLAGS) -c $<

io_teacher.o: io_teacher.f90 parameters.o
	$(FC) $(FFLAGS) -c $<

rdf_module.o: rdf_module.f90 parameters.o force.o
	$(FC) $(FFLAGS) -c $<

main.o: main.f90 parameters.o force.o motion.o io_teacher.o rdf_module.o
	$(FC) $(FFLAGS) -c $<

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f *.o *.mod $(TARGET)

