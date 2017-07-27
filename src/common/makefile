# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Modules
ModelPrecision.o : ModelPrecision.f90
	${FC} ${OPT} -c ModelPrecision.f90 -o $@

ConstantsDictionary.o : ModelPrecision.o ConstantsDictionary.f90
	${FC} ${OPT} -c ConstantsDictionary.f90 -o $@

LinkedList_Class.o : LinkedList_Class.f90
	${FC} ${OPT} -c LinkedList_Class.f90 -o $@

KeyRing_Class.o : CommonRoutines.o KeyRing_Class.f90
	${FC} ${OPT} -c KeyRing_Class.f90 -o $@

HashTable_Class.o : ConstantsDictionary.o LinkedList_Class.o HashTable_Class.f90
	${FC} ${OPT} -c HashTable_Class.f90 -o $@

CommonRoutines.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.f90
	${FC} ${OPT} -c CommonRoutines.f90 -o $@

Timing.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Timing.f90
	${FC} ${OPT} -c Timing.f90 -o $@
