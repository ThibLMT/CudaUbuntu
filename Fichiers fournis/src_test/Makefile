#############################################################################
# Makefile : demGCE code
# Discrete element modeling of Granular materials for Civil Engineering
# Developped by the Civil and Environmental Engineering Departement, MINES DOUAI
# Makefile developped by P. PIZETTE - last modif 26 june 2012 by
#############################################################################

# Compiler and options 
CC = gcc
FLAGS =  -O2 -W -Wall  -openmp 
FLAGSLINK =  -O2 -lm -openmp

# Names and paths of the executables and object and header directories
CIBLE = ./demGCE
CIBLE_BIN = ~/bin/demGCE
OBJ_DIR  = obj
H_DIR = .


# generate source and dependance files
SOURCES = $(shell ls *.c 2> /dev/null || true)
OBJETS = $(SOURCES:%.c=$(OBJ_DIR)/%.o)
DEPENDANCES = $(SOURCES:%.c=$(OBJ_DIR)/%.d)


# Staring compilation rules
all: $(OBJETS)
	$(CC) $^ -o $(CIBLE) $(FLAGSLINK)
	@echo "Generated executable :" $(CIBLE)

bin: $(OBJETS)
	$(CC) $^ -o $(CIBLE_BIN) $(FLAGSLINK)
	@echo "Generated executable :" $(CIBLE_BIN)

obj/%.o: %.c
	$(CC) -I$(H_DIR) $(FLAGS) -c -o $@ $<

%.d: %.c
	$(CC) -MM -I$(H_DIR) -MD -o $@ $<

-include $(DEPENDANCES)

# Clean options
clean:
	@echo "Remove object files"
	rm -f $(OBJ_DIR)/*.o
	@rm -f $(OBJ_DIR)/*.d 

cleanall: clean
	@echo "Remove generated excutables"
	rm -f $(CIBLE)
	rm -f $(CIBLE_TEST)

# Help rule to print help message
help:
	@echo "******************************************************"
	@echo "*      Available commands to compile demGCE code     *"
	@echo "******************************************************"
	@echo " make --> Compile files, then generate binary (demGCE)"
	@echo " make bin --> Create a binary in ~/bin" 
	@echo " make clean --> remove .o files"
	@echo " make cleanall --> remove .o files and binaries"
	@echo " make help --> Print help message"
	@echo " make doc --> Generate documentation files (not available yet)"
	@echo "******************************************************"


