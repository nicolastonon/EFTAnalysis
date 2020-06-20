# Makefile for COMBINE codes

# $< = name of first dependancy
# $@ = name of target ("main")
# -c ask the compiler not to link
# -o specify target
#

#Define variables
CC=g++ #Compiler used
DEBUG = -g
CFLAGS = -c  $(DEBUG) #Compilation options
LFLAGS = $(DEBUG) -lTMVA
ROOTFLAGS = `root-config --glibs --cflags`

SCRIPT_PRODUCER = GenerateAllDatacards.exe
GENERATE_DATACARDS = GenerateTemplateDatacard.exe

#Instructions
all: $(GENERATE_DATACARDS) $(SCRIPT_PRODUCER)

$(GENERATE_DATACARDS):	Generate_Template_Datacard.cxx
	@g++ $< -o $@ $(ROOTFLAGS) $(LFLAGS)
	@echo "-- Created executable ./$(GENERATE_DATACARDS) --"

$(SCRIPT_PRODUCER):	ScriptProducer_GenerateAllDatacards.cxx
	@g++ $< -o $@ $(ROOTFLAGS) $(LFLAGS)
	@echo "-- Created executable ./$(SCRIPT_PRODUCER) --"

#Erase all objects and executable
clean:
	rm *.exe *.sh *~ Template_Datacard.txt
	rm -r datacards_TemplateFit
