include $(ROOTSYS)/etc/Makefile.arch

CXXFLAGS += -I../.. -I.

ifeq ($(PLATFORM),macosx)
CXXFLAGS += -std=c++11
endif


TARGET=libTreeAnalysis.so

ASTREEINT   = AnalysisSupport/TreeInterface/src
ASUTILITIES = AnalysisSupport/Utilities/src
TREEREADING = TreeAnalyzer/src
#DATAFORMATS = DataFormats/src

SOURCE = $(wildcard $(ASTREEINT)/*.cc) $(wildcard $(ASUTILITIES)/*.cc) $(wildcard $(TREEREADING)/*.cc)
OBJ=$(join $(addsuffix ../obj/, $(dir $(SOURCE))), $(notdir $(SOURCE:.cc=.o))) 
DEPENDS=$(join $(addsuffix ../.dep/, $(dir $(SOURCE))), $(notdir $(SOURCE:.cc=.d)))


all:  $(TARGET)
	@true
clean:
	@-rm -f $(TARGET) $(OBJ) $(DEPENDS)
        
        
$(TARGET): $(OBJ)
	@echo "============="
	@echo "Linking the target $@"
	@echo "============="
ifeq ($(PLATFORM),macosx)
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
	@echo -- Link finished --
        
%.o : %.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@

$(ASTREEINT)/../obj/%.o : $(ASTREEINT)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(ASTREEINT)/../.dep/%.d: $(ASTREEINT)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(ASTREEINT)/../obj/$*.o^" > $@'
	
$(ASUTILITIES)/../obj/%.o : $(ASUTILITIES)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(ASUTILITIES)/../.dep/%.d: $(ASUTILITIES)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(ASUTILITIES)/../obj/$*.o^" > $@'

$(TREEREADING)/../obj/%.o : $(TREEREADING)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(TREEREADING)/../.dep/%.d: $(TREEREADING)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(TREEREADING)/../obj/$*.o^" > $@'
#	
#$(DATAFORMATS)/../obj/%.o : $(DATAFORMATS)/%.cc
#	@mkdir -p $(dir $@)
#	@echo "============="
#	@echo "Compiling $<"
#	$(CXX) $(CXXFLAGS) -c $<  -o $@
#$(DATAFORMATS)/../.dep/%.d: $(DATAFORMATS)/%.cc
#	@mkdir -p $(dir $@)
#	@echo "============="
#	@echo Building dependencies file for $*.o
#	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(DATAFORMATS)/../obj/$*.o^" > $@'
	
	

-include $(DEPENDS)