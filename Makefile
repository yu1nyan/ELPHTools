OBJDIR ?= objs
DEPDIR ?= $(OBJDIR)/dep
LIBS = `root-config --libs`
OPT ?= 2
DEBUG ?= -g
STD ?= gnu++0x

GENDEPFLAGS = -MMD -MP -MF $(DEPDIR)/$(@F).d
CXXFLAGS += $(CXXDEFS)
CXXFLAGS += -O$(OPT)
CXXFLAGS += -march=native
CXXFLAGS += $(DEBUG)
CXXFLAGS += -Wall -Wextra
CXXFLAGS += $(patsubst %,-I%,$(EXTRAINCDIRS))
CXXFLAGS += -std=$(STD)
CXXFLAGS += `root-config --cflags`
MATHLIB = -lm
LIBS += $(MATHLIB)
ALL_CXXFLAGS = $(GENDEPFLAGS) $(CXXFLAGS)
LDFLAGS = $(LIBS)

CXX := g++
LD := g++
REMOVE := rm -f
REMOVEDIR := rm -rf


TARGET = gap_check_3E run_calib run_proto run_1cube run_2cubes run_9cubes run_3Dprint run_proto_gapcheck run_proto_my hodo_cellHitRate
SRC = gap_check_3E.cc run_calib.cc run_proto.cc run_1cube.cc run_2cubes.cc run_9cubes.cc run_3Dprint.cc ToScintiCh.cc run_proto_gapcheck.cc run_proto_my.cc hodo_cellHitRate.cc
OBJ = $(OBJDIR)/gap_check_3E.o $(OBJDIR)/run_calib.o $(OBJDIR)/run_proto.o $(OBJDIR)/run_1cube.o $(OBJDIR)/run_2cubes.o $(OBJDIR)/run_9cubes.o $(OBJDIR)/run_3Dprint.o $(OBJDIR)/ToScintiCh.o $(OBJDIR)/run_gapcheck.o $(OBJDIR)/run_proto_my.o $(OBJDIR)/hodo_cellHitRate.o

all: $(TARGET)

hodo_cellHitRate: $(OBJDIR)/hodo_cellHitRate.o
	$(LD) -o $@ $^ $(LDFLAGS)

run_proto_my: $(OBJDIR)/run_proto_my.o $(OBJDIR)/ToScintiCh.o
	$(LD) -o $@ $^ $(LDFLAGS)

run_proto_gapcheck: $(OBJDIR)/run_proto_gapcheck.o
	$(LD) -o $@ $^ $(LDFLAGS)

gap_check_3E: $(OBJDIR)/gap_check_3E.o
	$(LD) -o $@ $^ $(LDFLAGS)

run_calib: $(OBJDIR)/run_calib.o $(OBJDIR)/ToScintiCh.o
	$(LD) -o $@ $^ $(LDFLAGS)

run_proto: $(OBJDIR)/run_proto.o $(OBJDIR)/ToScintiCh.o
	$(LD) -o $@ $^ $(LDFLAGS)

run_1cube: $(OBJDIR)/run_1cube.o
	$(LD) -o $@ $^ $(LDFLAGS)

run_2cubes: $(OBJDIR)/run_2cubes.o $(OBJDIR)/ToScintiCh.o
	$(LD) -o $@ $^ $(LDFLAGS)

run_9cubes: $(OBJDIR)/run_9cubes.o $(OBJDIR)/ToScintiCh.o
	$(LD) -o $@ $^ $(LDFLAGS)

run_3Dprint: $(OBJDIR)/run_3Dprint.o $(OBJDIR)/ToScintiCh.o
	$(LD) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o : %.cc
	$(CXX) -c $(ALL_CXXFLAGS) $<  -o  $@

clean:
	-$(REMOVE) $(TARGET)
	-$(REMOVE) $(OBJ)
	-$(REMOVEDIR) $(DEPDIR)
	-$(REMOVEDIR) $(OBJDIR)

$(shell mkdir $(OBJDIR) 2>/dev/null)
-include $(shell mkdir $(DEPDIR) 2>/dev/null) $(wildcard $(DEPDIR)/*.d)

.PRECIOUS: $(OBJ)
