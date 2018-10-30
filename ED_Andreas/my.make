INC=LAPSO_Req
CPPFLAGS=-ILAPSO_Req -I.
#CXXFLAGS= -O4 -Wall -fopenmp #-pg
CXXFLAGS= -g -Wall -fopenmp # debug flags
LDFLAGS=#-pg

OBJ= LaPSO.o anyoption.o CpuTimer.o VolVolume.o ED.o

build: lapso # $(TEST) 

lapso: main_test.o $(OBJ)
	$(CXX) $(CXXFLAGS) -o lapso main_test.o $(OBJ)


main.o: main.cpp $(INC)/LaPSO.hpp $(INC)/VolVolume.hpp

LaPSO.o: LaPSO.cpp $(INC)/LaPSO.hpp $(INC)/anyoption.h

clean:
	-rm -f *.o 


depend:
	makedepend -Y $(CPPFLAGS) -f my.make *.cpp >& /dev/null

# DO NOT DELETE

anyoption.o: LAPSO_Req/anyoption.h
CpuTimer.o: LAPSO_Req/CpuTimer.h
ED.o: LAPSO_Req/Random.h LAPSO_Req/VolVolume.hpp LAPSO_Req/LaPSO.hpp
ED.o: LAPSO_Req/CpuTimer.h ED.h
edge_disjoint.o: LAPSO_Req/LaPSO.hpp LAPSO_Req/CpuTimer.h
edge_disjoint.o: LAPSO_Req/VolVolume.hpp LAPSO_Req/Random.h
LaPSO.o: LAPSO_Req/LaPSO.hpp LAPSO_Req/CpuTimer.h LAPSO_Req/Random.h
LaPSO.o: LAPSO_Req/anyoption.h
main_test.o: ED.h LAPSO_Req/LaPSO.hpp LAPSO_Req/CpuTimer.h
