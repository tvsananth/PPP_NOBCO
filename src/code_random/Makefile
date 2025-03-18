CCC = g++ -O3

all: main_FR_random_instances

main_FR_random_instances: func_cost.o func_SR.o selfish_routing_leblanc.o func_leblanc.o  main_FR_random_instances.o
	$(CCC) -o main_FR_random_instances func_cost.o func_SR.o selfish_routing_leblanc.o func_leblanc.o  main_FR_random_instances.o

main_FR_random_instances.cpp.o: ./main_FR_random_instances.cpp
	$(CCC) -c ./main_FR_random_instances.cpp

func_cost.o: ./func_cost.cpp ./func_cost.h
	$(CCC) -c ./func_cost.cpp

func_SR.o: ./func_SR.cpp ./func_SR.h
	$(CCC) -c ./func_SR.cpp

selfish_routing_leblanc.o: ./selfish_routing_leblanc.cpp ./selfish_routing_leblanc.h
	$(CCC) -c ./selfish_routing_leblanc.cpp

func_leblanc.o: ./func_leblanc.cpp ./func_leblanc.h
	$(CCC) -c ./func_leblanc.cpp

clean:
	rm -f *.o main_FR_random_instances
