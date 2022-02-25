all: search cluster unit_test

CC = g++ -O2
FLAGS = -c #-g -std=c++14 #-Wall


OBJS1 =  main_search.o Grid.o Utilities.o Dataset.o HashTable.o LSH.o Hypercube.o
OBJS2 =  main_cluster.o Utilities.o Dataset.o HashTable.o LSH.o Grid.o Cluster.o BinaryTree.o
OBJS3 = Unit_test.o Utilities.o

Headers = Grid.h Utilities.h Dataset.h HashTable.h LSH.h Hypercube.h Cluster.h BinaryTree.o


EXEC1 = search
EXEC2 = cluster
EXEC3 = unit_test

# nearest neighbor
$(EXEC1): $(OBJS1) $(Headers) 
	$(CC)  $(OBJS1) Continuous_Frechet_github/*.cpp -o $(EXEC1)

# cluster
$(EXEC2): $(OBJS2) $(Headers) 
	$(CC)  $(OBJS2) Continuous_Frechet_github/*.cpp -o $(EXEC2)

# Cunit test
$(EXEC3): $(OBJS3) $(Headers) 
	$(CC)  $(OBJS3) Continuous_Frechet_github/*.cpp -o $(EXEC3) -lcunit

Utilities.o: Utilities.cpp
	$(CC) $(FLAGS) Utilities.cpp

Dataset.o: Dataset.cpp
	$(CC) $(FLAGS) Dataset.cpp

Grid.o: Grid.cpp
	$(CC) $(FLAGS) Grid.cpp

HashTable.o: HashTable.cpp
	$(CC) $(FLAGS) HashTable.cpp

LSH.o: LSH.cpp
	$(CC) $(FLAGS) LSH.cpp

Hypercube.o: Hypercube.cpp
	$(CC) $(FLAGS) Hypercube.cpp

Cluster.o: Cluster.cpp
	$(CC) $(FLAGS) Cluster.cpp

BinaryTree.o: BinaryTree.cpp
	$(CC) $(FLAGS) BinaryTree.cpp


main_search.o: main_search.cpp
	$(CC) $(FLAGS) main_search.cpp

main_cluster.o: main_cluster.cpp
	$(CC) $(FLAGS) main_cluster.cpp

Unit_test.o: Unit_test.cpp
	$(CC) $(FLAGS) Unit_test.cpp


clean:
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) $(OBJS1) $(OBJS2) $(OBJS3)



# execute nearest neighbor with input files
#./search –i <input file> –q <query file> –k <int> -L <int> -M <int> -probes <int> -ο <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete or continuous | only for –algorithm Frechet> -delta <double> 
run_NN_LSH:
	./search -i Datasets/Datasets/nasd_input.csv -q Datasets/Datasets/nasd_query.csv -o output_NN_LSH.txt -algorithm LSH 

run_NN_Hypercube:
	./search -i Datasets/Datasets/nasd_input.csv -q Datasets/Datasets/nasd_query.csv -o output_NN_Hypercube.txt -algorithm Hypercube 	

run_NN_Frechet_Discrete:
	./search -i Datasets/Datasets/nasd_input.csv -q Datasets/Datasets/nasd_query.csv -o output_NN_DFD.txt -algorithm Frechet -metric discrete -delta 2

run_NN_Frechet_Continuous:
	./search -i Datasets/Datasets/nasd_input.csv -q Datasets/Datasets/nasd_query.csv -o output_CFD.txt -algorithm Frechet -metric continuous -delta 1



# execute clustering with input files
run_cluster_LSH_frechet:
	./cluster -i Datasets/Datasets/nasd_input.csv -c cluster.conf -o output_cluster_LSH_Frechet.txt -update Mean Frechet -assignment LSH -complete -silhouette

run_cluster_Classic_frechet:
	./cluster -i Datasets/Datasets/nasd_input.csv -c cluster.conf -o output_cluster_Classic_Frechet.txt -update Mean Frechet -assignment Classic -complete -silhouette

run_cluster_Classic_Vector:
	./cluster -i Datasets/Datasets/nasd_input.csv -c cluster.conf -o output_cluster_Classic_Vector.txt -update Mean Vector -assignment Classic -complete -silhouette

run_cluster_LSH_Vector:
	./cluster -i Datasets/Datasets/nasd_input.csv -c cluster.conf -o output_cluster_LSH_Vector.txt -update Mean Vector -assignment LSH -complete -silhouette

run_cluster_Hypercube_Vector:
	./cluster -i Datasets/Datasets/nasd_input.csv -c cluster.conf -o output_cluster_Hypercube_Vector.txt -update Mean Vector -assignment Hypercube -complete -silhouette



run_unit_test:
	./unit_test
