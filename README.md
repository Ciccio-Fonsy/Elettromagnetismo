g++ generator.cpp particle.cpp $(root-config --cflags --libs) -o generator

g++ analysis.cpp $(root-config --cflags --libs) -o analysis
