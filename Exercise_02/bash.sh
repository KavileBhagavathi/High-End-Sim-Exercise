rm -rf vtufiles/
cd src
rm *.vtu *.pvd | true
rm main
g++ --std=c++20 -o main main.cpp

./main

tar -cvf particle.tar.gz *.vtu *.pvd
mkdir ../vtufiles
mv *.vtu *.pvd ../vtufiles/
mv *.tar.gz ../