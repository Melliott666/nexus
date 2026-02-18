rm -rf build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../ ../ -DCMAKE_PREFIX_PATH="/opt/homebrew/opt/qt@5";
make -j8
make install
cd ..
