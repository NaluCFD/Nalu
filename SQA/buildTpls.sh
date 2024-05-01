BASE_DIR=/home/whoAmI/gitHubWork/scratch_build/packages

echo ............OPERATTING ON: $nalu_install_dir (part of your module paths)
echo ............SCRIPT EXPECTS /tarFiles with zlib, hdf5, pnetcdf, and netcdf tars

cd $BASE_DIR
echo ............create a previous directory for old TPLs
mkdir previous

echo YAML ............CLEAN
cd yaml-cpp/
rm -rf build
mkdir build
cd build
echo YAML ............CONFIG
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=$nalu_install_dir/yaml/0.6.2 ..
echo YAML ............BUILD
make -j 32
echo YAML ............INSTALL
make install

cd $BASE_DIR
echo Zlib ............CLEAN
mv zlib-1.2.11 previous
cp tarFiles/zlib-1.2.11.tar.gz ./
tar -zxvf zlib-1.2.11.tar.gz
rm -rf zlib-1.2.11.tar.gz
cd zlib-1.2.11
echo Zlib ............CONFIG
CC=gcc CXX=g++ CFLAGS=-O3 CXXFLAGS=-O3 ./configure --prefix=$nalu_install_dir/zlib/1.2.11
echo Zlib ............BUILD
make -j 16
echo Zlib ............INSTALL
make install
echo Hdf5 ............CLEAN

cd $BASE_DIR
echo Hdf5 ............CLEAN
mv hdf5-1.10.10 previous
cp tarFiles/hdf5-1.10.10.tar.gz ./
tar -zxvf hdf5-1.10.10.tar.gz
rm -rf hdf5-1.10.10.tar.gz
cd hdf5-1.10.10/
echo Hdf5 ............CONFIG
./configure CC=mpicc FC=mpif90 CXX=mpicxx CXXFLAGS="-fPIC -O3" CFLAGS="-fPIC -O3" FCFLAGS="-fPIC -O3" --enable-parallel --with-zlib=$nalu_install_dir/zlib/1.2.11 --prefix=$nalu_install_dir/hdf5/1.10.10
echo Hdf5 ............BUILD
make -j 32
echo Hdf5 ............INSTALL
make install

cd $BASE_DIR
echo Pnetcdf ............CLEAN
mv pnetcdf-1.12.2 previous
cp tarFiles/pnetcdf-1.12.2.tar.gz ./
tar -zxvf pnetcdf-1.12.2.tar.gz
rm -rf pnetcdf-1.12.2.tar.gz
cd pnetcdf-1.12.2
echo Pnetcdf ............CONFIG
./configure --prefix=$nalu_install_dir/pnetcdf/1.12.2 CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I$nalu_install_dir/pnetcdf/1.12.2/include -O3" LDFLAGS=-L$nalu_install_dir/pnetcdf/1.12.2/lib --disable-fortran
echo Pnetcdf ............BUILD
make -j 32
echo Pnetcdf ............BUILD
make install

cd $BASE_DIR
echo Netcdf ............CLEAN
mv netcdf-c-4.8.1.tar.gz previous
cp tarFiles/netcdf-c-4.8.1.tar.gz ./
tar -zxvf netcdf-c-4.8.1.tar.gz
rm -rf netcdf-c-4.8.1.tar.gz
cd netcdf-c-4.8.1
echo Netcdf ............CONFIG
./configure --prefix=$nalu_install_dir/netcdf/4.8.1 CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I$nalu_install_dir/hdf5/1.10.10/include -I$nalu_install_dir/pnetcdf/1.12.2/include -O3" CPPFLAGS=${CFLAGS} LDFLAGS="-L$nalu_install_dir/hdf5/1.10.10/lib -L$nalu_install_dir/pnetcdf/1.12.2/lib -L$nalu_install_dir/zlib/1.2.11/lib -Wl,--rpath=$nalu_install_dir/hdf5/1.10.10/lib" --enable-pnetcdf --enable-parallel-tests --enable-netcdf-4 --disable-shared --disable-fsync --disable-cdmremote --disable-dap --disable-doxygen --disable-v2
echo Netcdf ............BUILD
make
echo Netcdf ............INSTALL
make install

