!#/bin/sh

#este script irá baixar, configurar e instalar as bibliotecas
# Desenvolvedor: Vitoria regia contato: vitoria.regia.016@ufrn.edu.br
#data de criação: 23/08/2023
#Dta de modificação

# ls - list
#pwd - observar no diretorio
#cd .. -voltar diretoprio so um 
#cd ~ 


#variavel espaço na memoria da maquina 
dir=/home/aluno/cac3002/cdo

#-dowload da biblioteca ZLIB
#wget https://www.zlib.net/fossils/zlib-1.2.8.tar.gz
#cd zlib-1.2.8/
#./configure --prefix=/home/aluno/cac3002/cdo    #--prefiz=$dir
#cd ..
#make
#make chek
#make install
#cd $dir


#-dowload da biblioteca HDF
#wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.13/bin/linux-x86_64/hdf5-1.8.13-linux-x86_64-shared.tar.gz -O hdf5-1.8.13.tar.gz
#tar -xvf hdf5-1.8.13.tar.gz
#cd hdf5-1.8.13
#./configure -with-zlib=$dir --prefix=$dir CFLAGS=-fPIC
#make && make check && make install
#make all
#make 
#make install

#-download da biblioteca Netcdf
#wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.5.0.tar.gz -O netcdf-4.5.0.tar.gz

#-download da biblioteca jasper
#wget http://www.ece.uvic.ca/~frodo/jasper/software/jasper-1.900.0.zip

#-download da biblioteca grip_api
#wget https://src.fedoraproject.org/lookaside/pkgs/grib_api/grib_api-1.24.0-Source.tar.gz/sha512/11d6992714880b5855224e706a71921c25ffffa154892b9231bd4f21dec175d6c2b3a7b921864e339e50fc3eb3a9d4744bb506b8b8ed663b2fb6d0687e200649/grib_api-1.24.0-Source.tar.gz -O grib_api-1.24.0.tar.gz

#-download da biblioteca CDO
#wget https://code.mpimet.mpg.de/attachments/download/15653/cdo-1.9.1.tar.gz
