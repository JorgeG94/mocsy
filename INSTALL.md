# How to build 

Mocsy has several build systems. 

## CMake (preferred, default) 

Ensure you have a functioning Fortran compiler and a CMake of version at least 3.22, then simply do:
```
mkdir build
cd build
cmake ..
make -j 
```

To specifiy where to install mocsy use: `cmake -DCMAKE_INSTALL_PREFIX=$MOCSY_PREFIX ../` and use `make -j install`

## Make 

Simply ensure you have a Fortran compiler and use `make libmocsy.a` 


## FPM 

Download and install the [Fortran Package Manager](https://fpm.fortran-lang.org/). You can install the latest version by doing:

```
git clone https://github.com/fortran-lang/fpm
cd fpm
./install.sh
```

This will put the `fpm` executable in your `~/.local/bin` if you don't have this in your `$PATH` simply do: `export PATH=$PATH:$HOME/.local/bin`

To test the mocsy build use: `fpm test --profile release`

To install mocsy to a defined location use: `fpm install --prefix $MOCSY_PREFIX --profile release`

*NOTE*: The FPM install only supports using double precision for mocsy. 





