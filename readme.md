## Building
First clone HElib with git.
```shell
git clone git@github.com:homenc/HElib.git
cd HElib
git checkout 3e337a6
```
Then, copy the `bgvboot-final.patch` file into `HElib/` and apply the changes by
```shell
git apply bgvboot-final.patch
```
You need to replace the path `your_path_to_main_dir` in `HElib/src/extractDigits.cpp` to the directory path of `fatboot.cpp`.

Then build and install HElib following the instructions in [INSTALL.md](https://github.com/homenc/HElib/blob/master/INSTALL.md), some dependencies like NTL and gmp may be needed. Remember to build in either `Release` or `RelWithDebInfo` mode.

Finally, switch to the directory of `fatboot.cpp` and build the test program by
```shell
cmake -S . -B build
cmake --build build
```
If CMake cannot final HElib, you may need to change the `helib_DIR` variable in `CMakeLists.txt` to the correct install path of HElib.

## Running
If the building is successful, it will produce a binary file `build/fatboot`.
Running `build/fatboot -h` will display the following help message
```
Usage: build/fatboot [i=<arg>] [h=<arg>] [t=<arg>] [newbts=<arg>] [newks=<arg>] [thick=<arg>] [repeat=<arg>]
  i      index of the chosen parameter set [ default=0 ]
  h      hwt of encapsulated key [ default=0 ]
  t      parameter t used in new bts [ default=0 ]
  newbts if new bts is used [ default=0 ]
  newks  if new ks is used [ default=1 ]
  thick  if thick bts is used [ default=0 ]
  repeat number of tests [ default=5 ]
```
+ Five parameter sets are available (given in Table 3 of the paper with ID I, II, III, IV and V). By setting the argument `i` to an integer in 0 to 4, the corresponding parameter set is chosen.
+ The Hamming weight of the encapsulated secret key needs to be assigned through the argument `h`. Again, please refer to Table 3 for the values of `h` with 80 or 128 bits of security.
+ Argument `t` controls the $t$ in our paper, i.e., $t=v_p(\Delta)$. Setting `t=0` tells the program to decide a positive $t$ that is as small as possible (type-A parameters), while setting `t=-1` uses $\Delta=\Delta_0$, i.e., $t=0$ (type-B parameters). Setting `t` to any positive integer will also set $t$ to the same value.
+ Set `newbts=1` if you want to test our optimized digit removal. Set `newbts=0` to test the native implementation of HElib.
+ Set `thick=1` to test general bootstrapping. Set `thick=0` to test thin bootstrapping.
+ Argument `repeat` indicates the number of tests to run. The performance data will be averaged over these tests. 