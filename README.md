## NAME
##  M_blas -- BLAS routines in a Fortran module format

The BLAS (Basic Linear Algebra Subprograms) are well-established routines
that provide standard building blocks for performing basic vector and
matrix operations.

There are several versions of the library available; this is based off
the Public Domain
[FORTRAN reference implementation of BLAS Basic Linear Algebra Subprograms](http://www.netlib.org/blas/)
provided by the National Science Foundation.

It has been converted to a free-format Fortran module and is in the
process of being changed to use more recent Fortran features primarily
for investigating the significance of what benefits this can have, such
as an __automatic calling interface__ and __name collision avoidance__
provided by being modularized, __INTENT specifications__, __array syntax__
and so on.

For compatibility purposes the wrapper file src/compatible.f90 is
included, which allows for an old non-interface calling style so if
you are using a BLAS library compatible with the reference version
you should be able to use this library (and the accompanying module)
as a drop-in replacement at the cost of loosing the strict interface
the module provides and the additional call overhead to the wrapper.

If that works with your existing code, the next step is to then include
USE statements in the calling source and to remove src/compatible.f90
and rebuild.

If there are significant benefits the idea is to then proceed with
a similar transformation of LAPACK, but I am just beginning this and
unless someone else is excited about it and wants to join in this is a
rainy-day project at this point, and __should not be considered stable__.

I would love to hear from anyone who tries this as a replacement for
their current libblas.a or libregblas.a.

## DOWNLOAD AND BUILD

### gmake ![gmake](docs/images/gnu.gif)

   ```bash
       git clone https://github.com/urbanjost/M_blas.git
       cd M_blas/src
       # change Makefile compiler information as appropriate in Makefile
       make # for gfortran
   ```
   This will compile the M_blas module and test programs.

### fpm ![fpm](docs/images/fpm_logo.gif)

   Alternatively, download the github repository and build it with
   fpm ( as described at [Fortran Package Manager](https://github.com/fortran-lang/fpm) )

   ```bash
        git clone https://github.com/urbanjost/M_blas.git
        cd M_blas
        # currently the tests are run from a script because they need redirection
        # and I want some timing information and to test multiple compilers. You
        # will have to make your own run script or run the tests manually if you
        # do not have bash(1).
        #
        # Take out or add the compilers from run.sh you wish to run then execute
        # the run.sh(1) script, twice if you just want run-times:
        bash run.sh
        # or just run "fpm build".
   ```

   or just list it as a dependency in your fpm.toml project file.

```toml
        [dependencies]
        M_blas        = { git = "https://github.com/urbanjost/M_blas.git" }
```
## References

- There is an `fpm` version of the
  [reference BLAS library]( https://github.com/brocolis/blas.git)
  available.

- See the homepage for [BLAS](http://www.netlib.org/blas)
  for more detailed information.
---
## [CHANGELOG](CHANGELOG.md)
---
