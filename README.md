# mcd-recurse

The recursor is part of our toolkit to tackle the 2020 Computational Geometry:
Solving Hard Optimization Problems (CG:SHOP) challenge of finding a minimal
convex decomposition of a pointset.

# Obtaining the source code

Clone the git repository:
The logging class is included via a git submodule, so clone the source using

    git clone  https://github.com/cgalab/mcd-recurse &&
    (cd mcd-recurse && git submodule update --init src/easyloggingpp)

You will also need [Shewchuk's triangle][triangle], which you can download from
http://www.netlib.org/voronoi/triangle.zip.  Then unpack it in `src/triangle/`.

    cd mcd-recurse/src/triangle &&
    wget http://www.netlib.org/voronoi/triangle.zip &&
    sha256sum triangle.zip | grep 1766327add038495fa3499e9b7cc642179229750f7201b94f8e1b7bee76f8480 && unzip triangle.zip &&
    cd ../..

[triangle]: https://www.cs.cmu.edu/~quake/triangle.html

# Building

To build the tool, run cmake and make:

    mkdir build &&
    cd build &&
    cmake -DCMAKE_BUILD_TYPE=Release .. &&
    make

# Running the tool

`mcd-recurse` accepts a number of optional parameteres, see the `--help` output
for a brief overview.  It needs a pointset as input (first argument or via
stdin), which by default is simply a series of vertex coordinates, `x y`, one
per line.

# License

`mcd-recurse` is free software.  You may redistribute it and/or modify
it under the terms of the GNU General Public License (v3).
