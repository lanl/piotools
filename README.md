# PIO Tools (LANL C21102)
## Author: Sriram Swaminarayan, sriram@lanl.gov
## See [Copyright.txt](Copyright.txt) for copyright and licensing information

This is a set of utilities written in Python and C++ for reading and
manipulating PIO files. These files are intended to be building blocks
for other scripts / programs that will use them to do great things.
The PIO format itself is described in LA-UR-05-7425 at page 102 and
embodied in code LA-CC-05-052.  In addition to reading the raw PIO
files, the bundled utilities will expand variable that use a
compressed sparse row notation.

## C++ Files
* `PIO`: The base PIO class contained in header file `pio.cpp` that
	enables you to load meta data contained in a raw PIO file and read
	variables from it.  If you want to conserve memory, you will want
	to write your code based on this class and create your own
	variable readers based on the code in class `PioInterface`.
* `PioInterface`: This class, contained in files `pioInterface.hpp`
      and `pioInterface.cpp`, provides a nicer interface to class
      `PIO` with utilities that will read in cell variables and expand
      compressed variables
	  
## Python Files
* `pio.py`: This file is an amalgamation of the code contained in the
   C++ classes `PIO` and `PioInterface`.  It enables you to read a PIO
   file, query variables from it, and also write out a new PIO file
   with a new set of variables.  As a demo, If you run it on its own
   with a PIO file as an argument, it will create a file
   `bigfile-dmp000000` that includes all variables from the original
   file and adds to them all the compressed variables expanded to be
   cell arrays so that the ParaView PIO reader can pull them in.


## Examples directory
* `addVar.py`: Demonstrates the use of the Python `pio` class by
  adding the processor ID as an additional variable to the file.  It
  also creates a `tmp.pio` to ease opening the new file with
  ParaView.
* `testPio.cpp`: A simple program to show how to use the C++ `PIO`
  class to read in a variable from the file.

