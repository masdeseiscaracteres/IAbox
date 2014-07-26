IAbox
========
What is IAbox?
--------------
The Interference Alignment MATLAB/Octave Toolbox (IAbox) is a set of 
MATLAB/Octave functions implementing results and algorithms that have 
emerged from studies on the interference aligment (IA) idea and other
related literature.

Other useful functionalities, not necessarily from the IA literature, but 
from the wireless communications world, have also been implemented.

This is a work in progress. It may (and will) contain errors, bugs, misconceptions and even some working stuff.
Critics, bug reports, collaborations and any sort of feedback is welcome.

Implemented Algorithms
--------------------------------------------
[Listed by publication date](src/Algorithms/README.md)

Implemented DoF / Feasibility Results
-----------------------------------------------------
[Listed by publication date](src/DoFandFeasibility/README.md)

Naming / style conventions
--------------------------
The code in this toolbox tries to stick to the MATLAB Programming Style 
Guide Wiki that can be found at
[https://sites.google.com/site/matlabstyleguidelines/home][1]

Input arguments to algorithms/functions are passed in the following order

`nT,nR,H,U,V,D,W`

Contributors
------------

List of people who have contributed to this toolbox:

 - Óscar González
 - Christian Lameiro
 - Jacobo Fanjul
 
License
-------
This repository is maintained by [Óscar González][2]. It is distributed under the terms of the BSD (3-Clause) License.  In short, this means that everyone is free to use it, to modify it and to redistribute it on a free basis. It is not in the public domain; it is copyrighted and there are restrictions on its distribution (see [LICENSE.txt](LICENSE.txt)).


  [1]: https://sites.google.com/site/matlabstyleguidelines/home
  [2]: http://gtas.unican.es/people/oscargf