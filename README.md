Interference alignment feasibility test
=========================================
What you can find here is a Matlab/Octave implementation of a feasibility test for interference alignment in MIMO interference channels with constants coefficients. It is a direct implementation of the floating point test proposed in the paper

> Ó. González, C. Beltrán, I. Santamaría, "*A Feasibility Test for Linear Interference Alignment in MIMO Channels with Constant Coefficients*," IEEE Transactions on Information Theory, to be published, 2014.

If you do not feel like downloading the code and want to see how it works in a simple scenario, you are invited to try our [online demo][1]. Note that our demo is executed in a backend server running Octave and, due to server limitations, it is constrained to a total number of 450 antennas and 100 streams. Besides, depending on the server workload these limits can sometimes be more restrictive. If you want to check the feasibility of larger systems you will need the code herein.

The code has been developed and copyrighted © 2014 by [Óscar González][2]. It is distributed under the terms of the BSD (3-Clause) License.  In short, this means that everyone is free to use it, to modify it and to redistribute it on a free basis. It is not in the public domain; it is copyrighted and there are restrictions on its distribution (see LICENSE.txt).

We welcome contributions, bug reports and feedback. Let us know what you think! If you use the code for your research, please cite the [paper][3]

    @article{GBS2014,
    title={A Feasibility Test for Linear Interference Alignment in MIMO Channels with Constant Coefficients},
    author={Óscar González and Carlos Beltrán and Ignacio Santamaría},
    journal={IEEE Transactions on Information Theory, to be published},
    year={2014},
    }

----------

Brief description
-------
This software contains the following Matlab/Octave files:

 - [`check_IA_feasibility.m`](check_IA_feasibility.m) is an independent function implementing the aforementioned floating point feasibility test. It is the only function you will need if you want to execute the test from your own code.
 - [`system_str2vec.m`](system_str2vec.m) is a convenience function that takes a text string representing an interference channel and converts it to a set of vectors which will serve as inputs of `check_IA_feasibility`.
 - [`main.m`](main.m) is a self-explanatory script which shows the use of both  `check_IA_feasibility` and `system_str2vec` by means of a simple example.
 


  [1]: http://gtas.unican.es/IAtest
  [2]: http://gtas.unican.es/people/oscargf
  [3]: http://gtas.unican.es/pub/336