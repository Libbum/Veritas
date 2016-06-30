<img src="http://veritas.readthedocs.io/en/latest/_images/veritas.svg" width=300px />

**V**lasov **E**ule**RI**an **T**ool for **A**cceleration **S**tudies

A relativistic Vlasov-Maxwell solver implemented using a fourth order finite volume scheme on a block structured adaptive mesh; upwind biased WENO-type reconstruction of face values; a Flux Corrected Transport algorithm to enhance stability and Runge-Kutta time advancement.

Further information can be found on the [Veritas Web Page](http://ft.nephy.chalmers.se/veritas/).

Documentation
-------------

A quick start guide can be found over at [Read The Docs](http://veritas.readthedocs.io/).

Installation
------------

If you don't care about specifics, from the build directory run 

    $ cmake ..
    $ make
    $ make install
  
This generates a `veritas` executable in the bin directory and sets up the environment ready to start running straight away.

License
-------

The project is licensed under the MIT license.
