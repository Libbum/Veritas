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

```bash
$ cmake ..
$ make
$ make install
```

This generates a `veritas` executable in the bin directory and sets up the environment ready to start running straight away.

Acknowledging Veritas
---------------------

If you use Veritas in your research, please reference the following article:

Svedung Wettervik et al., *European Physical Journal D* **71** 157 (2017)

DOI: [10.1140/epjd/e2017-80102-2](http://dx.doi.org/10.1140/epjd/e2017-80102-2)

```
@Article{SvedungWettervik2017,
  author   = {Svedung Wettervik, Benjamin and DuBois, Timothy C. and Siminos, Evangelos and F{\"u}l{\"o}p, T{\"u}nde},
  title    = {Relativistic {V}lasov-{M}axwell modelling using finite volumes and adaptive mesh refinement},
  journal  = {The European Physical Journal D},
  year     = {2017},
  volume   = {71},
  number   = {6},
  pages    = {157},
  issn     = {1434-6079},
  abstract = {The dynamics of collisionless plasmas can be modelled by the Vlasov-Maxwell system of equations. An Eulerian approach is needed to accurately describe processes that are governed by high energy tails in the distribution function, but is of limited efficiency for high dimensional problems. The use of an adaptive mesh can reduce the scaling of the computational cost with the dimension of the problem. Here, we present a relativistic Eulerian Vlasov-Maxwell solver with block-structured adaptive mesh refinement in one spatial and one momentum dimension. The discretization of the Vlasov equation is based on a high-order finite volume method. A flux corrected transport algorithm is applied to limit spurious oscillations and ensure the physical character of the distribution function. We demonstrate a speed-up by a factor of 7 {\texttimes} in a typical scenario involving laser pulse interaction with an underdense plasma due to the use of an adaptive mesh.},
  doi      = {10.1140/epjd/e2017-80102-2},
}
```

License
-------

The project is licensed under the MIT license.

[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2FLibbum%2FVeritas.svg?type=large)](https://app.fossa.io/projects/git%2Bgithub.com%2FLibbum%2FVeritas?ref=badge_large)
