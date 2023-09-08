Initial README.md for a Computational Electromagnetics (CEM), numerical physics project.

* The project will explore electromagnetics in context of rough surface scattering.

* Mainly based on the following material:
    * I. Simonsen, _Optics of Surface Disordered Systems: A Random Walk Through Rough Surface Scattering Phenomena_. Eur. Phys. J. Special Topics 181, 1 (2010). http://web.phys.ntnu.no/~ingves/Science/Publications/Abstracts/Published/Abst058.php

    * I. Simonsen, and A. A. Maradudin, _Numerical simulation of electromagnetic wave scattering from planar dielectric films deposited on rough perfectly conducting substrates_. Opt. Commun. 162, 99 (1999). http://web.phys.ntnu.no/~ingves/Science/Publications/Abstracts/Published/Abst002.php

### Theory
* Previous studies has used a varying the relative EM permittivity, $\varepsilon$. In this project we will vary the relative EM permeability, $\mu$.

* The central equations are derived from Maxwells equations.

* The reduced Rayleigh condition is used. 

* The theory may also touch on the Bloch theorem for periodic surfaces.

### Code

* The focus of the code will be high performance computing for the resulting matrix system of integral equations.

* Prototyping in `Julia` or `Python`. Look at the possibility of porting to `C` or `C++` for better scaleability and parallelization.