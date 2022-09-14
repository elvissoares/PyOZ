# PyOZ
A python implementation of a solver to the Ornstein–Zernike equation

In the integral theory of fluid [^1], the Ornstein–Zernike equation [^2] connects the radial distribution function with the direct correlation function through that 

$$ h(\boldsymbol{r}) = c(\boldsymbol{r}) + \rho \int d\boldsymbol{s}\ c(\boldsymbol{r}-\boldsymbol{s}) h(\boldsymbol{s})$$

where $h(r) = g(r)-1$ is the total correlation function, $g(r)$ is the radial distribution function, and $c(r)$ is the direct correlation function.

To solve this equation we need a closure relation. We can define a auxiliary quantity $\gamma(r)=h(r)-c(r)$, such that

$$c(r) = e^{-\beta u(r)+\gamma(r)+b(r)}-\gamma(r)-1$$

and $b(r)$ is the bridge-function. The following closure relations gives different bridge functions:
- [x] **P**ercus-**Y**evick (**PY**) - $b(r) = \log(1+\gamma(r))-\gamma(r)$
- [x] **H**iper**N**etted **C**hain (**HNC**) - $b(r) = 0$
- [ ] **M**ean-**S**pherical **A**pproximation (**MSA**) - 

# Examples

On the folder 'examples' you can find different applications of the OZ solver. 

## Hard-Sphere Fluids
|![Figure1](https://github.com/elvissoares/PyOZ/blob/main/examples/radialdistributionfunction-hardspheres.png)|![Figure2](https://github.com/elvissoares/PyOZ/blob/main/examples/contactvalue-rdf-hardspheres.png)|
|:--:|:--:|
| <b>Fig.1 - The radial distribution function of a pure hard-sphere fluid for three different densities. The symbols represent MC data. </b>| <b>Fig.2 - The contact value of the radial distribution function of a pure hard-sphere fluid as a function of the bulk density. The symbols represent MC data. </b>|

# References
[^1]: Hansen, Jean-Pierre, and Ian Ranald McDonald. *Theory of simple liquids: with applications to soft matter.* Academic press, 2013.

[^2]: Ornstein, L.S.; Zernike, F. [*Accidental deviations of density and opalescence at the critical point of a single substance.*](https://web.archive.org/web/20210206222100if_/https://www.dwc.knaw.nl/DL/publications/PU00012727.pdf) Proceedings of the Royal Netherlands Academy of Arts and Sciences. 17 (1914): 793.