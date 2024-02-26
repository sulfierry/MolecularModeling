# Formalizing the mathematical expressions for both methods in Markdown format

### Linear Interpolation Method (Assumed for `gmx sham`/`gmx wham`)

Interpolation in the context of calculating free energy landscapes often involves estimating values between known data points to construct a continuous representation of the energy landscape. While specific details of the interpolation method used in `gmx sham` or `gmx wham` are not detailed in the provided source code, a common approach is linear interpolation, which can be represented mathematically as follows:

Given two known points \((x_1, y_1)\) and \((x_2, y_2)\), the linear interpolation formula to find a value \(y\) at a point \(x\) is given by:

$$
y = y_1 + \frac{(x - x_1) \cdot (y_2 - y_1)}{x_2 - x_1}
$$

This method is straightforward but may not capture the complexities of free energy landscapes, especially in higher dimensions or with complex energy barriers.

**Function Location:** Specific interpolation functions or methodologies, like weighted histogram analysis method (WHAM), are often implemented within the source code of molecular dynamics tools but are not explicitly detailed in the provided `gmx_sham.cpp` or `gmx_wham.cpp` files.

### Kernel Density Estimation (KDE) Method

The KDE method provides a way to estimate the probability density function (PDF) of a random variable in a non-parametric way. In the context of your script for calculating the Free Energy Landscape, KDE is used to estimate the density of states, which can then be converted into free energy using the Boltzmann relation. The KDE for a set of \(n\) points \(\{x_i\}\) can be mathematically represented as:

$$
\hat{f}(x) = \frac{1}{n \cdot h} \sum_{i=1}^{n} K\left( \frac{x - x_i}{h} \right)
$$

where \(\hat{f}(x)\) is the estimated density at point \(x\), \(K\) is the kernel function (e.g., Gaussian), and \(h\) is the bandwidth, a parameter that controls the smoothness of the density estimate.

The conversion from the estimated density to free energy is typically done using the relation:

$$
G(x) = -k_B T \ln(\hat{f}(x))
$$

where \(G(x)\) is the free energy at point \(x\), \(k_B\) is the Boltzmann constant, and \(T\) is the temperature.

**Function Location:** Your script implements KDE in the `freeEnergyLandscape.py` file, using functions from libraries such as `numpy` and `scipy` for numerical operations and density estimation.

This KDE approach offers advantages in terms of smoothness and adaptability to complex data distributions, making it particularly suited for capturing detailed features of free energy landscapes.
