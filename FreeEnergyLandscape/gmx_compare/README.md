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

# GMX SHAM Mechanism Description

The `gmx_sham` tool is part of the GROMACS suite, designed to analyze free energy landscapes from simulation data. The core functionality revolves around creating multidimensional histograms from the provided data and calculating the free energy landscape. Here's a detailed breakdown of the process, including the improved explanation of "accumulated probability of all data points":

## Data Preparation

1. **Extremes Determination**: For each eigenvector, the code calculates minimum and maximum values across all data points, adding a small buffer to ensure all data is encompassed (`find_extremes` function).

2. **Normalization and Volume Correction**: Data undergo normalization and volume correction if necessary, based on user-specified dimensions. This step adjusts the probability calculations for spatial dimensions, correcting for volume effects that increase with distance in 2D and 3D spaces, as seen in the `correct_for_volume_effects` method.

## Histogram Creation and Free Energy Calculation

1. **Binning**: Data is sorted into multidimensional bins, each representing specific intervals of the eigenvectors (`binning_data` function). Points are allocated to bins based on their eigenvector values.

2. **Accumulated Probability Calculation**: The accumulated probability \(P_{acum}(bin)\) for each bin is calculated as the number of data points within the bin (\(N_{bin}\)) divided by the total number of data points (\(N_{total}\)), reflecting the proportion of occurrences for each bin relative to the entire dataset.

   \[ P_{acum}(bin) = \frac{N_{bin}}{N_{total}} \]

3. **Free Energy Calculation**: Free energy \(G\) for each bin is determined using the inverse Boltzmann relation, where the accumulated probability informs the energy level:

   \[ G(bin) = -kT \ln(P_{acum}(bin)) \]

   Here, \(k\) is Boltzmann's constant, and \(T\) is the system's temperature. This step is executed within the `calculate_free_energy` method.

4. **Probability Normalization and Energy Adjustment**: Finally, the probability in each bin is normalized, and the free energy values are adjusted so the minimum free energy across the landscape is set to zero, facilitating easier interpretation and visualization.

## Analogy for Understanding Accumulated Probability

Consider throwing darts at a target with different zones, where each zone represents a bin. Each dart hit corresponds to a data point in a bin. After many throws, comparing the number of darts per zone to the total throws gives the probability of hitting each zone. This scenario mirrors the accumulation of data points in bins and their conversion into a free energy landscape, offering insights into the system's stability and energetically favorable states.

## Conclusion

This mechanism enables `gmx_sham` to provide a detailed analysis of free energy landscapes, crucial for understanding molecular dynamics simulations. The tool's ability to dissect complex data into understandable energy landscapes makes it invaluable for researchers and scientists in the field of computational chemistry and molecular dynamics.


