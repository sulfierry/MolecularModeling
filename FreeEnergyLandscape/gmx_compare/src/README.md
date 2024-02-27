## Kernel Density Estimation (KDE) vs. Histogram-based Methods

Kernel Density Estimation (KDE) offers several advantages over histogram-based methods, such as those used in SHAM and WHAM, for estimating the probability density function (PDF) of a continuous random variable.

### Advantages of KDE

1. **Smoothness and Continuity**: KDE produces a smooth PDF, eliminating the discontinuities of histograms.
   $$P_{KDE}(x) = \frac{1}{Nh}\sum_{i=1}^{N} K\left(\frac{x-x_i}{h}\right)$$
   - \(K\) is the kernel function, typically Gaussian.
   - \(h\) is the bandwidth, determining the smoothness.
   - \(x_i\) are the data points, and \(N\) is the number of data points.

2. **Bandwidth Selection**: KDE utilizes methods like Scott's rule to optimize bandwidth.
   $$\text{Bandwidth (Scott's Rule)} = \sigma \cdot n^{-1/5}$$
   - \(\sigma\) is the standard deviation of the sample.
   - \(n\) is the sample size.

3. **Handling Sparse Data**: KDE can effectively estimate densities in sparse data areas, a challenge for histogram-based methods.

### Normalization in KDE

Normalization ensures the integral of the KDE over all states equals one, making it a valid probability distribution:
$$\int P(x) \, dx = 1$$

### Piecewise Constant Approximation in SHAM and WHAM

Histogram-based methods approximate the PDF using a piecewise constant function.
# AQUI
```math

P_{\text{histogram}}(x) = \frac{1}{N \cdot \Delta x} \sum_{i=1}^{N} \mathbf{1}_{[x_i, x_{i+1})}(x)
```

$$P_{\text{histogram}}(x) = \frac{1}{N \cdot \Delta x} \sum_{i=1}^{N} 1_{[x_i, x_{i+1})}(x)$$


- $N$ is the total number of bins.
- $\Delta x$ is the width of each bin.
-  $1_{[x_i, x_{i+1})}(x)$ is the indicator function, where $\mathbf{1}_{[a,b)}(x)$ equals 1 if $x$ is within the interval $[a,b)$ and 0 otherwise.


### Normalization of Histograms

Normalization ensures that the area under the histogram equals one, making it a valid PDF:
$$\int_{-\infty}^{\infty} P_{histogram}(x) \, dx = 1$$

### Cost-Benefit Analysis of KDE

- **Pros**: Provides a statistically superior and more nuanced estimation of the PDF.
- **Cons**: More computationally intensive, especially for large datasets or high-dimensional spaces.

KDE's smooth and accurate PDF estimates offer valuable insights into the energy landscape that histogram-based methods may miss.

## Probability Calculation and Normalization

The probability \(P\) of a system's state occurring within a specific bin is calculated as:

$$P(bin) = \frac{N_{bin}}{N_{total}}$$

where \(N_{bin}\) is the number of occurrences in the bin, and \(N_{total}\) is the total number of observations.

### Normalization of Probability

To ensure the probabilities sum up to 1 across all bins, normalization is performed as follows:

$$P_{normalized}(bin) = \frac{P(bin)}{\sum_{bins} P(bin)}$$

This step is crucial for the accurate calculation of the free energy landscape.


## Free Energy Landscape Calculation with Python Script

This document outlines the use of Kernel Density Estimation (KDE) and Boltzmann Inversion by the `freeEnergyLandscape.py` script to estimate the free energy landscape of a system from its collective variables.

### Kernel Density Estimation (KDE) Method

KDE smooths the distribution of collective variables to enhance the estimation of probability densities. It employs Scott's rule for bandwidth selection, ensuring an appropriately smooth density estimate for the data's scale.

$$
\text{Bandwidth (Scott's Rule)} = \sigma \cdot n^{-1/5}
$$

- $\sigma$ is the standard deviation of the sample.
- $n$ is the sample size.

KDE for a set of points is given by:

$$
\hat{f}(x) = \frac{1}{n \cdot h} \sum_{i=1}^{n} K\left( \frac{x - x_i}{h} \right)
$$

- $K$ is the kernel function, typically Gaussian.
- $h$ is the bandwidth.

### Boltzmann Inversion

The free energy $G(x)$ of a state is calculated from the probability density $\hat{f}(x)$, obtained via KDE:

$$
G(x) = -k_BT \ln(\hat{f}(x))
$$

- $k_B$ is the Boltzmann constant.
- $T$ is the temperature.

### Probability Normalization

To ensure accurate energy calculation, the probability densities are normalized over the entire range of possible values, ensuring they sum to one across the configurational space:

$$
\int_{-\infty}^{\infty} \hat{f}(x) \, dx = 1
$$

### Interpolation

The smoothly interpolated probability density function provided by KDE enables calculating the free energy landscape over a continuous range of collective variable values, highlighting energetically favorable states and barriers between them.


## SHAM (Self-Consistent Histogram Analysis Method)

SHAM is a method used for data analysis in molecular dynamics simulations, especially useful for calculating free energy differences between states. The iterative nature of SHAM, which seeks to find a self-consistent set of free energies, influences its computational complexity.

### Time Complexity
- **Depends on the number of states (or histograms) analyzed**, denoted as \(N\), and the number of iterations required for convergence.
- **Time Complexity**: \(O(IN)\), assuming the work done per state per iteration is constant.
- **Quadratic Behavior**: Under certain conditions, such as when the number of iterations \(I\) approaches or exceeds \(N\), the time complexity can approach \(O(N^2)\), indicating a potential quadratic behavior in cases where convergence is slow or difficult to achieve.

### Space Complexity
- **Influenced by the storage of histograms and auxiliary data structures** needed for analysis.
- **Space Complexity**: \(O(NC)\), where \(C\) represents the constant size of each histogram.

## WHAM (Weighted Histogram Analysis Method)

WHAM is employed to calculate free energy differences from multiple histograms, combining them in a weighted manner to optimize the free energy estimate.

### Time Complexity
- **Similar to SHAM**, depends on the number of histograms (\(N\)) and the number of iterations (\(I\)) to achieve convergence.
- **Time Complexity**: \(O(IN)\), where each iteration involves detailed calculations over each histogram.
- **Quadratic Behavior**: Like SHAM, WHAM can exhibit quadratic time complexity (\(O(N^2)\)) under conditions where the number of iterations is large relative to the number of histograms, particularly when convergence criteria are stringent or the dataset is complex.

### Space Complexity
- **Affected by the storage of multiple histograms and the associated weights**.
- **Space Complexity**: \(O(NC)\), with \(C\) representing the constant size of each histogram, similarly to SHAM.


## Python Script for Free Energy Landscape Analysis

This Python script utilizes Kernel Density Estimation (KDE) to analyze free energy landscapes from molecular dynamics simulations. The computational complexity of the script is influenced by the data size and the nature of the KDE implementation.

### Time Complexity
- **Efficient KDE Implementation**: The script's time complexity is primarily determined by the KDE computation. Leveraging the `scipy` library's KDE, known for efficient handling of large datasets, the average case time complexity is generally \(O(n \log n)\), where \(n\) is the number of data points.
- **Parallel Processing Enhancements**: Given that the script is designed to utilize parallel processing capabilities, this can further reduce the effective computational time, especially on multi-core systems, making the actual run time lower than the theoretical complexity might suggest.
- **High-Dimensional Data Caution**: While the KDE computation avoids naive pairwise distance calculations, the complexity can still be higher for very large or high-dimensional datasets due to the intrinsic challenges of high-dimensional space (the "curse of dimensionality"). However, this does not necessarily imply quadratic complexity in all scenarios.


### Space Complexity
- **Data and Computation Requirements**: The space complexity is influenced by the requirements for storing the dataset and any intermediate structures needed for KDE computation. This is typically \(O(n)\), where \(n\) is the number of data points, with additional considerations for the storage of KDE results.
- **Optimized Data Handling**: Assuming the implementation makes use of efficient data structures and algorithms, the space complexity remains manageable even for large datasets.

### Performance Consideration
The actual performance of this script can vary significantly depending on the specifics of the `scipy` library's KDE implementation and the dataset characteristics. For datasets that are large or high-dimensional, employing strategies to optimize data handling and computation, like using more efficient algorithms or data structures for KDE, can significantly enhance performance.



