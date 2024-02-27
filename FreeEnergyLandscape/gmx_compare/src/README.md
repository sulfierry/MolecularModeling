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
$$P_{histogram}(x) = \frac{1}{N \cdot \Delta x}\sum_{i=1}^{N} \mathbf{1}_{[x_i, x_{i+1})}(x)$$
- \(N\) is the total number of bins.
- \(\Delta x\) is the width of each bin.
- \(\mathbf{1}_{[x_i, x_{i+1})}(x)\) is the indicator function.

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

## Free Energy Calculation

The free energy \(G\) for each bin, after normalization, is calculated using the Boltzmann relation:

$$G(bin) = -kT \ln(P_{normalized}(bin))$$

Here, \(k\) is the Boltzmann constant, \(T\) is the temperature, and \(P_{normalized}(bin)\) is the normalized probability. The free energy landscape is adjusted so that the minimum free energy value is set to zero for ease of interpretation.

