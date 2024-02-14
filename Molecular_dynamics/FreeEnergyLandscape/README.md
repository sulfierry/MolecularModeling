# Free Energy Landscape Analysis

This Python class, `FreeEnergyLandscape`, is designed to analyze and visualize the free energy landscape of molecular systems using collective variables (CVs), such as angles and distances, obtained from molecular dynamics simulations. The approach leverages the concept of Boltzmann inversion to transform probability distributions into free energy surfaces, providing insights into the thermodynamic and kinetic properties of the system.

## Theoretical Background

The free energy landscape is a conceptual and computational tool used to understand the energetics and dynamics of molecular systems. It is particularly useful in the study of complex processes such as protein folding, chemical reactions, and phase transitions.

### Collective Variables (CVs)

Collective Variables (CVs) are a set of coordinates that describe the macroscopic state of a system. They are used to reduce the complexity of molecular systems by focusing on the relevant degrees of freedom. Examples include the distance between two atoms, angles, dihedrals, and more complex descriptors.

### Boltzmann Inversion

The Boltzmann inversion method is used to calculate the free energy landscape from the probability distribution of CVs. The free energy, \(G\), of a state is related to its probability, \(P\), by the Boltzmann equation:

\[G = -k_B T \ln(P)\]

where \(k_B\) is the Boltzmann constant (\(1.380649 \times 10^{-23} J/K\)), \(T\) is the temperature, and \(P\) is the probability density of the CVs. This relationship allows us to convert a histogram of CV values into a free energy surface.

## Implementation Details

### Class `FreeEnergyLandscape`

- **Initialization**: Sets up the paths to data files for CV1 and CV2, the temperature, and the Boltzmann constant. A custom color map is also defined for visualizing the energy levels.

- **Data Loading**: Reads the CV data from the specified files, assuming that the data of interest is in the second column.

- **Boltzmann Inversion**: For a given set of CV data, calculates the histogram (probability distribution), converts it to free energy using the Boltzmann equation, and visualizes the result.

- **Energy Landscape Visualization**: Combines the CV data to construct a 2D free energy landscape. This involves estimating the joint probability density function using a Gaussian kernel density estimate, converting this density into free energy, and plotting the landscape with contour lines representing different energy levels.

### Usage

To use this class, instantiate it with paths to the CV1 and CV2 data files, then call the `main` method. This will load the data, perform the Boltzmann inversion for each CV, and plot the 2D free energy landscape.

## Example

```python
if __name__ == "__main__":
    cv1_path = './proj1Out.txt'
    cv2_path = './proj2Out.txt'
    fel = FreeEnergyLandscape(cv1_path, cv2_path)
    fel.main()
