# Free Energy Landscape Analysis

This Python class, `FreeEnergyLandscape`, is designed to analyze and visualize the free energy landscape of molecular systems using collective variables (CVs), such as angles and distances, obtained from molecular dynamics simulations. The approach leverages the concept of Boltzmann inversion to transform probability distributions into free energy surfaces, providing insights into the thermodynamic and kinetic properties of the system.

## Theoretical Background

The free energy landscape is a conceptual and computational tool used to understand the energetics and dynamics of molecular systems. It is particularly useful in the study of complex processes such as protein folding, chemical reactions, and phase transitions.

### Collective Variables (CVs)

Collective Variables (CVs) are a set of coordinates that describe the macroscopic state of a system. They are used to reduce the complexity of molecular systems by focusing on the relevant degrees of freedom. Examples include the distance between two atoms, angles, dihedrals, and more complex descriptors.

### Boltzmann Inversion

The Boltzmann inversion method is used to calculate the free energy landscape from the probability distribution of CVs. The free energy, $\( G \)$, of a state is related to its probability, $\( P \)$, by the Boltzmann equation:

$$\Delta G = -k_B T [\ln(P) - \ln(P_{\text{max}})] \$$

where:

- $k_B$ is the Boltzmann constant $\( 8.314 \times 10^{-3} \ \text{KJ/mol.K} \)$
- $T$ is the temperature
- $P$ is the probability distribution of the data, obtained from a histogram of molecular dynamics (MD) data
- $P_{\text{max}}$ is the maximum value of this probability distribution, representing the most likely state or the minimum of free energy


This relationship allows us to convert a histogram of CV values into a free energy surface and this formulation also adjusts the free energy values so that the minimum energy associated with the highest probability $P_{\text{max}}$ is set to 0. Because the value of $\Delta G$ for the state with $P_{\text{max}}$ will be 0 because $\ln P_{\text{max}} - \ln P_{\text{max}} = 0\$.

This approach is useful for highlighting the relative differences in free energy between different states or conformations in an MD simulation, facilitating the identification of free energy minima and the relative comparison between different states. By setting the minimum of free energy to 0, you create a clear reference point to evaluate the relative stability of other states compared to the most stable state.


### Gaussian Kernel Density Estimation (KDE)

This statistical method is crucial for estimating the probability density function (PDF) of a dataset without assuming any predefined distribution shape. In the context of molecular dynamics and simulations, it allows us to visualize and analyze the distribution of free energy across different states or configurations defined by collective variables (CVs).

The `gaussian_kde` leverages a Gaussian (normal) distribution, placing it at each point in the dataset and summing these distributions to approximate the overall data's PDF. This technique is adept at capturing the underlying structure of the data, providing a smooth, continuous representation of the free energy landscape. The smoothness of the KDE is controlled by the bandwidth parameter, which determines the width of the Gaussian kernels used.

Mathematically, the density estimation at a point $`x`$ is calculated as follows:

$$\ f(x) = \frac{1}{n \cdot h} \sum_{i=1}^{n} K\left(\frac{x - x_i}{h}\right) \$$

where:
- $`n`$ is the number of data points
- $`h`$ is the bandwidth
- $`x_i`$ are the data points
- $`K(u)`$ represents the Gaussian kernel function


$$ K(u) = \frac{1}{\sqrt{2\pi}} e^{-\frac{1}{2}u^2} $$

where:
- $u^2$ is the square of this normalized distance, which serves to weight the contribution of each data point $x_i$ to the density estimate at $x$, based on how far $x_i$ is from $x$, adjusted by the bandwidth $h$.
- $e^{-\frac{1}{2}u^2}$ decreases rapidly as $u$ increases, meaning that points further away from $x$ will have less influence on the density estimate at $x$.
- $frac{1}{\sqrt{2\pi}}$ is a normalization term that ensures the Gaussian kernel function integrates to $1$, keeping it as a valid probability distribution.

In the context of kernel density estimation: 

$$ u = \frac{x - x_i}{h} $$

where:
- $u$ is the normalized distance
- $x$ is where the density is being estimated
- $x_i$ are the data points
- $h$ is the bandwidth
  
This variable $u$ is utilized within the Gaussian kernel function $K$, which is a probability density function.


Thus, $u$ is crucial in determining how the distance between $x$ and the data points $x_i$ affects the density estimate at $x$, with the bandwidth $h$ controlling the sensitivity of this influence.

### Application in the Script

The `FreeEnergyLandscape` class employs `gaussian_kde` in several key areas:

1. **Free Energy Distribution Estimation:** By applying `gaussian_kde` to the CVs collected from molecular simulations, we obtain a continuous estimate of the free energy landscape. This estimated landscape is crucial for identifying stable configurations (minima) and understanding the transition pathways (energy barriers) between different molecular states.

2. **Visualization:** The KDE result is used to generate visual representations of the free energy landscape, including both static 3D plots and animated GIFs. These visualizations allow for an intuitive exploration of the energy landscape, facilitating the identification of significant energy features that influence molecular behavior.

### Importance in Molecular Studies

The Gaussian KDE method provides a sophisticated approach to model the complex free energy landscapes encountered in molecular dynamics studies. It enables researchers to visualize the distribution of energy states without the constraints of parametric models, offering insights into molecular stability, transitions, and the energetics of molecular interactions. By incorporating `gaussian_kde` into our analysis, we enhance our ability to decipher the intricate energy landscapes that govern molecular systems, contributing significantly to the fields of computational chemistry and biophysics.


## Implementation Details

### Class `FreeEnergyLandscape`

- **Initialization**: Sets up the paths to data files for CV1 and CV2, the temperature, and the Boltzmann constant. A custom color map is also defined for visualizing the energy levels.

- **Data Loading**: Reads the CV data from the specified files, assuming that the data of interest is in the second column.

- **Boltzmann Inversion**: For a given set of CV data, calculates the histogram (probability distribution), converts it to free energy using the Boltzmann equation, and visualizes the result.

- **Energy Landscape Visualization**: Combines the CV data to construct a 2D free energy landscape. This involves estimating the joint probability density function using a Gaussian kernel density estimate, converting this density into free energy, and plotting the landscape with contour lines representing different energy levels.

## Required Libraries

This project relies on several key Python libraries for numerical computations, image processing, and plotting capabilities. Ensure you have the following libraries installed along with their specified versions to guarantee compatibility and proper functionality of the scripts:

- **NumPy** (1.23.5): A fundamental package for scientific computing with Python, providing support for large, multi-dimensional arrays and matrices, along with a collection of mathematical functions to operate on these arrays.
- **ImageIO** (2.34.0): A library for reading and writing a wide range of image, video, scientific, and volumetric data formats. It is used in this project for creating and manipulating images and GIFs.
- **Matplotlib** (3.7.4): A comprehensive library for creating static, animated, and interactive visualizations in Python. It is used for plotting the free energy landscapes.
- **SciPy** (1.10.1): An open-source Python library used for scientific computing and technical computing. It contains modules for optimization, linear algebra, integration, interpolation, special functions, FFT, signal and image processing, and more.

To install these libraries, you can use the following command:

```bash
pip install numpy==1.23.5 imageio==2.34.0 matplotlib==3.7.4 scipy==1.10.1
```
These dependencies also are listed in the requirements.txt file. To install them, run the following command in your terminal:

```bash
pip install -r requirements.txt
```

### Usage

To use this class, instantiate it with paths to the CV1 and CV2 data files, then call the `main` method. This will load the data, perform the Boltzmann inversion for each CV, and plot the 2D free energy landscape.

## Example

```bash
python freeEnergyLandscape.py proj1Out.txt proj2Out.txt 
```

## Visualizations Generated

The `FreeEnergyLandscape` class generates three key visualizations to aid in the analysis of the molecular system's free energy landscape. Each figure provides unique insights into the system's thermodynamic and kinetic behavior.

### Figure 1: Free Energy as a Function of CV1 (Angle)

This figure displays the free energy landscape as a function of the first collective variable (CV1), which represents an angle in the molecular system. It is generated using the Boltzmann inversion method from the histogram of CV1 values.

![Alt text da image](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/fel_angle.png)

- **X-Axis**: Represents the angle (CV1) in degrees or radians, depending on the system being studied.
- **Y-Axis**: Shows the free energy in kJ/mol, calculated from the probability distribution of CV1 using the Boltzmann equation.
- **Visualization**: A plot of free energy versus CV1, highlighting the energy barriers and minima that correspond to different conformational states of the molecule.

### Figure 2: Free Energy as a Function of CV2 (Distance)

Similar to Figure 1, this visualization plots the free energy landscape but as a function of the second collective variable (CV2), typically representing a distance within the molecular system.

![Alt text da image](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/fel_distance.png)

- **X-Axis**: Denotes the distance (CV2) in appropriate units (e.g., Ångströms).
- **Y-Axis**: Indicates the free energy, with the same calculation method and units as Figure 1.
- **Visualization**: Demonstrates how the free energy changes with variations in CV2, providing insights into the significance of certain distances for the system's stability and transitions.

### Figure 3: 2D Free Energy Landscape (CV1 vs. CV2)

The third figure combines both CV1 and CV2 to produce a two-dimensional free energy landscape, offering a comprehensive view of the system's energetics over the considered collective variables.

![Alt text da image](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/fel_angle_distance.png)

- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 1.
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 2.
- **Visualization**: A contour plot or a heatmap showing the free energy levels across the CV1 and CV2 space. The color gradient represents different energy levels, with cooler colors indicating low-energy regions (minima) and warmer colors highlighting high-energy barriers. This visualization is crucial for identifying transition states, stable conformations, and understanding the molecular system's behavior under various conditions.

### Figure 4: 3D Free Energy Landscape

This figure provides a three-dimensional visualization of the free energy landscape, combining both collective variables, CV1 and CV2, along with the calculated free energy values to offer a dynamic perspective on the system's energetics.

![Alt Text](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/energy_landscape_3D.gif)


- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 1.
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 2.
- **Z-Axis**: Shows the free energy in kJ/mol, calculated from the probability distribution of CV1 and CV2 using the Boltzmann equation.

**Visualization**: A 3D surface plot illustrates the variations in free energy across the space defined by CV1 and CV2. The color gradient enhances the visualization of energy barriers and minima, aiding in the identification of stable conformations and transition states.

**Animated 3D Free Energy Landscape**
The animated GIF provides a rotating view of the 3D free energy landscape, offering an immersive exploration of the energy barriers and minima across the collective variable space. The animation helps in visualizing the landscape's depth and complexity from multiple angles, enhancing the understanding of the molecular system's energetics.

The 3D and animated visualizations of the free energy landscape extend the analysis beyond two dimensions, offering a richer and more nuanced understanding of the system's energetics. By exploring the landscape in three dimensions, researchers can better identify and analyze the regions of interest, such as low-energy conformations and high-energy transition states, which are vital for deciphering the molecular mechanisms underlying biological processes and material behaviors.

## Interpretation

Together, these figures provide a multi-faceted view of the molecular system's free energy landscape. Analyzing these visualizations helps in understanding how variations in critical structural parameters (angles and distances) influence the stability and dynamics of the system, which is vital for predicting reaction pathways, designing drugs, and engineering materials with desired properties.
