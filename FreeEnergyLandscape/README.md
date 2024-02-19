# Free Energy Landscape Analysis

This Python class `FreeEnergyLandscape`, is designed to analyze and visualize the free energy landscape of molecular systems using collective variables (CVs), such as angles and distances, obtained from molecular dynamics simulations. The approach leverages the concept of Boltzmann inversion to transform probability distributions into free energy surfaces, providing insights into the thermodynamic and kinetic properties of the system.

## Theoretical Background

The free energy landscape is a conceptual and computational tool used to understand the energetics and dynamics of molecular systems. It is particularly useful in the study of complex processes such as protein folding, chemical reactions, and phase transitions.

### Collective Variables (CVs)

Collective Variables (CVs) are a set of coordinates that describe the macroscopic state of a system. They are used to reduce the complexity of molecular systems by focusing on the relevant degrees of freedom. Examples include the distance between two atoms, angles, dihedrals, and more complex descriptors.

CVs in biomolecular systems are mathematically represented as functions of the atomic coordinates. These functions are designed to capture the essential features of the system's configuration that are relevant to its macroscopic properties or behaviors. Below are examples of commonly used CVs and their mathematical formulations:

1. **Distance**: The distance $d$ between two atoms $i$ and $j$ with positions $\vec{r}_i$ and $\vec{r}_j$ can be calculated using the Euclidean distance formula:

   First, determine the position vectors of atoms $i$ and $j$:
   $$\vec{r}_i = (x_i, y_i, z_i)$$
   $$\vec{r}_j = (x_j, y_j, z_j)$$

   Then, calculate the distance $d$ as:
   $$d = |\vec{r}_i - \vec{r}_j| = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}$$

   - $\vec{r}_i$ and $\vec{r}_j$ are the position vectors of atoms $i$ and $j$, respectively.
   - $(x_i, y_i, z_i)$ and $(x_j, y_j, z_j)$ denote the Cartesian coordinates of atoms $i$ and $j$.
   - $|\vec{r}_i - \vec{r}_j|$ represents the magnitude of the vector difference between $\vec{r}_i$ and $\vec{r}_j$, giving the direct distance between the two atoms.

2. **Angle**: The angle $\theta$ formed by three atoms $i$, $j$, and $k$, where $j$ is the vertex, can be calculated using the dot product:

   First, determine the vectors $\vec{r}_ji$ and $\vec{r}_jk$:
   $$\vec{r}_ji = \vec{r}_i - \vec{r}_j$$
   $$\vec{r}_jk = \vec{r}_k - \vec{r}_j$$

   Then, calculate the angle $\theta$ as:
   $$\cos(\theta) = \frac{\vec{r}_ji \cdot \vec{r}_jk}{|\vec{r}_ji| |\vec{r}_jk|}$$
   $$\theta = \arccos\left(\frac{\vec{r}_ji \cdot \vec{r}_jk}{|\vec{r}_ji| |\vec{r}_jk|}\right)$$

   - $\vec{r}_ji$ and $\vec{r}_jk$ are vectors pointing from atom $j$ to atoms $i$ and $k$, respectively.
   - $\cdot$ denotes the dot product between the vectors $\vec{r}_ji$ and $\vec{r}_jk$.
   - $|\vec{r}_ji|$ and $|\vec{r}_jk|$ represent the magnitudes of the vectors $\vec{r}_ji$ and $\vec{r}_jk$, respectively.
   - $\arccos$ is the inverse cosine function, used to find the angle $\theta$ from the cosine value.

These representations allow for a simplified description of the system's state, facilitating the study of its behavior and properties through the manipulation of a reduced set of variables rather than the full set of atomic coordinates.

### Kernel Density Estimation (KDE)

This statistical method is crucial for estimating the probability density function (PDF) of a dataset without assuming any predefined distribution shape. In the context of molecular dynamics and simulations, it allows us to visualize and analyze the distribution of free energy across different states or configurations defined by collective variables (CVs).

The `gaussian_kde` leverages a Gaussian (normal) distribution, placing it at each point in the dataset and summing these distributions to approximate the overall data's PDF. This technique is adept at capturing the underlying structure of the data, providing a smooth, continuous representation of the free energy landscape. The smoothness of the KDE is controlled by the bandwidth parameter, which determines the width of the Gaussian kernels used.

Mathematically, the density estimation at a point $`x`$ is calculated as follows:

$$\ f(x) = \frac{1}{n \cdot h} \sum_{i=1}^{n} K(u) \$$   

where:
- $K(u)$ represents the Gaussian kernel function

$$ K(u) = \frac{1}{\sqrt{2\pi}} e^{-\frac{1}{2}u^2} $$

- $u$ is the normalized distance

$$ u = \frac{x - x_i}{h} $$

- $n$ is the number of data points
- $h$ is the bandwidth
- $x$ is where the density is being estimated
- $x_i$ are the data points
  
note that:
- $u^2$ is the square of this normalized distance, which serves to weight the contribution of each data point $x_i$ to the density estimate at $x$, based on how far $x_i$ is from $x$, adjusted by the bandwidth $h$.
- $e^{-\frac{1}{2}u^2}$ decreases rapidly as $u$ increases, meaning that points further away from $x$ will have less influence on the density estimate at $x$.
- $\frac{1}{\sqrt{2\pi}}$ is a normalization term that ensures the Gaussian kernel function integrates to $1$, keeping it as a valid probability distribution.

In the context of KDE the variable $u$ is utilized within the Gaussian kernel function $K$, which is a probability density function. Thus, $u$ is crucial in determining how the distance between $x$ and the data points $x_i$ affects the density estimate at $x$, with the bandwidth $h$ controlling the sensitivity of this influence.

The Gaussian KDE method provides a sophisticated approach to model the complex free energy landscapes encountered in molecular dynamics studies. It enables researchers to visualize the distribution of energy states without the constraints of parametric models, offering insights into molecular stability, transitions, and the energetics of molecular interactions.

### Boltzmann distribution

The Boltzmann distribution relates the energy of a state to its probability in a canonical ensemble, this provides a fundamental link between the microstate probabilities of a thermodynamic system and its macroscopic properties. It expresses the probability $P$ of a system being in a state with a certain free energy $\Delta G$ at a specific temperature $T$:

$$P \propto e^{-\frac{\Delta G}{k_B T}}$$

where:
- $P$ is the probability of observing the state,
- $\Delta G$ represents the free energy difference of that state relative to a chosen reference state,
- $k_B$ is the Boltzmann constant, a fundamental physical constant that relates energy scales to temperature,
- $T$ is the absolute temperature of the system, measured in Kelvin.

### Normalization of Probability

For the concept of probability to be meaningful in the context of statistical mechanics, it's imperative that the probabilities of all conceivable states sum to unity. This requirement ensures that the predicted behavior encompasses all possible configurations of the system. Consider a state that represents the maximum free energy, denoted as $\Delta G_{\text{max}}$, its corresponding maximum probability, $P_{\text{max}}$, can be analogously expressed through the Boltzmann factor:

$$P_{\text{max}} \propto e^{-\frac{\Delta G_{\text{max}}}{k_B T}}$$

### From Proportionality to Quantitative Relationship

To quantitatively relate $\Delta G$ with $P$ and $P_{\text{max}}$, we proceed as follows:

1. Express $P$ as a function of $\Delta G$:

   $$P = e^{-\frac{\Delta G}{k_B T}}$$

2. Similarly, express $P_{\text{max}}$ as a function of $\Delta G_{\text{max}}$:

   $$P_{\text{max}} = e^{-\frac{\Delta G_{\text{max}}}{k_B T}}$$

Taking the natural logarithm on both sides of these expressions yields:

$$\ln(P) = -\frac{\Delta G}{k_B T}$$

$$\ln(P_{\text{max}}) = -\frac{\Delta G_{\text{max}}}{k_B T}$$

Subtracting the equation for $\ln(P_{\text{max}})$ from the equation for $\ln(P)$ removes the dependency on the reference state $\Delta G_{\text{max}}$, leading to:

$$\ln(P) - \ln(P_{\text{max}}) = -\frac{\Delta G}{k_B T} + \frac{\Delta G_{\text{max}}}{k_B T}$$

Given the interest in determining the free energy difference $\Delta G$ relative to the state with maximum probability $P_{\text{max}}$, we can rearrange this relationship to solve explicitly for $\Delta G$:

$$\Delta G = -k_B T (\ln(P) - \ln(P_{\text{max}}))$$

The Gibbs free energy provides a measure of the amount of "useful work" that can be obtained from a thermodynamic system in a process at constant temperature and pressure. Here, $\Delta G$ is expressed as a function of the difference in the logarithms of the probabilities of finding the system in any state versus the state of highest probability (or lowest free energy), multiplied by the system's absolute temperature and the Boltzmann constant. This relationship shows how the free energy difference between two states can be calculated from their relative probabilities, providing a direct bridge between statistical thermodynamics and experimental observations or computational simulations.

This relationship allows us to convert a histogram of CV values into a free energy surface and this formulation also adjusts the free energy values so that the minimum energy associated with the highest probability $P_{\text{max}}$ is set to 0. Because the value of $\Delta G$ for the state with $P_{\text{max}}$ will be 0 because $\ln (P_{\text{max}}) - \ln (P_{\text{max}}) = 0\$.

This expression quantitatively links the probability distribution of states within a thermodynamic system to their respective free energy differences, providing a foundation for analyzing the system's behavior at the molecular level and this approach is useful for highlighting the relative differences in free energy between different states or conformations in an MD simulation, facilitating the identification of free energy minima and the relative comparison between different states. By setting the minimum of free energy to 0, you create a clear reference point to evaluate the relative stability of other states compared to the most stable state.


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
These dependencies also are listed in the `requirements.txt` file. To install them, run the following command in your terminal:

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
