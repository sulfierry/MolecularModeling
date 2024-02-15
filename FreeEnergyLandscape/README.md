# Free Energy Landscape Analysis

This Python class, `FreeEnergyLandscape`, is designed to analyze and visualize the free energy landscape of molecular systems using collective variables (CVs), such as angles and distances, obtained from molecular dynamics simulations. The approach leverages the concept of Boltzmann inversion to transform probability distributions into free energy surfaces, providing insights into the thermodynamic and kinetic properties of the system.

## Theoretical Background

The free energy landscape is a conceptual and computational tool used to understand the energetics and dynamics of molecular systems. It is particularly useful in the study of complex processes such as protein folding, chemical reactions, and phase transitions.

### Collective Variables (CVs)

Collective Variables (CVs) are a set of coordinates that describe the macroscopic state of a system. They are used to reduce the complexity of molecular systems by focusing on the relevant degrees of freedom. Examples include the distance between two atoms, angles, dihedrals, and more complex descriptors.

### Boltzmann Inversion

The Boltzmann inversion method is used to calculate the free energy landscape from the probability distribution of CVs. The free energy, \( G \), of a state is related to its probability, \( P \), by the Boltzmann equation:

\[ G = -k_B T \ln(P) \]

where \( k_B \) is the Boltzmann constant (\( 1.380649 \times 10^{-23} \, \text{J/K} \)), \( T \) is the temperature, and \( P \) is the probability density of the CVs. This relationship allows us to convert a histogram of CV values into a free energy surface.

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
    t = 300         # Temperature in K
    kB = 8.314e-3   # Boltzmann constant in kJ/(mol·K)

    cv1_path = './proj1Out.txt'
    cv2_path = './proj2Out.txt'

    fel = FreeEnergyLandscape(cv1_path, cv2_path, t, kB)
    fel.main()
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

![Alt text da image](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/fel3D.png)


- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 1.
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 2.
- **Z-Axis**: Shows the free energy in kJ/mol, calculated from the probability distribution of CV1 and CV2 using the Boltzmann equation.

**Visualization**: A 3D surface plot illustrates the variations in free energy across the space defined by CV1 and CV2. The color gradient enhances the visualization of energy barriers and minima, aiding in the identification of stable conformations and transition states.

**Animated 3D Free Energy Landscape**
The animated GIF provides a rotating view of the 3D free energy landscape, offering an immersive exploration of the energy barriers and minima across the collective variable space. The animation helps in visualizing the landscape's depth and complexity from multiple angles, enhancing the understanding of the molecular system's energetics.

**Animation**: The GIF loops through a series of rotations around the Z-axis, presenting a 360-degree view of the free energy landscape. This continuous motion provides a comprehensive perspective on the distribution of energy states.

![Alt Text](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/energy_landscape_3D.gif)

**Purpose**: The animated visualization aids in grasping the multidimensional nature of the free energy landscape, which is crucial for understanding the thermodynamics and kinetics of molecular systems.
Interpretation

The 3D and animated visualizations of the free energy landscape extend the analysis beyond two dimensions, offering a richer and more nuanced understanding of the system's energetics. By exploring the landscape in three dimensions, researchers can better identify and analyze the regions of interest, such as low-energy conformations and high-energy transition states, which are vital for deciphering the molecular mechanisms underlying biological processes and material behaviors.

## Interpretation

Together, these figures provide a multi-faceted view of the molecular system's free energy landscape. Analyzing these visualizations helps in understanding how variations in critical structural parameters (angles and distances) influence the stability and dynamics of the system, which is vital for predicting reaction pathways, designing drugs, and engineering materials with desired properties.
