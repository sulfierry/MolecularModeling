# Kinase Target Prediction System

This system is designed for predicting kinase targets and groups for given ligands using their chemical structures. The prediction is based on machine learning models that are trained on molecular fingerprints.

## Modules Description

### 1. Data Balance Checking (`1_check_data_balance.py`)

- **`CheckBalance` class**: Manages the loading and processing of chemical data, and checks for balance in dataset classification labels.
  
  - **`load_data(filepath)`**: Loads specific columns from a TSV file necessary for further processing and analysis.
  
  - **`analyze_class_balance(data, activity_threshold)`**: Analyzes and visualizes the balance of classes within the dataset both overall and by kinase group.
  
  - **`analyze_class_metrics(data, activity_threshold)`**: Calculates statistical metrics such as class ratio, entropy of label distribution, and coefficient of variation to assess the balance in the dataset.
  
  - **`save_output(data)`**: Saves the filtered and processed dataset to a local file after removing specified kinase group entries.

### 2. Data Preprocessing (`2_ml_input.py`)

- **`load_data(filepath)`**: Loads data with predefined columns including chemical IDs, SMILES notation, kinase information, and labels from a file.
  
- **`smiles_to_fingerprints(smiles, radius, n_bits)`**: Converts SMILES strings to molecular fingerprints, handling missing values and generating bit vector representations.
  
- **`convert_smiles_batch(smiles_list, radius, n_bits)`**: Utilizes multiple processors to convert batches of SMILES strings into fingerprints, improving computational efficiency.
  
- **`preprocess_data(data)`**: Manages batch processing of SMILES to fingerprints, prepares labels from categorical to binary format, and integrates all preprocessing steps.
  
- **`split_data(data)`**: Splits the dataset into training and testing sets, based on processed fingerprints and binary labels, facilitating model training and evaluation.

### 3. Model Training and Evaluation (`3_modelling.py`)

- **`Classifiers` class**: Handles loading of data and training of machine learning models for kinase prediction.
  
  - **`load_data()`**: Loads training and testing sets from a `.npz` file containing molecular fingerprints and binary labels.
  
  - **`evaluate_model(model, name)`**: Conducts cross-validation on the training data to compute several performance metrics such as accuracy, precision, recall, F1 score, balanced accuracy, and geometric mean.
  
  - **`train_and_evaluate()`**: Trains multiple predefined models, evaluates them, and saves the top-performing models based on accuracy. The method visualizes performance metrics across all models.

## Usage

Each script is designed to be run sequentially to ensure data is processed, balanced, and used for training machine learning models efficiently. Ensure all prerequisites such as required libraries (pandas, NumPy, scikit-learn, RDKit, joblib) are installed and data paths are correctly set in the scripts.
