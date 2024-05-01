import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import joblib
import time
import gc
import os
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA  # Import PCA
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, make_scorer, balanced_accuracy_score
from sklearn.metrics import roc_auc_score
from imblearn.metrics import geometric_mean_score  # Import from imbalanced-learn
from sklearn.pipeline import Pipeline
from sklearn.utils import parallel_backend

# Importing various classifiers
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier

class Classifiers:
    def __init__(self, data_path, batch_size=128):
        self.data_path = data_path
        self.batch_size = batch_size
        self.load_data()
        self.kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        self.save_directory = './model_outputs'  # Directory to save models and results
        os.makedirs(self.save_directory, exist_ok=True)  # Ensure the directory exists

    def load_data(self):
        loader = np.load(self.data_path, allow_pickle=True)
        self.X_train, self.X_test = loader['X_train'], loader['X_test']
        self.y_train, self.y_test = loader['y_train'], loader['y_test']

    def evaluate_model(self, model, name):
"""     
        Definição de Métricas: Define um dicionário de métricas de avaliação para serem calculadas durante a validação cruzada. Cada métrica ajuda a entender diferentes aspectos do desempenho do modelo:
        Accuracy: Proporção de predições corretas.
        Precision: Proporção de predições positivas corretas em relação a todas as predições positivas.
        Recall: Proporção de predições positivas corretas em relação a todas as amostras positivas reais.
        F1: Média harmônica de precisão e recall.
        Balanced Accuracy: Média de recall obtida em cada classe, útil para dados desbalanceados.
        Geometric Mean: Média geométrica de recall para todas as classes.
"""
        metrics = {
            'accuracy': 'accuracy',
            'precision': make_scorer(precision_score, average='macro'),
            'recall': make_scorer(recall_score, average='macro'),
            'f1': make_scorer(f1_score, average='macro'),
            'balanced_accuracy': 'balanced_accuracy',
            'geometric_mean': make_scorer(geometric_mean_score, average='macro')
        }
        scores = cross_validate(model, self.X_train, self.y_train, cv=self.kfold,
                                scoring=metrics, return_train_score=False,
                                return_estimator=True, n_jobs=1)
        return {metric: np.mean(scores[f'test_{metric}']) for metric in metrics}, scores['estimator'][-1]

    def train_and_evaluate(self):
        start_time = time.time()
        results = []
        model_details = []
        num_batches = max(1, len(self.X_train) // self.batch_size)

        models = [
            ("Nearest Neighbors", KNeighborsClassifier(3)),
            ("Linear SVM", SVC(kernel="linear", C=0.025, probability=True)),
            ("RBF SVM", SVC(gamma=2, C=1, probability=True)),
            ("Gaussian Process", GaussianProcessClassifier(1.0 * RBF(1.0))),
            ("Decision Tree", DecisionTreeClassifier(max_depth=5)),
            ("Random Forest", RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)),
            ("Neural Net", MLPClassifier(alpha=1, max_iter=1000)),
            ("AdaBoost", AdaBoostClassifier()),
            ("Naive Bayes", GaussianNB(var_smoothing=1e-8)),
            ("QDA", QuadraticDiscriminantAnalysis()),
            ("Logistic Regression", LogisticRegression(max_iter=1000)),
            ("Gradient Boosting", GradientBoostingClassifier())
        ]

        for name, clf in tqdm(models, desc="Processing classifiers"):
            model = Pipeline([
                ('scaler', StandardScaler()),
                ('pca', PCA(n_components=0.95)),  # Reducing dimensionality to retain 95% of variance
                ('classifier', clf)
            ])
            batch_results = []
            for batch in np.array_split(np.arange(len(self.X_train)), num_batches):
                X_batch, y_batch = self.X_train[batch], self.y_train[batch]
                batch_result, fitted_model = self.evaluate_model(model, name)
                batch_results.append(batch_result)
                del X_batch, y_batch
                gc.collect()

            avg_results = {metric: np.mean([r[metric] for r in batch_results]) for metric in batch_results[0]}
            results.append((name, avg_results))
            model_details.append((name, avg_results, fitted_model))

            joblib.dump(fitted_model, os.path.join(self.save_directory, f'model_{name}.joblib'))
            np.save(os.path.join(self.save_directory, f'results_{name}.npy'), avg_results)

        end_time = time.time()
        print(f"Total execution time: {end_time - start_time:.2f} seconds")
        return results

def plot_metrics(directory):
    metrics = ['accuracy']
    fig, ax = plt.subplots(figsize=(10, 6))
    names = []
    accuracies = []

    for result_file in os.listdir(directory):
        if result_file.startswith('results_') and result_file.endswith('.npy'):
            name = result_file.replace('results_', '').replace('.npy', '')
            results = np.load(os.path.join(directory, result_file), allow_pickle=True).item()
            accuracy = results['accuracy']
            names.append(name)
            accuracies.append(accuracy)

    ax.barh(names, accuracies, color='skyblue')
    ax.set_title('Accuracy Score')
    ax.set_xlabel('Score')
    ax.set_xlim(0, 1.0)
    plt.tight_layout()
    plt.show()

def main():
    classifier = Classifiers('./split_data.npz', batch_size=128)
    results = classifier.train_and_evaluate()
    plot_metrics(classifier.save_directory)

if __name__ == '__main__':
    main()
