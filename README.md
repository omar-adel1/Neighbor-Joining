# Neighbor-Joining Phylogenetic Tree Construction

This project implements the Neighbor-Joining (NJ) algorithm for constructing phylogenetic trees based on distance matrices. It includes two main components: a script (`main.py`) for generating NJ and UPGMA trees using BioPython and an interactive web application (`app.py`) built with Streamlit for visualizing and analyzing phylogenetic trees.

## Installation

### Requirements

- Python 3.7+
- Required Python libraries: BioPython, NumPy, pandas, Matplotlib, Streamlit

### Setup

1. Clone the repository
2. Set up a virtual environment (optional but recommended):
    ```bash
     python -m venv venv
    source venv/bin/activate  # On Windows, use venv\Scripts\activate  
     ```
3. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```
## Usage

### Using `main.py`

`main.py` implements the Neighbor-Joining and UPGMA algorithms using predefined distance matrices. It generates phylogenetic trees and outputs tree details to `tree_details.txt`.

Run the script:

```bash
python main.py
```
### Using `app.py`

`app.py` provides an interactive web interface for constructing and visualizing phylogenetic trees based on user-provided distance matrices. It uses Streamlit to display the initial matrix and the resulting phylogenetic tree.

Run the application:

```bash
streamlit run app.py
```

## Input Requirements

To successfully run the Neighbor-Joining algorithm using the provided scripts (`main.py` and `app.py`), ensure the following input requirements are met:

- **Distance Matrix**: The distance matrix should be provided in a square format.
- **Matrix Format**: Each row of the matrix should represent distances from a single sequence to all other sequences, separated by spaces or commas.
- **Diagonal Elements**: Ensure all diagonal elements of the matrix are set to zero, as they represent the distance from a sequence to itself.
- **Symmetry**: The values below the diagonal should mirror those above it, ensuring symmetric distances between sequences.
