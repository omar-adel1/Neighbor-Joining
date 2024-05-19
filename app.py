import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO


# Set page configuration
st.set_page_config(page_title="Neighbor Joining", page_icon='🧬', layout='wide')

# Function to generate text content for the download file
def generate_text_content(matrix, tree, labels):
    content = "Phylogenetic tree\n\n"
    
    # Add information about the distances between sequences
    computed_pairs = set()  # Keep track of computed pairs
    content += "Distances between sequences:\n"
    for idx, label in enumerate(labels):
        for idx2 in range(idx + 1, len(labels)):
            pair = tuple(sorted([label, labels[idx2]]))
            if pair not in computed_pairs:
                distance = tree.distance(label, labels[idx2])
                content += f"Distance from {label} to {labels[idx2]}: {distance}\n"
                computed_pairs.add(pair)
    
    # Add information about the distances between nodes and other sequences
    computed_node_pairs = set()  # Keep track of computed node-sequence pairs
    for clade in tree.find_clades(order="level"):
        if clade.name.startswith("node"):
            node_label = clade.name.split()[1]
            for label in labels:
                pair = tuple(sorted([f"node {node_label}", label]))
                if pair not in computed_node_pairs:
                    distance = tree.distance(f"node {node_label}", label)
                    content += f"Distance from node {node_label} to {label}: {distance}\n"
                    computed_node_pairs.add(pair)
    
    return content


# Function to draw the phylogenetic tree using BioPython
def draw_tree(tree):
    fig, ax = plt.subplots(figsize=(10, 8))
    Phylo.draw(tree, axes=ax)
    ax.axis('off')  # Turn off axis
    st.pyplot(fig)

# Recursive function for neighbor joining
def neighbor_joining_recursive(matrix, labels, node_label='a', tree=None, parent_node=None):
    n = matrix.shape[0]
    
    if n == 2:
        i, j = 0, 1
        distance = matrix[i, j] / 2
        clade_i = Phylo.Newick.Clade(name=labels[i], branch_length=distance)
        clade_j = Phylo.Newick.Clade(name=labels[j], branch_length=distance)
        clade = Phylo.Newick.Clade(name=f'node {node_label}', clades=[clade_i, clade_j])
        if parent_node:
            parent_node.clades.append(clade)
        else:
            tree = Phylo.Newick.Tree(root=clade)
        return tree
    
    # Calculate the J matrix
    row_sums = np.sum(matrix, axis=1)
    J_matrix = np.zeros_like(matrix)
    for i in range(n):
        for j in range(n):
            if i != j:
                J_matrix[i, j] = (n - 2) * matrix[i, j] - row_sums[i] - row_sums[j]
            else:
                J_matrix[i, j] = 0

    # Find the minimum value in the J matrix
    min_value = np.min(J_matrix[np.nonzero(J_matrix)])
    min_indices = np.where(J_matrix == min_value)
    i, j = min_indices[0][0], min_indices[1][0]

    # Calculate the distance from each sequence to the new node
    new_distances = []
    for k in range(n):
        if k != i and k != j:
            new_distance = 0.5 * (matrix[i, k] + matrix[j, k] - matrix[i, j])
            new_distances.append(new_distance)

    # Create a new node in the tree
    distance_i = 0.5 * (matrix[i, j] + (row_sums[i] - row_sums[j]) / (n - 2))
    distance_j = matrix[i, j] - distance_i
    clade_i = Phylo.Newick.Clade(name=labels[i], branch_length=distance_i)
    clade_j = Phylo.Newick.Clade(name=labels[j], branch_length=distance_j)
    new_clade = Phylo.Newick.Clade(name=f'node {node_label}', clades=[clade_i, clade_j])

    if tree is None:
        tree = Phylo.Newick.Tree(root=new_clade)
    else:
        parent_node.clades.append(new_clade)

    # Update the distance matrix
    indices_to_keep = [idx for idx in range(n) if idx != i and idx != j]
    new_matrix = matrix[indices_to_keep, :][:, indices_to_keep]
    new_matrix = np.vstack([new_matrix, new_distances])
    new_col = np.append(new_distances, 0)
    new_matrix = np.column_stack([new_matrix, new_col])

    # Update sequence labels
    updated_labels = [labels[k] for k in indices_to_keep] + [f'node {node_label}']

    # Recursively call the function with the updated matrix and next node label
    return neighbor_joining_recursive(new_matrix, updated_labels, chr(ord(node_label) + 1), tree, new_clade)

st.title("Neighbor Joining Distance Matrix")

# Sidebar for inputs
with st.sidebar:
    st.header("Inputs")
    
    # Step 1: User
    # input for the distance matrix
    st.subheader("Step 1: Enter the distance matrix")
    matrix_input = st.text_area("Paste the distance matrix here (rows separated by newlines and values by spaces or commas):")

if matrix_input:
    # Convert the input string into a numpy array
    try:
        matrix = np.array([list(map(float, row.split())) for row in matrix_input.strip().split('\n')])
        
        # Validate that the matrix is square
        if matrix.shape[0] != matrix.shape[1]:
            st.error("The provided matrix is not square. Please ensure you enter a valid distance matrix.")
         # Validate that the diagonal elements are all zeros
        elif not np.all(np.diag(matrix) == 0):
            st.error("The diagonal elements of the matrix should all be zeros. Please ensure you enter a valid distance matrix.")
        # Validate that the values below the diagonal are symmetric with the values above the diagonal
        elif not np.allclose(matrix, matrix.T, atol=1e-8):
            st.error("The values below the diagonal are not symmetric with the values above the diagonal. Please ensure you enter a valid distance matrix.")
        else:
            col1, col2 = st.columns(2)
            with col1:
                # Display the initial distance matrix
                st.subheader("The initial distance matrix is:")
                df = pd.DataFrame(matrix, columns=[f'Seq {i+1}' for i in range(matrix.shape[0])], index=[f'Seq {i+1}' for i in range(matrix.shape[0])])
                st.dataframe(df)
                
                # Initialize labels
                labels = [f'Seq {i+1}' for i in range(matrix.shape[0])]
                
                # Start the neighbor joining process recursively
                tree = neighbor_joining_recursive(matrix, labels)

            with col2:
                # Draw the final tree
                if tree:
                    st.subheader("Phylogenetic Tree")
                    draw_tree(tree)
                
            # Generate text content for download
            text_content = generate_text_content(matrix, tree, labels)

            # Add a download button
            st.sidebar.download_button("Download Tree Distances", data=text_content, file_name='tree_distances.txt', mime='text/plain')

    except ValueError:
        st.error("There was an error processing the matrix. Please ensure all values are numerical and the matrix format is correct.")
