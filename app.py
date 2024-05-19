import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

# Set page configuration
st.set_page_config(page_title="Neighbor Joining", page_icon='🧬', layout='wide')

# Function to generate text content for the download file
def generate_text_content(matrix, J_matrix, min_value, min_indices, new_distances, updated_labels):
    content = "Tree Distances and Groupings:\n\n"
    
    # Add information about the J matrix
    content += "J Matrix:\n"
    content += str(pd.DataFrame(J_matrix, columns=[f'Seq {i+1}' for i in range(len(J_matrix))], index=[f'Seq {i+1}' for i in range(len(J_matrix))])) + "\n\n"
    
    # Add information about the minimum value in the J matrix
    content += f"The minimum value in the J matrix is {min_value} at position ({min_indices[0][0]+1}, {min_indices[1][0]+1}).\n\n"
    
    # Add information about new distances
    content += "New Distances:\n"
    for idx, distance in enumerate(new_distances):
        content += f"Distance from {updated_labels[idx]} to new node {chr(97 + idx)}: {distance}\n"
    
    return content

# Recursive function for neighbor joining
def neighbor_joining_recursive(matrix, node_label='a'):
    n = matrix.shape[0]
    
    if n == 2:
        st.warning("The matrix size is 2x2. The recursion stops here.")
        return
    
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

    # Update the distance matrix
    indices_to_keep = [idx for idx in range(n) if idx != i and idx != j]
    new_matrix = matrix[indices_to_keep, :][:, indices_to_keep]
    new_matrix = np.vstack([new_matrix, new_distances])
    new_col = np.append(new_distances, 0)
    new_matrix = np.column_stack([new_matrix, new_col])

    # Update sequence labels
    updated_labels = [f'Seq {k+1}' for k in indices_to_keep] + [f'node {node_label}']
    
    # Display the updated distance matrix
    st.write("The updated distance matrix is:")
    new_df = pd.DataFrame(new_matrix, columns=updated_labels, index=updated_labels)
    st.dataframe(new_df)
    
    # Generate text content for download
    text_content = generate_text_content(matrix, J_matrix, min_value, min_indices, new_distances, updated_labels)
    
    # Add a download button
    st.sidebar.download_button("Download Tree Distances", data=text_content, file_name='tree_distances.txt', mime='text/plain')
    
    # Recursively call the function with the updated matrix and next node label
    neighbor_joining_recursive(new_matrix, chr(ord(node_label) + 1))

st.title("Neighbor Joining Distance Matrix")

# Sidebar for inputs
with st.sidebar:
    st.header("Inputs")
    
    # Step 1: User input for the distance matrix
    st.subheader("Step 1: Enter the distance matrix")
    matrix_input = st.text_area("Paste the distance matrix here (rows separated by newlines and values by spaces or commas):")

if matrix_input:
    # Convert the input string into a numpy array
    try:
        matrix = np.array([list(map(float, row.split())) for row in matrix_input.strip().split('\n')])
        
        # Validate that the matrix is square
        if matrix.shape[0] != matrix.shape[1]:
            st.error("The provided matrix is not square. Please ensure you enter a valid distance matrix.")
        else:
            # Display the initial distance matrix
            st.write("The initial distance matrix is:")
            df = pd.DataFrame(matrix, columns=[f'Seq {i+1}' for i in range(matrix.shape[0])], index=[f'Seq {i+1}' for i in range(matrix.shape[0])])
            st.dataframe(df)
            
            # Start the neighbor joining process recursively
            neighbor_joining_recursive(matrix)

    except ValueError:
        st.error("There was an error processing the matrix. Please ensure all values are numerical and the matrix format is correct.")
