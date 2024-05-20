from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio import Phylo

# Define the distance matrix
dm = DistanceMatrix(names=['Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon'],
                    matrix=[[0],
                            [0.182692307692, 0],
                            [0.115384615385, 0.115384615385, 0],
                            [0.0769230769231, 0.115384615385, 0.0769230769231, 0],
                            [0.0512820512821, 0.102564102564, 0.0769230769231, 0.0769230769231, 0]])

# Initialize the constructor
constructor = DistanceTreeConstructor()

# Construct the trees
nj_tree = constructor.nj(dm)
upgma_tree = constructor.upgma(dm)

# Draw the NJ tree in a window
print("Neighbor Joining Tree:")
Phylo.draw(nj_tree)

# Draw the UPGMA tree in a window
print("\nUPGMA Tree:")
Phylo.draw(upgma_tree)

# Open a text file to write the details
with open("tree_details.txt", "w") as file:
    # Write Neighbor Joining Tree details
    file.write("Neighbor Joining Tree:\n")
    file.write("Tree distances between sequences:\n")
    file.write(str(dm) + "\n\n")
    file.write("Grouping of nodes to parents:\n")
    for clade in nj_tree.find_clades():
        if clade.name:
            file.write(f"{clade.name} -> {clade.clades}\n")
    file.write("\n\n")
    
    # Write UPGMA Tree details
    file.write("UPGMA Tree:\n")
    file.write("Tree distances between sequences:\n")
    file.write(str(dm) + "\n\n")
    file.write("Grouping of nodes to parents:\n")
    for clade in upgma_tree.find_clades():
        if clade.name:
            file.write(f"{clade.name} -> {clade.clades}\n")
