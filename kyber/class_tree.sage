# Node class for the residue tree
class ResidueTreeNode:
    def __init__(self, residue):
        self.residue = residue  # Value stored in the node (a residue)
        self.left = None        # Left child node
        self.right = None       # Right child node

# Class representing the entire residue tree structure
class ResidueTree:
    def __init__(self, root_residue):
        self.root = ResidueTreeNode(root_residue)  # Initialize tree with root node

    # Inserts left and/or right child nodes under a given parent residue
    def insert(self, parent_residue, left_residue=None, right_residue=None):
        parent_node = self.find(self.root, parent_residue)
        if parent_node:
            if left_residue is not None:
                parent_node.left = ResidueTreeNode(left_residue)
            if right_residue is not None:
                parent_node.right = ResidueTreeNode(right_residue)

    # Recursive search function to locate a node with a given residue value
    def find(self, node, residue):
        if node is None:
            return None
        if node.residue == residue:
            return node
        return self.find(node.left, residue) or self.find(node.right, residue)

    # Displays the tree structure in a readable format
    def display(self, node=None, level=0):
        if node is None:
            node = self.root
        print("  " * level + f"({node.residue})")
        if node.left:
            self.display(node.left, level + 1)
        if node.right:
            self.display(node.right, level + 1)

# Recursively builds a binary residue tree using modular square roots and inverses
def build_residue_tree(tree, node, mod_value, depth, r_list):
    if depth == 2:
        r_list.append(node.residue)  
        return

    try:
        r1 = mod(-node.residue, mod_value).sqrt()  
        r_inv = R(r1)^-1  
        ri = (node.residue * r_inv) % mod_value  

        # Insert left and right children into the tree
        tree.insert(node.residue, left_residue=ri, right_residue=r1)

        # Recursively build left and right subtrees
        left_node = tree.find(tree.root, r1)
        right_node = tree.find(tree.root, ri)
        build_residue_tree(tree, left_node, mod_value, depth // 2, r_list)
        build_residue_tree(tree, right_node, mod_value, depth // 2, r_list)

    except ValueError:
        # Raised when sqrt does not exist in the field
        return

# Creates a residue tree with specified root (alpha), degree, and modulus q1
def build_tree(alpha, degree, q1):
    r_list = []
    tree = ResidueTree(alpha)
    build_residue_tree(tree, tree.root, q1, degree, r_list)
    return tree, r_list  

# Recursively inverts the residue values in the tree (modular inverse modulo q)
def invert_residue_tree(original_root, q):
    if original_root is None:
        return None

    try:
        inv_residue = original_root.residue^-1 % q 
    except ZeroDivisionError:
        return None  # Skip if inverse doesn't exist (zero residue)

    new_node = ResidueTreeNode(inv_residue)
    new_node.left = invert_residue_tree(original_root.left, q)
    new_node.right = invert_residue_tree(original_root.right, q)
    return new_node

# Builds a new tree where each residue is the modular inverse of the original tree's residues
def build_inverse_tree(original_tree, q, alpha):
    inverse_root = invert_residue_tree(original_tree.root, q)
    inverse_tree = ResidueTree(alpha)  
    inverse_tree.root = inverse_root
    return inverse_tree 

