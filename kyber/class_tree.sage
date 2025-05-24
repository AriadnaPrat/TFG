class ResidueTreeNode:
    def __init__(self, residue):
        self.residue = residue  
        self.left = None  
        self.right = None  

class ResidueTree:
    def __init__(self, root_residue):
        self.root = ResidueTreeNode(root_residue)

    def insert(self, parent_residue, left_residue=None, right_residue=None):
        parent_node = self.find(self.root, parent_residue)
        if parent_node:
            if left_residue is not None:
                parent_node.left = ResidueTreeNode(left_residue)
            if right_residue is not None:
                parent_node.right = ResidueTreeNode(right_residue)

    def find(self, node, residue):
        if node is None:
            return None
        if node.residue == residue:
            return node
        return self.find(node.left, residue) or self.find(node.right, residue)

    def display(self, node=None, level=0):
        if node is None:
            node = self.root
        print("  " * level + f"({node.residue})")
        if node.left:
            self.display(node.left, level + 1)
        if node.right:
            self.display(node.right, level + 1)

def build_residue_tree(tree, node, mod_value, depth, r_list):
    if depth == 2:
        r_list.append(node.residue)
        return

    try:
        r1 = mod(-node.residue, mod_value).sqrt()  
        r_inv = R(r1)^-1  # r1^-1
        ri = (node.residue * r_inv) % mod_value  
        tree.insert(node.residue, left_residue=ri, right_residue=r1)
        left_node = tree.find(tree.root, r1)
        right_node = tree.find(tree.root, ri)
        build_residue_tree(tree, left_node, mod_value, depth//2, r_list)
        build_residue_tree(tree, right_node, mod_value, depth//2, r_list)

    except ValueError:
        return
    
def build_tree(alpha, degree, q1):
    r_list = []
    tree = ResidueTree(alpha)
    build_residue_tree(tree, tree.root, q1, degree, r_list)
    return tree, r_list
