import os
import json
import copy
import numpy as np

import argparse

class Cell_Node(object):
    """
    A class to store all of the information (children, ancestors, name) related
    to a node on the dendrogram
    """

    def __init__(self, name, ancestors, level):
        """
        Parameters
        ----------
        name -- the name of the current node

        ancestors -- a list of the names of the ancestors of this node

        level -- the level (0=the common ancestor of the whole tree) of this node
        """
        self._name = name
        self._ancestors = copy.deepcopy(ancestors)
        if ancestors is None:
            self._ancestor_set = None
        else:
            self._ancestor_set = set(ancestors)
        self._level = level
        self._list_of_children = []

    @property
    def name(self):
        """
        The cell_set_accession of this node
        """
        return self._name

    @property
    def ancestors(self):
        """
        A list of the cell_set_accessions of the ancestors of this node
        """
        return copy.deepcopy(self._ancestors)

    @property
    def level(self):
        """
        The level of this node in the tree (level=0 is the common ancestor
        of the whole tree)
        """
        return self._level

    def descended_from(self, name):
        """
        Return a boolean answering the question "is this node descended
        from the node specified by 'name'?"
        """
        return name in self._ancestor_set

    def add_child(self, child):
        """
        Add a node's name to the list of this node's children
        """
        self._list_of_children.append(child)

    @property
    def children(self):
        """
        A list of the cell_set_accessions of all nodes descended from this tree
        """
        return copy.deepcopy(self._list_of_children)


def _build_tree(node_in, dendrogram_out, ancestors_in=None, level=0):
    """
    Add a node to the dendrogram

    node_in -- the node to be added to the tree, as realized in the original
    dict stored in dend.json

    dendgrogram_out -- the dict in which the tree is being constructed

    ancestors_in -- a list of the names of ancestors of this node (None if
    there are no ancestors)

    level -- the level of this node in the tree (level=0 means the common
    ancestor of the tree)
    """
    is_leaf = False
    if 'leaf_attributes' in node_in.keys():
        attr_key = 'leaf_attributes'
        is_leaf = True
    else:
        attr_key = 'node_attributes'
    name = node_in[attr_key][0]['cell_set_accession']
    tree_node = Cell_Node(name, ancestors_in, level)

    dendrogram_out[name] = tree_node
    level += 1

    if ancestors_in is None:
        ancestors = [name]
    else:
        dendrogram_out[ancestors_in[-1]].add_child(name)
        ancestors = copy.deepcopy(ancestors_in)
        ancestors.append(name)
    del tree_node

    if not is_leaf:
        for child in node_in['children']:
            _build_tree(child, dendrogram_out,
                        ancestors_in=ancestors, level=level)

def build_tree(dendrogram_in, dendrogram_out):
    """
    dendrogram_in is the dict that results from running

        import json
        with open('dend.json','r') as in_file:
            dendrogram_in = json.load(in_file)

    dendrogram_out is an empty dict into which the new tree will be written

    dendrogram_out is a dict whose keys are the cell_set_accession names of the
    nodes in dend.json. It's values are instances of the class CellNode defined
    above. Each CellNode has a list children and a list ancestors containing
    the cell_set_accession of the children and ancestors of that node
    """
    _build_tree(dendrogram_in, dendrogram_out,
                ancestors_in=None, level=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dend_name', type=str, default=None,
                        help='path to dend.json file')
    args = parser.parse_args()
    with open(args.dend_name, 'r') as in_file:
        raw_dendrogram = json.load(in_file)
    parsed_dendrogram = {}
    build_tree(raw_dendrogram, parsed_dendrogram)
