import os
import json
import copy
import numpy as np

class Cell_Node(object):

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
        self._ct = 0

    @property
    def name(self):
        return self._name

    @property
    def ancestors(self):
        return copy.deepcopy(self._ancestors)

    @property
    def level(self):
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
        return copy.deepcopy(self._list_of_children)

    @property
    def ct(self):
        return self._ct

    def add_ct(self, ii):
        self._ct += ii


def _build_tree(node_in, dendrogram_out, ancestors_in=None, level=0):
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
    _build_tree(node_in, dendrogram_out,
                ancestors_in=None, level=0)


if __name__ == "__main__":
    fname = 'dendrogram/dend.json'
    assert os.path.isfile(fname)
    with open(fname, 'r') as in_file:
        dendrogram = json.load(in_file)

    tree = {}
    build_tree(dendrogram, tree)

    for node in tree.values():
        ct = 0
        if len(node.children)==0:
            continue
        for cc in node.children:
            ct += tree[cc].ct
        assert(ct==node.ct)

    max_level = -1
    for node in tree.values():
        if node.level > max_level:
            max_level = node.level

    leaf_list = []
    for node in tree.values():
        if len(node.children) == 0:
            leaf_list.append(node)

    for level in range(max_level-1):
        out_dir = 'mouse_clusters/level%.3d' % level
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        for node in tree.values():
            if node.level != level:
                continue
            root_node = node
            out_name = os.path.join(out_dir, '%s.txt' % root_node.name)
            leaves = []
            leaf_dexes = []
            for leaf in leaf_list:
                if leaf.descended_from(root_node.name) or leaf.name==root_node.name:
                    ct_file_name = 'mouse_clusters/leaves/%s.txt' % leaf.name
                    with open(ct_file_name, 'r') as in_file:
                        census = in_file.readline().strip().split()
                        dex = int(census[2])
                        leaves.append(leaf.name)
                        leaf_dexes.append(dex)
            with open(out_name, 'w') as out_file:
                out_file.write('# n_cells %d\n' % node.ct)
                for n in leaves:
                    out_file.write('%s ' % n)
                out_file.write('\n')
                for d in leaf_dexes:
                    out_file.write('%d ' % d)
                out_file.write('\n')
                for cc in root_node.children:
                    out_file.write('%s -- %d\n' % (cc, tree[cc].ct))
            print('got %d -- %s' % (node.ct, out_name))

    out_dir = 'mouse_clusters/level999'
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    for leaf in leaf_list:
        out_name = os.path.join(out_dir,'%s.txt' % leaf.name)
        with open(out_name, 'w') as out_file:
            out_file.write('# n_cells %d\n' % leaf.ct)
            out_file.write('%s\n' % leaf.name)

    #### list example pairs ####
    parity = []
    two_to_one = []
    worse = []
    for node in tree.values():
        if len(node.children) == 0:
            continue
        child_list = node.children
        n_children = len(child_list)
        #print(node.name,n_children)
        for ic1, _c1 in enumerate(child_list):
            c1 = tree[_c1]
            n1 = np.log10(c1.ct)
            f1 = 'level%.3d/%s' % (c1.level, c1.name)
            for _c2 in child_list[ic1+1:]:
                c2 = tree[_c2]
                n2 = np.log10(c2.ct)
                f2 = 'level%.3d/%s' % (c2.level, c2.name)
                out_tuple = (f1, f2, c1.ct, c2.ct,
                             len(c1.children), len(c2.children))
                if np.abs(n1-n2)<0.2:
                    parity.append(out_tuple)
                elif np.abs(n1-n2)<0.6:
                    two_to_one.append(out_tuple)
                else:
                    worse.append(out_tuple)

    print('parity %d' % len(parity))
    print('two to one %d' % len(two_to_one))
    print('worse %d' % len(worse))
    #print(parity)
    print(len(tree))

    out_dir = 'mouse_clusters/examples'
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    name_list = ['parity.txt', 'two_to_one.txt', 'worse.txt']
    for name, data in zip(name_list, [parity, two_to_one, worse]):
        with open(os.path.join(out_dir, name), 'w') as out_file:
            out_file.write("# name_1 name_2 ct_1 ct_2 children_1 children_2\n")
            for pair in data:
                out_file.write('%s %s %d %d %d %d\n' %
                               (pair[0], pair[1],
                                pair[2], pair[3],
                                pair[4], pair[5]))
