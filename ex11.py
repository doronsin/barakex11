import itertools


class Node:
    def __init__(self, data, pos=None, neg=None):
        self.data = data
        self.positive_child = pos
        self.negative_child = neg

    def is_leaf(self):
        return self.positive_child is None


class Record:
    def __init__(self, illness, symptoms):
        self.illness = illness
        self.symptoms = symptoms


def parse_data(filepath):
    with open(filepath) as data_file:
        records = []
        for line in data_file:
            words = line.strip().split()
            records.append(Record(words[0], words[1:]))
        return records


class Diagnoser:
    def __init__(self, root):
        self.root = root

    ###################################### DIAGNOSE ########################################################################
    def diagnose_helper(self, root, symptoms):
        '''
        This recursive function receives the root of a tree (type: Node) and a list of symptoms (list of strings)
        and return the data in the leaf - i.e the diagnosed illness (string)
        '''

        if root.positive_child == None:  # if this is a 'leaf' return its data. i.e. illness
            return root.data
        if root.data in symptoms:  # if the data is in the list of symptoms then it means that we need to go deeper
            root = root.positive_child
        else:
            root = root.negative_child
        return self.diagnose_helper(root, symptoms)  # go recursively with the 'root'

    def diagnose(self, symptoms):
        '''
        This function receives a list of symptoms and returns the diagnosis
        :param symptoms: list of strings
        :return: the name of the illness (string)
        '''
        return self.diagnose_helper(self.root, symptoms)  # calling a recursuve helper with self.root)

    ############################ SUCCESS RATE ##############################################################################

    def calculate_success_rate(self, records):
        '''
        This function gets a list of records and return the diagnostic success rate (type: float)
        '''
        num_of_records = len(records)
        num_of_success = 0

        for record in records:
            # if the tree diagnoses according to the record, count a success
            if self.diagnose(record.symptoms) == record.illness:
                num_of_success += 1
        return num_of_success / num_of_records

    ############################### ALL ILLNESSES #########################################################################

    def dict_leaf(self, list_leaf):
        '''
        This is a 'util' function of all_illnesses function.
        :param list_leaf: a list of strings representing the values of the tree leaf
        :return: a dictionary of keys which represent the values of the tree leafs and the number of their occurrences
         in the tree.
        '''
        my_dict = {k: list_leaf.count(k) for k in list_leaf}
        return my_dict

    def dict_to_list(self, dict_leaf):
        '''
        This is a 'util' function of all_illnesses function. It receives a dictionary of illnesses and their occurrences,
        (dict_leaf) and turn them to a list of tuples, sort them according to the number of occurrences and creates a sorted list of
        illnesses
        :return: A sorted list of illnesses, from the most common to the most rare.
        '''
        tuple_list = list(dict_leaf.items())  # turns the dic. into a list of tuples
        tuple_list.sort(key=lambda x: x[1])  # sort the tuple list according to the second value of each tuple
        sorted_list = tuple_list[::-1]  # sort the tuple list from the highest to the lowest
        # now we isolate the first item from the tuples (i.e. the illness names) and creats a new list out pf the,
        new_list = []
        for i in sorted_list:
            new_list.append(i[0])
        return new_list

    def all_illnesses_helper(self, root):
        '''
        This is a recursive helper that reaches all the leafs of a tree and returns a list of their values (string)
        '''
        if not root.positive_child:  # if you are in a leaf then return its data.
            return [root.data]
        # else, run a recurtion both to the right and to the left branch and return the sum of them - a list of strings
        return self.all_illnesses_helper(root.positive_child) + self.all_illnesses_helper(root.negative_child)

    def all_illnesses(self):
        '''
        This function returns a list of all the illnesses in the tree sorted according to the number of their occurrence
        '''
        root = self.root
        list_leaf = self.all_illnesses_helper(root)
        dict_leaf = self.dict_leaf(list_leaf)
        return self.dict_to_list(dict_leaf)

    ############################### MOST RARW ILLNESS ######################################################################

    def most_rare_illness(self, records):
        '''
        This function gets a list og records and creates a diagnosis list using a tree. Then it returns the most rare
        illness diagnosed
        '''
        all_illnesses = self.all_illnesses()  # contains all illnesses in the tree
        for record in records:  # run over al the records
            all_illnesses.append(self.diagnose(record.symptoms))  # add diagnosis using the diagnose function
        # after this for-loop, all illness appears 1 or more times in the list "all_illnesses", and therefore if there
        # is a illness which is not diagnosed from the records, it would be selected as the most rare one (because it
        # appears 1 time exactly in all_illnesses
        dict_diag = self.dict_leaf(all_illnesses)  # creates a dictionary of occurrences from the list
        tuple_list = list(dict_diag.items())  # turn the dictionaly to the tuple list
        return min(tuple_list, key=lambda x: x[1])[0]  # find the key with lowest occurrence value and return the key

    ################################ PATH TO ILLNESS #######################################################################

    def path_to_illness_helper(self, root, illness, path_sublist):
        '''
        This function gets the current node of a tree and the illness as well as the path to the cuurent Node,
        and returns the sub paths
        '''
        if not root.positive_child:
            if root.data == illness:
                return [path_sublist]
            else:
                return []
        else:
            # returns all the paths in the positive direction
            pos_paths = self.path_to_illness_helper(root.positive_child, illness, path_sublist + [True])

            # returns all the paths in the negative directions
            neg_paths = self.path_to_illness_helper(root.negative_child, illness, path_sublist + [False])

            # sum all the paths from each direction.
            return pos_paths + neg_paths

    def paths_to_illness(self, illness):
        '''
        This function gets a illness name and return a list of all the paths to it.
        It calls a helper with the tree root to search and an empty list
        :param illness: a string
        :return: a list of boolean expressions
        '''

        return self.path_to_illness_helper(self.root, illness, [])


############################################### Build_Tree #############################################################

def place_illness(records, symptoms_on_path, symptoms):
    '''
    This is a helper function of build_tree function. It find the write illness to place on a leaf and returns the
    leaf Node
    :param symptoms_on_path: A list of symptoms (strings) which we went through until we got to the leaf
    :param symptoms: A list of symptoms that the tree has in all its branches
    :return: Node of illness
    '''
    dict_ill = {}  # create a dictionary of illnesses from the records
    for record in records:  # run through the records
        illness = record.illness
        rec_symptoms = record.symptoms

        # these are the criteria for the illnesses 'candidates'.
        if set(symptoms_on_path).issubset(rec_symptoms) and \
                (set(rec_symptoms) - set(symptoms_on_path)).intersection(set(symptoms)) == set():
            dict_ill[illness] = dict_ill.get(illness, 0) + 1  # if you pass the criteria then add it to the list
    if dict_ill == {}:  # no candidate found just return the current record's illness
        return Node(records[0].illness, None, None)
    else:  # if there are several candidates, then return the one with the max amount of occurrences in the list
        fre_ill = max(list(dict_ill.items()), key=lambda x: x[1])[0]
        return Node(fre_ill, None, None)


def build_tree_helper(symptoms, counter, records, path_sublist):
    '''
    This is a recursive helper for the build_tree function, who is responsible for building the tree.
    :param symptoms: list of strings
    :param counter: counter that responsible that the depth of the tree will be the number of the symptoms
    :param path_sublist: is a list of symptoms that we went through until we got to the leaf
    :return: the 'rood' (Node)
    '''
    if len(symptoms) == counter:
        return place_illness(records, path_sublist, symptoms)
    else:
        return Node(symptoms[counter],
                    build_tree_helper(symptoms, counter + 1, records, path_sublist + [symptoms[counter]]),
                    build_tree_helper(symptoms, counter + 1, records, path_sublist))


def build_tree(records, symptoms):
    '''
    This function creates a tree out of a list of symptoms and records
    :param records: list of 'record' objects
    :param symptoms: a list of string
    :return: the tree root (Node)
    '''
    counter = 0
    tree_node = build_tree_helper(symptoms, counter, records, [])
    return tree_node


############################################# OPTIMAL TREE #############################################################

def optimal_tree(records, symptoms, depth):
    '''
    This function creates many diagnostic trees with a specific depth and returns the one with the mose successful
    diagnostic success rate
    :param records:
    :param symptoms:
    :param depth: the number of symptoms the tree includes
    :return: a tree (root's Node)
    '''
    opt_tree = None
    opt_tree_rate = 0

    # creates a list of subsets of symptoms in the length of 'depth'
    sub_sets = list(itertools.combinations(symptoms, depth))
    for i in range(len(sub_sets)):  # for each subset
        temp_tree = build_tree(records, sub_sets[i])  # create a diagnostic tree
        temp_diagnoser = Diagnoser(temp_tree)  # create a diagnoser from the tree
        curr_rate = temp_diagnoser.calculate_success_rate(records)  # calculate it success rate

        # save the the tree with the best rate
        if curr_rate > opt_tree_rate:
            opt_tree_rate = curr_rate
            opt_tree = temp_tree
    return opt_tree


if __name__ == "__main__":

    # Manually build a simple tree.
    #                cough
    #          Yes /       \ No
    #        fever           healthy
    #   Yes /     \ No
    # influenza   cold

    flu_leaf = Node("influenza", None, None)
    l_cold_leaf = Node("cold", None, None)
    r_cold_leaf = Node("cold", None, None)
    healthy_leaf = Node("healthy", None, None)
    left_vertex = Node("headache", flu_leaf, l_cold_leaf)
    right_vertex = Node("headache", r_cold_leaf, healthy_leaf)
    root = Node("cough", left_vertex, right_vertex)

    diagnoser = Diagnoser(root)

    # Simple test
    diagnosis = diagnoser.diagnose(["cough", "fever"])
    if diagnosis == "influenza":
        print("Test passed")
    else:
        print("Test failed. Should have printed cold, printed: ", diagnosis)

# #simple test 2
# records = [Record('cold', ['coughh']), Record('cold', ['cough', 'headache']), Record('influenza', ['cough', 'fever'])]
# rate = diagnoser.calculate_success_rate(records)
# print(rate)
#
# #simple test 3
# print(diagnoser.all_illnesses())
#
# #simple test 4
# print(diagnoser.most_rare_illness(records))
#
# #simple test 5
# print(diagnoser.paths_to_illness("cold"))
#
# #simple test 6
# x= build_tree(records, ['cough', 'headache'])
#
#
# # simple tet 7
# records = [Record('flu',['cough', 'fever']), Record('cold', ['cough'])]
# x = optimal_tree(records,['cough','fever'],1)
# pass
