from ex11 import Record, Diagnoser, Node, build_tree, optimal_tree


def check(num, expected, result, input, checker=None):
    b = expected == result
    if checker is not None:
        b = checker(expected, result)
    if b:
        print("Test "+str(num)+" passed")
    else:
        print("Test "+str(num)+" failed. for "+str(input)+" Should have printed "+str(expected)+", printed: ", result)


def list_check(lst1, lst2):
    for a in lst1:
        if a not in lst2:
            return False
    for a in lst2:
        if a not in lst1:
            return False
    return True


def tree_check(tree1, tree2):
    if tree1.data != tree2.data:
        return False
    if tree1.is_leaf() or tree2.is_leaf():
        if tree1.is_leaf() and tree2.is_leaf():
            return True
        return False
    if not tree_check(tree1.positive_child, tree2.positive_child):
        return False
    if not tree_check(tree1.negative_child, tree2.negative_child):
        return False
    return True


def print_tree(root_tree, level=1, answer=''):
    print('    ' * (level - 1) + str(answer) + '+---' * (level > 0) + root_tree.data)
    if root_tree.positive_child:
        if root_tree.positive_child:
            print_tree(root_tree.positive_child, level + 1, 'P')
        if root_tree.negative_child:
            print_tree(root_tree.negative_child, level + 1, 'N')


def parse_data(filepath):
    with open(filepath) as data_file:
        records = []
        for line in data_file:
            words = line.split()
            records.append(Record(words[0], words[1:]))
        return records


if __name__ == "__main__":

    #   Manually build a simple tree.
    #                                         cough
    #                      Yes /                            \ No
    #                   sneezing                        sneezing
    #               Yes /       \ No                Yes /       \ No
    #             fever         fever               fever         fever
    #       Yes /    \ No    Yes/     \No       Yes /    \ No    Yes/     \No
    #       dead   cold   influenza   cold      dead   influenza cold   healthy

    leaf1 = Node("dead", None, None)
    leaf2 = Node("cold", None, None)
    fever1 = Node("fever", leaf1, leaf2)

    leaf3 = Node("influenza", None, None)
    leaf4 = Node("cold", None, None)
    fever2 = Node("fever", leaf3, leaf4)

    sneezing1 = Node("sneezing", fever1, fever2)

    leaf5 = Node("dead", None, None)
    leaf6 = Node("influenza", None, None)
    fever3 = Node("fever", leaf5, leaf6)

    leaf7 = Node("cold", None, None)
    leaf8 = Node("healthy", None, None)
    fever4 = Node("fever", leaf7, leaf8)

    sneezing2 = Node("sneezing", fever3, fever4)

    root = Node("cough", sneezing1, sneezing2)

    diagnoser = Diagnoser(root)
    # Simple test
    expected = "cold"
    diagnosis = diagnoser.diagnose(["cough", "sneezing"])

    check(1, expected, diagnosis, ["cough", "sneezing"])
    # test 2
    records = []
    records.append(Record("healthy", []))
    records.append(Record("influenza", ["cough", "fever", "sneezing"]))
    records.append(Record("dead", ["cough", "fever", "sneezing"]))
    records.append(Record("dead", ["cough", "fever"]))
    records.append(Record("healthy", ["cough", "fever", "sneezing"]))
    records.append(Record("dead", ["cough", "fever", "sneezing"]))
    records.append(Record("influenza", ["cough", "fever"]))
    records.append(Record("influenza", ["cough", "fever"]))
    records.append(Record("cold", ["cough", "sneezing"]))
    records.append(Record("influenza", ["cough"]))
    records.append(Record("influenza", ["cough"]))
    records.append(Record("healthy", ["cough"]))
    records.append(Record("healthy", ["cough"]))
    records.append(Record("cold", ["cough"]))
    records.append(Record("cold", ["cough"]))
    records.append(Record("cold", ["cough"]))
    records.append(Record("dead", ["fever", "sneezing"]))
    records.append(Record("influenza", ["fever", "sneezing"]))
    records.append(Record("dead", ["fever", "sneezing"]))
    records.append(Record("influenza", ["sneezing"]))
    records.append(Record("cold", ["fever"]))
    records.append(Record("dead", ["fever"]))
    records.append(Record("dead", ["fever"]))
    records.append(Record("healthy", ["fever"]))
    records.append(Record("healthy", ["fever"]))
    records.append(Record("dead", ["fever"]))
    records.append(Record("cold", ["fever"]))
    records.append(Record("cold", ["fever"]))
    records.append(Record("cold", ["fever"]))
    records.append(Record("cold", []))
    records.append(Record("healthy", []))
    records.append(Record("cold", []))
    records.append(Record("cold", []))
    records.append(Record("healthy", []))
    records.append(Record("healthy", []))
    records.append(Record("healthy", []))

    success = diagnoser.calculate_success_rate(records)
    check(2, 20 / 36, success, "")

    # test 3
    expected = ["cold", "influenza", "dead", "healthy"]
    result = diagnoser.all_illnesses()
    check(3, expected, result, None, list_check)

    # test 4
    expected = "influenza"
    result = diagnoser.most_rare_illness(records)
    check(4, expected, result, "")

    # test 5
    expected = [[True, False, True], [False, True, False]]
    result = diagnoser.paths_to_illness("influenza")
    check(5, expected, result, "influenza", list_check)

    # test 6
    result = build_tree(records, ["cough", "sneezing", "fever"])
    check(6, root, result, "", tree_check)
    if not tree_check(result, root):
        print("expected: ")
        print_tree(root)
        print("result: ")
        print_tree(result)

    # test 7
    result = optimal_tree(records, ["cough", "sneezing", "fever"], 2)
    new_leaf1 = Node("influenza", None, None)
    new_leaf2 = Node("cold", None, None)
    new_fever1 = Node("fever", new_leaf1, new_leaf2)

    new_leaf3 = Node("dead", None, None)
    new_leaf4 = Node("healthy", None, None)
    new_fever2 = Node("fever", new_leaf3, new_leaf4)

    new_root = Node("cough", new_fever1, new_fever2)

    if not tree_check(result, new_root):
        print("expected: ")
        print_tree(new_root)
        print("result: ")
        print_tree(result)
    new_records=parse_data(r"big_data.txt")
    check(7, new_root, result, "", tree_check)
    result = optimal_tree(new_records, ["cough", "fever", "headache", "nausea", "fatigue"], 3)
    new_leaf1 = Node("influenza", None, None)
    new_leaf2 = Node("meningitis", None, None)
    new_lvl11 = Node("fatigue", new_leaf1, new_leaf2)

    new_leaf3 = Node("influenza", None, None)
    new_leaf4 = Node("strep", None, None)
    new_lvl12 = Node("fatigue", new_leaf3, new_leaf4)

    new_lvl21 = Node("headache", new_lvl11 , new_lvl12)

    new_leaf5 = Node("mono", None, None)
    new_leaf6 = Node("cold", None, None)
    new_lvl13 = Node("fatigue", new_leaf5, new_leaf6)

    new_leaf7 = Node("mono", None, None)
    new_leaf8 = Node("healthy", None, None)
    new_lvl14 = Node("fatigue", new_leaf7, new_leaf8)

    new_lvl22 = Node("headache", new_lvl13, new_lvl14)

    new_root = Node("fever", new_lvl21, new_lvl22)
    check(7.1, new_root, result, "", tree_check)

    if not tree_check(result, new_root):
        print("expected: ")
        print_tree(new_root)
        print("result: ")
        print_tree(result)
