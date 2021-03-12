import csv23
import numpy as np
import re
from Bio.PDB import *

aa_index = [['A', 0], ['C', 1], ['D', 2], ['E', 3], ['F', 4], ['G', 5], ['H', 6], ['I', 7], ['K', 8], ['L', 9],
            ['M', 10], ['N', 11], ['P', 12], ['Q', 13], ['R', 14], ['S', 15], ['T', 16], ['V', 17], ['W', 18],
            ['Y', 19]]


def read_matrices():
    reader = csv23.open_reader('CCIM.csv')
    matrix = []
    with reader as reader:
        for row in reader:
            matrix.append(row)
    global charge_matrix
    charge_matrix = np.array(matrix, dtype=float)

    reader = csv23.open_reader('HCIM.csv')
    matrix = []
    with reader as reader:
        for row in reader:
            matrix.append(row)
    global hydro_matrix
    hydro_matrix = np.array(matrix, dtype=float)

    reader = csv23.open_reader('SCIM.csv')
    matrix = []
    with reader as reader:
        for row in reader:
            matrix.append(row)
    global size_matrix
    size_matrix = np.array(matrix, dtype=float)


def read_pdb():
    parser = PDBParser()
    global structure
    structure = parser.get_structure('HK885 RR426', '3DGE.pdb')
    print(structure)


def aa_min_distance():
    residues = structure.get_residues()
    res_list = Selection.unfold_entities(structure, "R")
    interface_residues = [260, 263, 266, 267, 268, 270, 271, 272, 274, 275, 278, 279, 282, 283, 287, 290, 291, 294, 298,
                          387, 437, 438, 439, 494, 495, 497, 498, 501, 502, 505, 536, 537, 538, 565, 586, 587, 588, 589,
                          590, 591]

    aa_atoms = [['ALA', ['CA', 'CB']], ['CYS', ['CA', 'CB', 'SG']], ['ASP', ['CA', 'CB', 'CG', 'OD1', 'OD2']],
                ['GLU', ['CA', 'CB', 'CG', 'CD', 'OE1', 'OE2']], ['GLY', ['CA']], ['SER', ['CA', 'CB', 'OG']],
                ['ARG', ['CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2']], ['LYS', ['CA', 'CB', 'CG', 'CD', 'CE', 'NZ']],
                ['PHE', ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']],
                ['TYR', ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']],
                ['ASN', ['CA', 'CB', 'CG', 'OD1', 'ND2']],
                ['GLN', ['CA', 'CB', 'CG', 'CD', 'OE1', 'NE2']],
                ['HIS', ['CA', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']],
                ['VAL', ['CA', 'CB', 'CG1', 'CG2']],
                ['TRP', ['CA', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']],
                ['PRO', ['CA', 'CB', 'CG', 'CD']],
                ['MET', ['CA', 'CB', 'CG', 'SD', 'CE']],
                ['LEU', ['CA', 'CB', 'CG', 'CD1', 'CD2']],
                ['THR', ['CA', 'CB', 'OG1', 'CG2']],
                ['ILE', ['CA', 'CB', 'CG1', 'CG2', 'CD1']]]

    residues = [r for r in structure.get_residues()]

    for i in range(0, 22):
        for j in range(23, len(interface_residues)):

            x = interface_residues[i] - 245
            y = (interface_residues[j] + 38) - 245

            # alpha_carbon_1 = residues[x]['CA']
            # alpha_carbon_2 = residues[y]['CA']

            aa_name_1 = str(residues[x])[9:12]
            aa_name_2 = str(residues[y])[9:12]

            # distance = alpha_carbon_1 - alpha_carbon_2
            min_distance =  40
            for k1 in range(len(aa_atoms)):
                if aa_atoms[k1][0] == aa_name_1:
                    # print(aa_name_1, aa_atoms[k1][0], len(aa_atoms[k1][1]))

                    for k2 in range(len(aa_atoms)):
                        if aa_atoms[k2][0] == aa_name_2:
                            # print(aa_name_2, aa_atoms[k2][0], len(aa_atoms[k2][1]))

                            for i2 in range(len(aa_atoms[k1][1])):
                                for j2 in range(len(aa_atoms[k2][1])):

                                    alpha_carbon_1 = residues[x][aa_atoms[k1][1][i2]]
                                    alpha_carbon_2 = residues[y][aa_atoms[k2][1][j2]]

                                    distance = alpha_carbon_1 - alpha_carbon_2
                                    # print(distance)
                                    # print('{} {} {} - {} {} {} = {}'.format(aa_name_1, x + 245, aa_atoms[k1][1][i2],
                                    #                                         aa_name_2, ((y - 38) + 245) - 481,
                                    #                                         aa_atoms[k2][1][j2],
                                    #                                         np.linalg.norm(
                                    #                                             alpha_carbon_1 - alpha_carbon_2)))
                                    # print(f'Minimum distance {min_distance} > distance {distance}')
                                    if min_distance > distance:
                                        # print('yes')
                                        min_distance = distance
                                        # print(min_distance)

                                        min_aa_1 = aa_name_1
                                        min_aa_pos_1 = x + 245
                                        min_aa_atom_1 = aa_atoms[k1][1][i2]

                                        min_aa_2 = aa_name_2
                                        min_aa_pos_2 = ((y - 38) + 245) - 481
                                        min_aa_atom_2 = aa_atoms[k2][1][j2]

                                        # print('{} {} {} - {} {} {} = {}'.format(aa_name_1, x + 245, aa_atoms[k1][1][i2],
                                        #                                         aa_name_2, ((y - 38) + 245) - 481,
                                        #                                         aa_atoms[k2][1][j2],
                                        #                                         np.linalg.norm(
                                        #                                             alpha_carbon_1 - alpha_carbon_2)))
                            if min_distance < 6:
                                print('MIN {} {} {} - {} {} {} = {}'.format(min_aa_1, min_aa_pos_1, min_aa_atom_1,
                                                                  min_aa_2, min_aa_pos_2, min_aa_atom_2,
                                                                  min_distance))
                                # line = 'MIN {} {} {} - {} {} {} = {} \n '.format(min_aa_1, min_aa_pos_1, min_aa_atom_1,
                                #                                   min_aa_2, min_aa_pos_2, min_aa_atom_2,
                                #                                   min_distance)
                                # f = open("interface neighbors.txt", "a")
                                # f.write(line)
                                # f.close()


def aa_min_distances():
    residues = structure.get_residues()
    res_list = Selection.unfold_entities(structure, "R")
    interface_residues = [260, 263, 266, 267, 268, 270, 271, 272, 274, 275, 278, 279, 282, 283, 287, 290, 291, 294, 298,
                          387, 437, 438, 439, 494, 495, 497, 498, 501, 502, 505, 536, 537, 538, 565, 586, 587, 588, 589,
                          590, 591]

    aa_atoms = [['ALA', ['CA', 'CB']], ['CYS', ['CA', 'CB', 'SG']], ['ASP', ['CA', 'CB', 'CG', 'OD1', 'OD2']],
                ['GLU', ['CA', 'CB', 'CG', 'CD', 'OE1', 'OE2']], ['GLY', ['CA']], ['SER', ['CA', 'CB', 'OG']],
                ['ARG', ['CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2']], ['LYS', ['CA', 'CB', 'CG', 'CD', 'CE', 'NZ']],
                ['PHE', ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']],
                ['TYR', ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']],
                ['ASN', ['CA', 'CB', 'CG', 'OD1', 'ND2']],
                ['GLN', ['CA', 'CB', 'CG', 'CD', 'OE1', 'NE2']],
                ['HIS', ['CA', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']],
                ['VAL', ['CA', 'CB', 'CG1', 'CG2']],
                ['TRP', ['CA', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']],
                ['PRO', ['CA', 'CB', 'CG', 'CD']],
                ['MET', ['CA', 'CB', 'CG', 'SD', 'CE']],
                ['LEU', ['CA', 'CB', 'CG', 'CD1', 'CD2']],
                ['THR', ['CA', 'CB', 'OG1', 'CG2']],
                ['ILE', ['CA', 'CB', 'CG1', 'CG2', 'CD1']]]

    residues = [r for r in structure.get_residues()]

    counter_x = 1
    print(len(interface_residues))

    for i in range(0, len(interface_residues)):

        counter_y = 1

        if i <= 22:
            a = - 245
        else:
            a = 38 - 245

        for j in range(0, len(interface_residues)):

            if j <= 22:
                b = - 245
                c = 0
            else:
                b = - 245
                c = 38

            x = interface_residues[i] + a
            y = (interface_residues[j] + c) + b

            aa_name_1 = str(residues[x])[9:12]
            aa_name_2 = str(residues[y])[9:12]

            min_distance = 150
            for k1 in range(len(aa_atoms)):
                if aa_atoms[k1][0] == aa_name_1:
                    # print(aa_name_1, aa_atoms[k1][0], len(aa_atoms[k1][1]), 'x', counter_x, 'y', counter_y)

                    for k2 in range(len(aa_atoms)):
                        if aa_atoms[k2][0] == aa_name_2:
                            # print(aa_name_2, aa_atoms[k2][0], len(aa_atoms[k2][1]))

                            for i2 in range(len(aa_atoms[k1][1])):
                                for j2 in range(len(aa_atoms[k2][1])):

                                    alpha_carbon_1 = residues[x][aa_atoms[k1][1][i2]]
                                    alpha_carbon_2 = residues[y][aa_atoms[k2][1][j2]]

                                    distance = alpha_carbon_1 - alpha_carbon_2

                                    if i <= 22:
                                        a_2 = 245

                                    else:
                                        a_2 = - 38 + 245 - 481

                                    if j <= 22:
                                        b_2 = 245
                                    else:
                                        b_2 = - 38 + 245 - 481

                                    # print('{} {} {} - {} {} {} = {}'.format(aa_name_1, x + a_2, aa_atoms[k1][1][i2],
                                    #                                         aa_name_2, y + b_2,
                                    #                                         aa_atoms[k2][1][j2],
                                    #                                         np.linalg.norm(
                                    #                                             alpha_carbon_1 - alpha_carbon_2)))

                                    if min_distance > distance:
                                        # print('yes')
                                        min_distance = distance
                                        # print(min_distance)

                                        rounding = (str(min_distance))[0:7]

                            print('{},{},{}'.format(counter_x, counter_y, rounding))
                            line = '{},{},{}\n'.format(counter_x, counter_y, float(rounding))
                            f = open("interface distances.csv", "a")
                            f.write(line)
                            f.close()
            # print('x', counter_x, 'y', counter_y)
            # print('{},{},{}'.format(counter_x, counter_y, rounding))

            counter_y += 1

        counter_x += 1


def aa_input ():
    print('Enter a pair of amino acid')
    global answ
    answ = str(input())

    verifier = re.fullmatch('[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]', answ)
    # print(verifier)

    while verifier == None:
        print('Invalid pair.')
        print('Enter a valid pair')
        answ = str(input())
        verifier = re.fullmatch('[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]', answ)


def aa_compatibility():

    for i in range(len(aa_index)):
        if aa_index[i][0] == answ[0]:
            x = aa_index[i][1]
        if aa_index[i][0] == answ[1]:
            y = aa_index[i][1]

    global pair_compatibility
    pair_compatibility = []

    matrices_values = [size_matrix[x][y], hydro_matrix[x][y], charge_matrix[x][y]]
    pair_compatibility.append(matrices_values)
    print(pair_compatibility)


def candidate_pairs():
    # print(len(size_matrix),len(hydro_matrix),len(charge_matrix))
    global candidate_pair_index
    candidate_pair_index = []
    for i in range(len(size_matrix)):

        for j in range(len(size_matrix)):
            counter = 0

            # print('size',size_matrix[i][j], 'hydropacithy', hydro_matrix[i][j])
            if size_matrix[i][j] >= pair_compatibility[0][0]:
                # print('entered', aa_index[i][0], aa_index[j][0])
                size_value = size_matrix[i][j]
                counter +=1
                # print('counter', counter)

                if hydro_matrix[i][j] >= pair_compatibility[0][1]:
                    # print('entered', aa_index[i][0], aa_index[j][0])
                    hydro_value = hydro_matrix[i][j]
                    counter +=1
                    # print('counter', counter)

                    if charge_matrix[i][j] >= pair_compatibility[0][2]:
                        # print('entered', aa_index[i][0], aa_index[j][0])
                        charge_value = charge_matrix[i][j]
                        counter +=1
                        # print('counter', counter)

                    else:
                        continue

                else:
                    continue

            else:
                continue

            if counter == 3:
                a_value = size_value
                b_value = hydro_value
                c_value = charge_value
                values = [a_value, b_value, c_value, aa_index[i][0]+aa_index[j][0]]
                candidate_pair_index.append(values)
                counter = 0

    print(candidate_pair_index)
    print(len(candidate_pair_index))


def plot():
    from matplotlib \
    import pyplot
    from mpl_toolkits.mplot3d import Axes3D

    sequence_containing_x_vals = []
    sequence_containing_y_vals = []
    sequence_containing_z_vals = []
    pair_labels = []

    for i in range(len(candidate_pair_index)):
        sequence_containing_x_vals.append(candidate_pair_index[i][0])
        sequence_containing_y_vals.append(candidate_pair_index[i][1])
        sequence_containing_z_vals.append(candidate_pair_index[i][2])
        pair_labels.append(candidate_pair_index[i][3])

    fig = pyplot.figure()
    ax = Axes3D(fig)

    ax.set_xlabel('Size')
    ax.set_ylabel('Hydropacity')
    ax.set_zlabel('Charge')

    label_list = []
    position_list = []
    for j in range(len(pair_labels)):

        label_position = pair_labels[j]
        label_list.append(label_position)
        for k in range(len(label_list)):
                x = label_list[k][1] + label_list[k][0]
                # print(x, label_list[k])
                if x == label_list[k]:
                    # print('yup')
                    pos = 0.27
                    position_list.append(pos)
                else:
                    pos = -0.05
                    position_list.append(pos)

    for jj in range(len(pair_labels)):
        ax.text(sequence_containing_x_vals[jj], sequence_containing_y_vals[jj], sequence_containing_z_vals[jj] - position_list[jj],
                pair_labels[jj], color='red', fontsize=7)


    ax.scatter(sequence_containing_x_vals, sequence_containing_y_vals, sequence_containing_z_vals)
    pyplot.show()


if __name__ == '__main__':
    # read_matrices()
    # aa_input()
    # aa_compatibility()
    # candidate_pairs()
    # plot()
    read_pdb()
    aa_min_distances()