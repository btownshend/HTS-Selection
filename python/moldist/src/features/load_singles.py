import csv
import os
from rdkit.Chem import AllChem


def load_singles():
    fold = []  # hits[i][j][2] is CI for aptamer j in the presence of target i
    targets = []  # (names,SMILES) of the targets
    aptamers = []  # Naseqs of the aptamers
    with open(os.getenv("DATA") + '/Hits/singles.csv') as csv_file:  # From matlab:   TRPSummary.dumpsingles()
        csv_reader = csv.reader(csv_file, delimiter='\t')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                # print(f'Column names are {", ".join(row)}')
                aptamers = row[2::2]
                line_count += 1
            else:
                line_count += 1
                fold.append([[float(x) for x in row[2::2]], [float(x) for x in row[3::2]]])
                mol = AllChem.MolFromSmiles(row[1])
                targets.append((row[0],mol))
        print(f'Processed {line_count - 1} targets.')
    return fold, aptamers, targets
