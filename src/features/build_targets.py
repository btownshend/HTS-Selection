# Set expected fold change of molecules as a function of the aptamer
import csv
import math


def build_fold(tested):
    fold = []  # fold[i][j] is the fold change for aptamer j in the presence of target i
    targets = []  # Names of the targets in the form, for example,  91A2  (PlateRowCol)
    aptamers = []  # Names of the aptamers
    with open('../data/raw/fold.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                aptamers = row[1:]
                line_count += 1
            else:
                line_count += 1
                fold.append([float(x) if math.isfinite(float(x)) else 1.0 for x in row[1:]])
                targets.append(row[0])
        print(f'Processed {line_count-1} targets.')
    print(targets)

    # Setup ML output as y
    y = []
    nofold = [1.0 for _ in fold[0]]
    found=[]
    for t in tested:
        name = t.GetProp("NAME")
        if name in targets:
            y.append(fold[targets.index(name)])
            found.append(name)
        else:
            y.append(nofold)
    if len(found) != len(fold):
        missing = [x for x in targets if x not in found]
        print("No match against molecules for: ", missing)
    return y, aptamers


def build_hit(tested):
    # Set categorization of molecules
    # Property, HIT, will be set to true iff we found an aptamer for the molecule
    allhits = []
    with open('../data/raw/hasaptamer.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                # print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                # print(f'\tPlate {row[0]}, well {row[1]} is a hit.')
                line_count += 1
                allhits.append("%s%s" % (row[0], row[1]))
        print(f'Processed {line_count - 1} hits.')
    # Setup ML output as y
    y = []
    found = []
    for t in tested:
        name = t.GetProp("NAME")
        if name in allhits:
            y.append(1.0)
            found.append(name)
        else:
            y.append(0.0)
    print("Matched %d hits" % len(found))
    if len(found) != len(allhits):
        missing = [x for x in allhits if x not in found]
        print("No match against molecules for: ", missing)
    return y
