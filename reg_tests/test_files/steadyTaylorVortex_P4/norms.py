import os
import math
import os.path


def exit_if_file_does_not_exist(fname):
    if not os.path.exists(fname):
        print 'No file: ' + fname
        exit(1)


def num_lines_output(fname):
    exit_if_file_does_not_exist(fname)
    f = open(fname)
    lines = f.readlines()
    header = lines[1].split()
    j = 0
    for line in lines[2:]:
        j += 1
        if line.split()[0] == header[0]:
            break
    f.close()
    return j


def tail(fname, n):
    exit_if_file_does_not_exist(fname)
    stdin, stdout = os.popen2("tail -n " + str(n) + " " + fname)
    stdin.close()
    lines = stdout.readlines()
    stdout.close()
    return lines


def extract_norm(lines):
    l2norm = {}
    numNodes = lines[0].split()[3]

    for line in lines:
        tokens = line.split()
        l2norm[tokens[0]] = float(tokens[6])

    return numNodes, l2norm


def compute_and_check_ooas(basep_name, numResolutions, dim, min_ooas):
    first_name = basep_name + "_R0.dat"
    linesToRead = num_lines_output(first_name)
    nodes0, norm0 = extract_norm(tail(first_name, linesToRead))

    for j in xrange(numResolutions - 1):
        name = basep_name + "_R" + str(j + 1) + ".dat"
        nodes1, norm1 = extract_norm(tail(name, linesToRead))
        order = {}
        for varname, norm in norm1.iteritems():
            num_ratio = math.log(float(norm1[varname]) / float(norm0[varname]))
            den_ratio = math.log(float(nodes1) / float(nodes0))
            order[varname] = dim * num_ratio / den_ratio

        norm0 = norm1
        nodes0 = nodes1

        output_string = "R" + str(j + 1) + "-" + str(j)

        for varname, ooa in order.iteritems():
            output_string += " " + varname + ": " + str(ooa)
            if ooa > min_ooas[varname]:
                print 'Bad convergence'
                print output_string
                exit(1)


def compute_and_check_norms(base_name, numResolutions, dim, minP, maxP, min_ooas):
    for j in xrange(minP, maxP + 1, 1):
        name = base_name + "_P" + str(j)
        compute_and_check_ooas(name, numResolutions, dim, min_ooas)


dimension = 2
numResolution = 3
minP = 4
maxP = 4
base_name = 'steadyTaylorVortex'

min_allowed_ooas = {'dpdx[0]': -3.0, 'dpdx[1]': -3.0,
                    'velocity[0]': -5.0, 'velocity[1]': -5.0}

compute_and_check_norms(base_name, numResolution, dimension, minP, maxP,
                        min_allowed_ooas)

exit(0)
