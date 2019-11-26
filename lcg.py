import sys

param = {}
sequence = []

# upper percentage points of distribution for alpha = 0.05
dist_points = [3.84, 5.99, 7.81, 9.49, 11.1, 12.6, 14.1, 15.5, 16.9, 18.3,
               19.7, 21.0, 22.4, 23.7, 25.0, 26.3, 27.6, 28.9, 30.1, 31.4]


def read_param_from_file(path):
    try:
        with open(path) as file:
            data = [int(x) for x in file.read().split(' ')]
            if len(data) != 5:
                print('Incorrect data.\nTemplate: [modulus] [multiplier] [increment] [start value] [size]')
                sys.exit(1)
            # formula: Xn+1 = (a * Xn + c) mod m
            param['modulus'] = data[0]      # m
            param['multiplier'] = data[1]   # a
            param['increment'] = data[2]    # c
            param['start value'] = data[3]  # X1
            param['size'] = data[4]         # sequence size
    except IOError:
        print('Error opening ' + path)


def save_result_to_file(path, data):
    try:
        with open(path, 'w') as file:
            file.write('Sequence: ')
            for e in data[0]:
                file.write(str(e) + ' ')
            file.write('\nPeriod: {}\nUniform distribution: {}'.format(data[1], data[2]))
    except IOError:
        print('Error opening ' + path)


def linear_congruential_generator():
    global sequence
    sequence = [param['start value']]
    for i in range(param['size'] - 1):
        # Xn+1 = (a * Xn + c) mod m
        value = (param['multiplier'] * sequence[i] + param['increment']) % param['modulus']
        sequence.append(value)


def get_period():
    for i in range(1, len(sequence)):   # except the first element
        if sequence[i] == sequence[0]:  # find the first match with the first element
            return i


def Pearsons_chi_squared_test():
    # distribute numbers into groups
    groups = {}
    for n in sequence:
        groups[n] = groups[n] + 1 if n in groups.keys() else 1

    if len(groups) > 20:
        return 'Unknown'
    p = 1 / len(groups)     # probability
    n = len(sequence)       # numbers quantity in a sequence

    x2 = sum([pow(groups[key] / n - p, 2) for key in groups])   # calculating the test-statistic
    x2 *= n / p
    return x2 < dist_points[len(groups) - 1]    # (len(groups) - 1) is degrees of freedom


if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print('Usage: lcg.py [path] Input.txt [path] Output.txt (optional)')
        sys.exit(1)

    # read parameters from file
    read_param_from_file(sys.argv[1])

    # generate sequence
    linear_congruential_generator()

    # print result
    result = [sequence, get_period(), Pearsons_chi_squared_test()]
    print('Sequence: {}\nPeriod: {}\nUniform distribution: {}'.format(result[0], result[1], result[2]))

    # save result
    output_path = 'output.txt' if len(sys.argv) == 2 else sys.argv[2]
    save_result_to_file(output_path, result)
    print('Saved result to ' + output_path)
