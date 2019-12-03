import sys


# upper percentage points of distribution for alpha = 0.05
dist_points = [3.84, 5.99, 7.81, 9.49, 11.1, 12.6, 14.1, 15.5, 16.9, 18.3,
               19.7, 21.0, 22.4, 23.7, 25.0, 26.3, 27.6, 28.9, 30.1, 31.4]


def read_param_from_file(path):
    try:
        with open(path) as file:
            data = [int(x) for x in file.read().split(' ')]
            if len(data) != 4:
                print('Incorrect data.\nTemplate: [multiplier] [start value] [increment] [modulus]')
                sys.exit(1)
        return data
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


def linear_congruential_generator(a, x0, c, N):
    sequence = [x0]
    for i in range(N - 1):
        sequence.append((a * sequence[i] + c) % N)
    return sequence


def get_period(sequence):
    max_period = 0
    for i in range(len(sequence) - 1):
        try:
            max_period = max(max_period, sequence[i+1:].index(sequence[i]) + 1)
        except ValueError:
            continue
    return max_period


def Pearsons_chi_squared_test(sequence):
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
    return x2, x2 < dist_points[len(groups) - 1]    # (len(groups) - 1) is degrees of freedom


if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print('Usage: lcg.py [path] Input.txt [path] Output.txt (optional)')
        sys.exit(1)

    # read parameters from file
    params = read_param_from_file(sys.argv[1])

    # generate sequence
    seq = linear_congruential_generator(params[0], params[1], params[2], params[3])

    # print result
    result = [seq, get_period(seq), Pearsons_chi_squared_test(seq)]
    print('Sequence: {}\nPeriod: {}\nUniform distribution: {}'.format(result[0], result[1], result[2]))

    # save result
    output_path = 'output.txt' if len(sys.argv) == 2 else sys.argv[2]
    save_result_to_file(output_path, result)
    print('Saved result to ' + output_path)
