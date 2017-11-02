from __future__ import print_function

import sys
import argparse
import os
import os.path
import errno
import random
import subprocess
import collections

SERIAL_TESTS = [
    'common/test/commonTest',
    'tasking/test/taskingTest',
    'kregistry/test/kregistryTest',
    'hmemo/test/hmemoTest',
    'blaskernel/test/blaskernelTest',
    'utilskernel/test/utilskernelTest',
    'sparsekernel/test/sparsekernelTest',
    'dmemo/test/dmemoTest',
    'lama/test/lamaTest',
    'lama/test/storage/lamaStorageTest',
    'lama/test/matrix/lamaMatrixTest',
    'partitioning/test/partitioningTest',
    'solver/test/solverTest'
]

PASSED = " [ PASSED ] "
FAILED = " [ FAILED ] "


Test = collections.namedtuple('Test', ['name', 'path'])


def ensure_directory_exists(dir_path):
    try:
        os.makedirs(dir_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def run_test(testname, path, output_dir):
    # Note: we are implicitly assuming unique test names here
    report_path = os.path.join(output_dir, "{}_report.txt".format(testname))
    log_path = os.path.join(output_dir, "{}_log.txt".format(testname))
    stdout_path = os.path.join(output_dir, "{}_stdout.txt".format(testname))
    stderr_path = os.path.join(output_dir, "{}_stderr.txt".format(testname))

    args = [
        path,
        "--report_level=detailed",
        "--report_format=HRF",
        "--report_sink={}".format(report_path),
        "--log_level=all",
        "--log_format=HRF",
        "--log_sink={}".format(log_path)
    ]

    ensure_directory_exists(output_dir)

    with open(stdout_path, 'w') as stdout, open(stderr_path, 'w') as stderr:
        retcode = subprocess.call(args, stdout=stdout, stderr=stderr)
        return retcode == 0


def run_serial_tests(executable_paths, output_dir):
    all_tests = [ Test(name=os.path.basename(path), path=path) for path in executable_paths ]
    failed_tests = []
    successful_tests = []
    status_messages = [ "  Running test #{} {} ... ".format(index + 1, test.name)
                       for (index, test) in enumerate(all_tests) ]
    longest_message = max(status_messages, key=len)

    print("Running {} serial tests ...".format(len(executable_paths)))

    for (message, test) in zip(status_messages, all_tests):
        padding = ' ' * ( len(longest_message) - len(message) )
        print(message + padding, end="")
        sys.stdout.flush()

        try:
            passed = run_test(test.name, test.path, output_dir)
        except Exception as e:
            passed = False;
            print("\nException when running test {}:\n{}\n".format(test.name, e))

        if passed:
            successful_tests.append(test)
            print(PASSED)
        else:
            failed_tests.append(test)
            print(FAILED)

    all_tests_passed = len(failed_tests) == 0
    result_label = PASSED if all_tests_passed else FAILED
    print("\nSerial tests completed. Results:")
    print("{}    {} tests run. {} successful tests. {} failed tests."
            .format(result_label, len(executable_paths), len(successful_tests), len(failed_tests)))

    return len(failed_tests) == 0


def main():
    parser = argparse.ArgumentParser(description='Run LAMA tests.')
    parser.add_argument('--output_dir', dest='output_dir', required=True,
                        help='The directory in which to store the standard output, test logs and test reports.')
    args = parser.parse_args()
    serial_result = run_serial_tests(SERIAL_TESTS, args.output_dir)
    print()

    return 0 if serial_result else 1

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
