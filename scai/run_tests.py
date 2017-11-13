from __future__ import print_function

import sys
import argparse
import os
import os.path
import errno
import random
import subprocess
import collections
import shutil
import urlparse, urllib
from glob import glob
from copy import deepcopy


Test = collections.namedtuple('Test', ['name', 'args', 'is_boost_test'])

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
TESTSUPPORT_DIR = os.path.join(SCRIPT_DIR, 'testsupport')

class colors:
    PASS = '\033[92m'
    FAIL = '\033[91m'
    NOCOLOR = '\033[0m'


NORMAL_TESTS = [
    Test('commonTest', [ 'common/test/commonTest' ], is_boost_test=True),
    Test('loggingTest', [ 'logging/test/test.sh' ], is_boost_test=False),
    Test('tracingTest', [ 'tracing/test/test.sh' ], is_boost_test=False),
    Test('taskingTest', [ 'tasking/test/taskingTest' ], is_boost_test=True),
    Test('kregistryTest', [ 'kregistry/test/kregistryTest' ], is_boost_test=True),
    Test('hmemoTest', [ 'hmemo/test/hmemoTest' ], is_boost_test=True),
    Test('blaskernelTest', [ 'blaskernel/test/blaskernelTest' ], is_boost_test=True),
    Test('utilskernelTest', [ 'utilskernel/test/utilskernelTest' ], is_boost_test=True),
    Test('sparsekernelTest', [ 'sparsekernel/test/sparsekernelTest' ], is_boost_test=True),
]

MPI_TESTS = [
    Test('dmemoTest', [ 'dmemo/test/dmemoTest' ], is_boost_test=True),
    Test('lamaTest', [ 'lama/test/lamaTest' ], is_boost_test=True),
    Test('lamaStorageTest', [ 'lama/test/storage/lamaStorageTest' ], is_boost_test=True),
    Test('lamaMatrixTest', [ 'lama/test/matrix/lamaMatrixTest' ], is_boost_test=True),
    Test('partitioningTest', [ 'partitioning/test/partitioningTest' ], is_boost_test=True),
    Test('solverTest', [ 'solver/test/solverTest' ], is_boost_test=True)
]

PASSED = " " + colors.PASS + "[ PASSED ]" + colors.NOCOLOR + " "
FAILED = " " + colors.FAIL + "[ FAILED ]" + colors.NOCOLOR + " "

def url_for_file_path(file_path):
    abs_path = os.path.abspath(file_path)
    return urlparse.urljoin('file:', urllib.pathname2url(abs_path))


def combine_reports(output_dir):
    report_glob = os.path.join(output_dir, '*_report.xml')
    report_path = os.path.join(output_dir, 'report.xml')
    with open(report_path, 'w') as output_file:
        output_file.write('<TestOutput>')
        for report_file in glob(report_glob):
            with open(report_file, 'r') as input_file:
                shutil.copyfileobj(input_file, output_file)
        output_file.write('</TestOutput>')

    stylesheet_path = os.path.join(TESTSUPPORT_DIR, 'test_summary.xslt')
    report_html_path = os.path.join(output_dir, 'report.html')

    xsltproc_args = [
        "xsltproc",
        "-o",
        report_html_path,
        stylesheet_path,
        report_path
    ]
    subprocess.check_call(xsltproc_args)
    return report_html_path

def ensure_directory_exists(dir_path):
    try:
        os.makedirs(dir_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def run_test(test, output_dir, prepend_args = []):
    # Note: we are implicitly assuming unique test names here
    stdout_path = os.path.join(output_dir, "{}_stdout.txt".format(test.name))
    stderr_path = os.path.join(output_dir, "{}_stderr.txt".format(test.name))

    args = prepend_args + deepcopy(test.args)

    if test.is_boost_test:
        args +=  [
            "--report_level=detailed",
            "--log_level=message",
            "--output_format=XML",
            "--output_dir={}".format(output_dir)
        ]

    ensure_directory_exists(output_dir)

    with open(stdout_path, 'w') as stdout, open(stderr_path, 'w') as stderr:
        retcode = subprocess.call(args, stdout=stdout, stderr=stderr)
        return retcode == 0


def run_tests(tests, output_dir, prepend_args = []):
    failed_tests = []
    successful_tests = []
    status_messages = [ "  Running test #{} {} ... ".format(index + 1, test.name)
                       for (index, test) in enumerate(tests) ]
    longest_message = max(status_messages or [ '' ], key=len)

    for (message, test) in zip(status_messages, tests):
        padding = ' ' * ( len(longest_message) - len(message) )
        print(message + padding, end="")
        sys.stdout.flush()

        try:
            passed = run_test(test, output_dir, prepend_args=prepend_args)
        except Exception as e:
            passed = False;
            print("\nException when running test {}:\n{}\n".format(test.name, e))

        if passed:
            successful_tests.append(test)
            print(PASSED)
        else:
            failed_tests.append(test)
            print(FAILED)

    return (successful_tests, failed_tests)


def run_mpi_tests(tests, output_dir, np):
    failed_tests = []
    passed_tests = []

    for n in np:
        print("Running {} MPI tests ({} processors) ...".format(len(MPI_TESTS), n))
        mpi_args = [ "mpirun", "-np", str(n), "--tag-output" ]
        (passed, failed) = run_tests(tests, output_dir, prepend_args=mpi_args)
        print()
        passed_tests += passed
        failed_tests += failed
    return (passed_tests, failed_tests)


def test_exists(name):
    all_tests = NORMAL_TESTS + MPI_TESTS
    all_test_names = [ test.name for test in all_tests ]
    return name in all_test_names


def main():
    parser = argparse.ArgumentParser(description='Run LAMA tests.')
    parser.add_argument('--output_dir', dest='output_dir', required=True,
                        help='The directory in which to store the standard output, test logs and test reports.')
    parser.add_argument('--mpi', dest='mpi', action='store_true',
                        help='Whether or not to use MPI for MPI-enabled tests.')
    parser.add_argument('--np', dest='np', type=int, nargs='+', required=False, default=[ 1 ],
                        help='The number of processors to use for MPI. Can supply multiple arguments (e.g. --np 1 2 3)')
    parser.add_argument('--tests', dest='tests', nargs='+', type=str, required=False,default=None,
                        help='A list of tests to run. If not specified, all tests are run.')
    args = parser.parse_args()
    output_dir = args.output_dir

    print("Running LAMA tests.")
    print("Output from individual tests (logs, reports, stdout/stderr) will be stored in the provided output directory.")
    print("Output directory: {}\n".format(output_dir))

    if args.tests:
        for name in args.tests:
            if not test_exists(name):
                print("ERROR: User requested to run test '{}', which does not exist. Aborting...".format(name))
                return -1

        print("Note: Only running the following requested tests:\n  {}\n".format(", ".join(args.tests)))

    # If MPI is not requested, run "mpi tests" as normal tests (e.g. without mpirun).
    normal_tests = deepcopy(NORMAL_TESTS) if args.mpi else NORMAL_TESTS + MPI_TESTS
    filtered_normal_tests = [ test for test in normal_tests if test.name in args.tests ] if args.tests else normal_tests
    filtered_mpi_tests = [ test for test in MPI_TESTS if test.name in args.tests ] if args.tests else MPI_TESTS

    print("Running {} normal tests ...".format(len(filtered_normal_tests)))
    (normal_passed, normal_failed) = run_tests(filtered_normal_tests, output_dir)
    print()

    mpi_passed = []
    mpi_failed = []
    if args.mpi:
        (mpi_passed, mpi_failed) = run_mpi_tests(filtered_mpi_tests, output_dir, args.np)
    else:
        print("MPI tests not requested. Skipping ...")

    all_tests = normal_passed + normal_failed + mpi_passed + mpi_failed
    passed = normal_passed + mpi_passed
    failed = normal_failed + mpi_failed
    all_tests_passed = len(failed) == 0
    result_label = PASSED if all_tests_passed else FAILED
    print()
    print("Tests completed. Results:")
    print("{}    {} tests run. {} successful tests. {} failed tests."
            .format(result_label, len(all_tests), len(passed), len(failed)))

    print()
    print("Generating test report...")
    html_report_path = combine_reports(output_dir)
    print("Report successfully generated. URL:\n")
    print(url_for_file_path(html_report_path))
    print()

    return len(failed)

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
