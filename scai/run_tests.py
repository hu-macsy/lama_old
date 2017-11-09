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


SERIAL_TESTS = [
    Test('commonTest', [ 'common/test/commonTest' ], is_boost_test=True),
    Test('loggingTest', [ 'logging/test/test.sh' ], is_boost_test=False),
    Test('tracingTest', [ 'tracing/test/test.sh' ], is_boost_test=False),
    Test('taskingTest', [ 'tasking/test/taskingTest' ], is_boost_test=True),
    Test('kregistryTest', [ 'kregistry/test/kregistryTest' ], is_boost_test=True),
    Test('hmemoTest', [ 'hmemo/test/hmemoTest' ], is_boost_test=True),
    Test('blaskernelTest', [ 'blaskernel/test/blaskernelTest' ], is_boost_test=True),
    Test('utilskernelTest', [ 'utilskernel/test/utilskernelTest' ], is_boost_test=True),
    Test('sparsekernelTest', [ 'sparsekernel/test/sparsekernelTest' ], is_boost_test=True),
    Test('dmemoTest', [ 'dmemo/test/dmemoTest' ], is_boost_test=True),
    Test('lamaTest', [ 'lama/test/lamaTest' ], is_boost_test=True),
    Test('lamaStorageTest', [ 'lama/test/storage/lamaStorageTest' ], is_boost_test=True),
    Test('lamaMatrixTest', [ 'lama/test/matrix/lamaMatrixTest' ], is_boost_test=True),
    Test('partitioningTest', [ 'partitioning/test/partitioningTest' ], is_boost_test=True),
    Test('solverTest', [ 'solver/test/solverTest' ], is_boost_test=True)
]

PASSED = " [ PASSED ] "
FAILED = " [ FAILED ] "


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


def run_test(test, output_dir):
    # Note: we are implicitly assuming unique test names here
    stdout_path = os.path.join(output_dir, "{}_stdout.txt".format(test.name))
    stderr_path = os.path.join(output_dir, "{}_stderr.txt".format(test.name))

    args = deepcopy(test.args)

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


def run_serial_tests(tests, output_dir):
    failed_tests = []
    successful_tests = []
    status_messages = [ "  Running test #{} {} ... ".format(index + 1, test.name)
                       for (index, test) in enumerate(tests) ]
    longest_message = max(status_messages, key=len)

    print("Running {} serial tests ...".format(len(tests)))

    for (message, test) in zip(status_messages, tests):
        padding = ' ' * ( len(longest_message) - len(message) )
        print(message + padding, end="")
        sys.stdout.flush()

        try:
            passed = run_test(test, output_dir)
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
            .format(result_label, len(tests), len(successful_tests), len(failed_tests)))

    print("Generating test report...")
    html_report_path = combine_reports(output_dir)
    print("Report successfully generated. Path to report:\n")
    print(url_for_file_path(html_report_path))
    print()

    return len(failed_tests) == 0


def main():
    parser = argparse.ArgumentParser(description='Run LAMA tests.')
    parser.add_argument('--output_dir', dest='output_dir', required=True,
                        help='The directory in which to store the standard output, test logs and test reports.')
    args = parser.parse_args()

    print("Running LAMA tests.")
    print("Output from individual tests (logs, reports, stdout/stderr) will be stored in the provided output directory.")
    print("Output directory: {}\n".format(args.output_dir))
    serial_result = run_serial_tests(SERIAL_TESTS, args.output_dir)
    print()

    return 0 if serial_result else 1

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
