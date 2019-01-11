###
 # @file run_tests.py
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief Main test runner for LAMA.
 # @author Andreas Longva
 # @date 02.11.2017
###
from __future__ import print_function

import sys
import argparse
import os
import os.path
import random
import subprocess
import collections
import time
import tempfile
from copy import deepcopy
from testsupport import testsupport
from testsupport.testsupport import colors


Test = collections.namedtuple('Test', ['name', 'args', 'is_boost_test'])


NORMAL_TESTS = [
    Test( 'commonTest', [ 'common/test/commonTest' ], is_boost_test=True ),
    Test( 'loggingTest', [ 'logging/test/test.sh' ], is_boost_test=False ),
    Test( 'tracingTest', [ 'tracing/test/test.sh' ], is_boost_test=False ),
    Test( 'taskingTest', [ 'tasking/test/taskingTest' ], is_boost_test=True ),
    Test( 'kregistryTest', [ 'kregistry/test/kregistryTest' ], is_boost_test=True ),
    Test( 'hmemoTest', [ 'hmemo/test/hmemoTest' ], is_boost_test=True ),
    Test( 'blaskernelTest', [ 'blaskernel/test/blaskernelTest' ], is_boost_test=True ),
    Test( 'utilskernelTest', [ 'utilskernel/test/utilskernelTest' ], is_boost_test=True ),
    Test( 'sparsekernelTest', [ 'sparsekernel/test/sparsekernelTest' ], is_boost_test=True ),
    Test( 'lamaStorageTest', [ 'lama/test/storage/lamaStorageTest' ], is_boost_test=True )
]

# ToDo: Test( 'ipbclsTest', [ 'ipbcls/test/lsbctest' ], is_boost_test=False )

GPU_ONLY_TESTS = [
    Test( 'commonCUDATest', [ 'common/test/cuda/commonCUDATest' ], is_boost_test=True ),
    Test( 'taskingCUDATest', [ 'tasking/test/cuda/taskingCUDATest' ], is_boost_test=True ),
    Test( 'hmemoCUDATest', [ 'hmemo/test/cuda/hmemoCUDATest' ], is_boost_test=True )
]

MPI_TESTS = [
    Test( 'dmemoTest', [ 'dmemo/test/dmemoTest' ], is_boost_test=True ),
    Test( 'lamaTest', [ 'lama/test/lamaTest' ], is_boost_test=True ),
    Test( 'lamaMatrixTest', [ 'lama/test/matrix/lamaMatrixTest' ], is_boost_test=True ),
    Test( 'partitioningTest', [ 'partitioning/test/partitioningTest' ], is_boost_test=True ),
    Test( 'solverTest', [ 'solver/test/solverTest' ], is_boost_test=True )
]

PASSED = " " + colors.PASS + "[ PASSED ]" + colors.NOCOLOR + " "
FAILED = " " + colors.FAIL + "[ FAILED ]" + colors.NOCOLOR + " "


def run_test(test, output_dir, temp_dir = None, prepend_args = [], tag = ""):
    # Note: we are implicitly assuming unique test names here
    stdout_path = os.path.join(output_dir, "{}{}_stdout.txt".format(test.name, tag))
    stderr_path = os.path.join(output_dir, "{}{}_stderr.txt".format(test.name, tag))

    args = prepend_args + deepcopy(test.args)

    if test.is_boost_test:
        args +=  [
            "--report_level=detailed",
            "--log_level=test_suite",
            "--output_format=XML",
            "--output_dir={}".format(output_dir)
        ]
        if temp_dir:
            args += [ "--temp_dir={}".format(temp_dir) ]

    testsupport.ensure_directory_exists(output_dir)

    if temp_dir:
        testsupport.ensure_directory_exists(temp_dir)

    with open(stdout_path, 'w') as stdout, open(stderr_path, 'w') as stderr:
        retcode = subprocess.call(args, stdout=stdout, stderr=stderr)
        return retcode == 0


def run_tests(tests, output_dir, temp_dir = None, prepend_args = [], tag = ""):
    failed_tests = []
    successful_tests = []
    status_messages = [ "  Running test #{} {} ... ".format(index + 1, test.name)
                       for (index, test) in enumerate(tests) ]
    longest_message = max(status_messages or [ '' ], key=len)

    for (message, test) in zip(status_messages, tests):
        padding = ' ' * ( len(longest_message) - len(message) )
        print(message + padding, end="")
        sys.stdout.flush()

        start = time.time()

        try:
            passed = run_test(test, output_dir, temp_dir=temp_dir, prepend_args=prepend_args, tag=tag)
        except Exception as e:
            passed = False;
            print("\nException when running test {}:\n{}\n".format(test.name, e))

        duration = time.time() - start

        if passed:
            successful_tests.append(test)
            print(PASSED, end="")
        else:
            failed_tests.append(test)
            print(FAILED, end="")

        print("( {:.2f} s )".format(duration))

    return (successful_tests, failed_tests)


def run_mpi_tests(tests, output_dir, np, intel, temp_dir=None):
    failed_tests = []
    passed_tests = []

    for n in np:
        print("Running {} MPI tests ({} processors) ...".format(len(tests), n))
        tag = "_mpi_{}".format(str(n))
        print(tag)
        if intel:
            mpi_args = [ "mpirun", "-np", str(n), "-l" ]
        else:
            mpi_args = [ "mpirun", "-np", str(n), "--tag-output" ]
        (passed, failed) = run_tests(tests, output_dir, temp_dir=temp_dir, prepend_args=mpi_args, tag=tag)
        print()
        passed_tests += passed
        failed_tests += failed
    return (passed_tests, failed_tests)

def run_gpu_tests(tests, output_dir, temp_dir=None):
    (passed_tests, failed_tests) = run_tests(tests, output_dir, temp_dir=temp_dir)
    return (passed_tests, failed_tests)

def test_exists(name):
    all_tests = NORMAL_TESTS + MPI_TESTS
    all_test_names = [ test.name for test in all_tests ]
    return name in all_test_names


def main():
    parser = argparse.ArgumentParser(description='Run LAMA tests.')
    parser.add_argument('--output_dir', dest='output_dir', default=None,
                        help='The directory in which to store the standard output, test logs and test reports.')
    parser.add_argument('--temp_dir', dest='temp_dir', default=None,
                        help='The directory in which to store temporary files used by tests.')
    parser.add_argument('--gpu', dest='gpu', action='store_true',
                        help='Whether or not to run tests on GPU device.')
    parser.add_argument('--mpi', dest='mpi', action='store_true',
                        help='Whether or not to use MPI for MPI-enabled tests.')
    parser.add_argument('--intel', dest='intel', action='store_true',
                        help='To be set when Intel MPI is used')
    parser.add_argument('--np', dest='np', type=int, nargs='+', required=False, default=[ 1 ],
                        help='The number of processors to use for MPI. Can supply multiple arguments (e.g. --np 1 2 3)')
    parser.add_argument('--tests', dest='tests', nargs='+', type=str, required=False,default=None,
                        help='A list of tests to run. If not specified, all tests are run.')
    args = parser.parse_args()

    # This needs to happen before we assign to output_dir, otherwise tempfile will
    # cause the warning to be displayed.
    if args.output_dir and os.path.isdir(args.output_dir):
        testsupport.warning("Output directory already exists. Note that this may cause the test report to contain "
                            "information from old tests. You may wish to delete the output directory before continuing.")

    output_dir = args.output_dir if args.output_dir else tempfile.mkdtemp(suffix='_lama_testoutput')

    print("Running LAMA tests.")
    print("Output from individual tests (logs, reports, stdout/stderr) will be stored in the following output directory.")
    print("Output directory:          {}".format(output_dir))
    print("Temporary files directory: {}\n".format(args.temp_dir if args.temp_dir else "Not specified. Determined by test binary."))
    print()

    #  For running tests on GPU we set the environment variable for context to CUDA
    #  Do not use cusparse library as entries in CSR data are not always sorted yet
    if args.gpu:
        os.environ[ "SCAI_CONTEXT" ] = "CUDA"
        os.environ[ "SCAI_CUDA_USE_CUSPARSE" ] = "0"

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
    (normal_passed, normal_failed) = run_tests(filtered_normal_tests, output_dir, temp_dir=args.temp_dir)
    print()

    mpi_passed = []
    mpi_failed = []
    if args.mpi:
        (mpi_passed, mpi_failed) = run_mpi_tests(filtered_mpi_tests, output_dir, args.np, args.intel, temp_dir=args.temp_dir)
    else:
        print("MPI tests not requested. Skipping ...")

    gpu_passed = []
    gpu_failed = []
    if args.gpu:
        (gpu_passed, gpu_failed ) = run_gpu_tests(GPU_ONLY_TESTS, output_dir, temp_dir=args.temp_dir)
    else:
        print("GPU/Cuda tests not requested. Skipping ...")

    all_tests = normal_passed + normal_failed + mpi_passed + mpi_failed + gpu_passed + gpu_failed
    passed = normal_passed + mpi_passed + gpu_passed
    failed = normal_failed + mpi_failed + gpu_failed
    all_tests_passed = len(failed) == 0
    result_label = PASSED if all_tests_passed else FAILED
    print()
    print("Tests completed. Results:")
    print("{}    {} tests run. {} successful tests. {} failed tests."
            .format(result_label, len(all_tests), len(passed), len(failed)))

    print()
    print("Generating test report...")
    html_report_path = testsupport.combine_reports(output_dir)
    print("Report successfully generated. URL:\n")
    print(testsupport.url_for_file_path(html_report_path))
    print()

    return len(failed)

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
