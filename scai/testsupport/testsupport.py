###
 # @file testsupport/testsupport.py
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
 # @brief Helper functions for the main, Python-based test runner.
 # @author Andreas Longva
 # @date 15.11.2017
###
from __future__ import print_function
import urlparse, urllib
import errno
import os
import subprocess
import shutil
import sys
from glob import glob


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
TESTSUPPORT_DIR = SCRIPT_DIR


class colors:
    PASS = '\033[92m'
    FAIL = '\033[91m'
    WARNING = '\033[93m'
    NOCOLOR = '\033[0m'


def warning(message):
    print(colors.WARNING + "WARNING: " + message + colors.NOCOLOR)
    sys.stdout.flush()


def ensure_directory_exists(dir_path):
    try:
        os.makedirs(dir_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def url_for_file_path(file_path):
    abs_path = os.path.abspath(file_path)
    return urlparse.urljoin('file:', urllib.pathname2url(abs_path))


def combine_reports(output_dir):
    report_glob = os.path.join(output_dir, '*_report.xml')
    log_glob = os.path.join(output_dir, '*_log.xml')
    report_path = os.path.join(output_dir, 'report.xml')
    with open(report_path, 'w') as output_file:
        output_file.write('<TestOutput>')
        for report_file in glob(report_glob):
            with open(report_file, 'r') as input_file:
                shutil.copyfileobj(input_file, output_file)
        for log_file in glob(log_glob):
            with open(log_file, 'r') as input_file:
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
