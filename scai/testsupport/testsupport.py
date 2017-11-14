import urlparse, urllib
import errno
import os
import subprocess
import shutil
from glob import glob


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
TESTSUPPORT_DIR = SCRIPT_DIR


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
