<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="html" indent="yes" />

<xsl:template match="TestOutput">
    <html>
    <head>
    <title>Test results</title>
    <style>
        body {
            font-family: "Helvetica";
        }

        table.summary {
            text-align: left;
            border-collapse:collapse;
        }

        .summary td, .summary th {
            padding: 0.3em;
            padding-right: 0.7em;
            padding-left: 0.7em;
        }

        .summary td {
            border-top: 1px solid white;
        }

        td.number {
            text-align: right;
        }

        .test_summary {
            padding: 20px;
            margin-top: 20px;
            margin-bottom: 20px;
            border: 1px solid black;
        }

        .test_summary {
            padding: 20px;
            margin-top: 20px;
            margin-bottom: 20px;
            border: 1px solid black;
        }

        .section {
            margin-top: 30px;
        }

        .summary_header {
            margin-top: 0;
        }

        .suite_anchor, .suite_anchor:visited, #nav-pane a, #nav-pane a:visited {
            text-decoration: none;
            color: #6489c4;
        }

        .passed {
            background-color: #5cb85c;
            color: white;
        }

        .passed-link:hover {
            background-color: #92D092;
            cursor: pointer;
        }

        .failed {
            background-color: #d9534f;
            color: white;
        }

        .failed-link:hover {
            background-color: #E68D78;
            cursor: pointer;
        }

        .indicator {
            width: 1em;
        }

        table.log_table tr:nth-child(odd) td {
            background-color: #eeeeee;
        }

        #nav-pane {
            position: fixed;
            width: 200px;
            height: 100%;
        }

        #nav-pane ul {
            list-style-type: none;
            margin: 0.2em;
            margin-left: 0.4em;
            padding: 0;
        }

        #main-content {
            margin-left: 220px;
        }

    </style>
    </head>
    <body>
        <div id="nav-pane">
            <h4>Quick navigation</h4>
            <ul>
                <li><a href="#test_summary">Top-level summary</a></li>
                <li><a href="#test_suite_details">Test suite details</a></li>
                <li><a href="#test_case_logs">Test case logs</a></li>
            </ul>
        </div>

        <div id="main-content">
            <h1 id="test_summary">Test summary</h1>
            <div class="test_summary">
                <h2 class="summary_header">Top-level test suite results</h2>
                <table class="summary">
                <thead>
                    <tr>
                        <th>Test suite name</th>
                        <th>Status</th>
                        <th>Passed tests</th>
                        <th>Failed tests</th>
                    </tr>
                </thead>
                <tbody>
                    <xsl:for-each select="TestResult/TestSuite">
                        <xsl:sort select="@result"/>
                        <xsl:sort select="translate(@name, 'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ')" order="ascending" />
                        <tr>
                            <xsl:call-template name="applyPassedOrFailed"/>
                            <xsl:call-template name="fullNameOnClick"/>

                            <td><xsl:value-of select="@name"/></td>
                            <td><xsl:value-of select="@result"/></td>
                            <td class="number"><xsl:value-of select="@test_cases_passed"/></td>
                            <td class="number"><xsl:value-of select="@test_cases_failed"/></td>
                        </tr>
                    </xsl:for-each>
                    </tbody>
                </table>
            </div>

            <h1 id="test_suite_details">Test suite details</h1>
            <xsl:apply-templates select="//TestSuite[ancestor::TestResult]"/>

            <h1 id="test_case_logs">Test suite logs</h1>
            <xsl:apply-templates select="//TestCase[ancestor::TestLog]">
                <!-- Make sure most relevant messages (i.e. errors) come first -->
                <xsl:sort select="Error" order="descending"/>
                <xsl:sort select="Warning" order="descending"/>
                <xsl:sort select="Message" order="descending"/>
            </xsl:apply-templates>
        </div>
    </body>
    </html>
</xsl:template>

<xsl:template match="*" mode="fullName">
    <xsl:if test="../@name">
        <xsl:apply-templates select=".." mode="fullName"/><xsl:text> :: </xsl:text>
    </xsl:if>
    <xsl:value-of select="@name"/>
</xsl:template>

<xsl:template match="TestSuite[ancestor::TestResult]">
    <div class="test_summary" style="clear: both;">
        <xsl:attribute name="id"><xsl:apply-templates select="current()" mode="fullName"/></xsl:attribute>
        <h2 class="summary_header">Test suite summary</h2>
        <table class="summary section">
            <tr>
                <td class="indicator" rowspan="4">
                    <xsl:choose>
                        <xsl:when test="@result = 'passed'">
                            <xsl:attribute name="class">indicator passed</xsl:attribute>
                        </xsl:when>
                        <xsl:otherwise>
                            <xsl:attribute name="class">indicator failed</xsl:attribute>
                        </xsl:otherwise>
                    </xsl:choose>
                </td>
                <td>Suite name</td>
                <td><xsl:value-of select="@name"/></td>
            </tr>
            <tr>
                <td>Full name</td>
                <td><xsl:apply-templates select="current()" mode="fullName"/></td>
            </tr>
            <tr>
                <td>Parent suite</td>
                <td>
                    <a class="suite_anchor">
                        <xsl:attribute name="href">#<xsl:apply-templates select=".." mode="fullName"/></xsl:attribute>
                        <xsl:apply-templates select=".." mode="fullName"/>
                    </a>
                </td>
            </tr>
            <tr>
                <td>Status</td>
                <td><xsl:value-of select="@result"/></td>
            </tr>
        </table>
        <table class="summary section">
            <thead>
                <tr>
                    <th>Passed tests</th>
                    <th>Failed tests</th>
                    <th>Assertions failed</th>
                    <th>Warnings failed</th>
                    <th>Expected failures</th>
                    <th>Passed with warning</th>
                    <th>Skipped tests</th>
                    <th>Aborted tests</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td class="number"><xsl:value-of select="@test_cases_passed"/></td>
                    <td class="number"><xsl:value-of select="@test_cases_failed"/></td>
                    <td class="number"><xsl:value-of select="@assertions_failed"/></td>
                    <td class="number"><xsl:value-of select="@warnings_failed"/></td>
                    <td class="number"><xsl:value-of select="@expected_failures"/></td>
                    <td class="number"><xsl:value-of select="@test_cases_passed_with_warnings"/></td>
                    <td class="number"><xsl:value-of select="@test_cases_skipped"/></td>
                    <td class="number"><xsl:value-of select="@test_cases_aborted"/></td>
                </tr>
            </tbody>
        </table>

        <xsl:if test="TestSuite">
            <table class="summary section">
                <thead>
                    <th>Suite name</th>
                    <th>Status</th>
                    <th>Passed</th>
                    <th>Failed</th>
                    <th>Skipped</th>
                    <th>Aborted</th>
                </thead>
                <tbody>
                    <xsl:for-each select="TestSuite">
                        <xsl:sort select="@result"/>
                        <xsl:sort select="translate(@name, 'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ')" order="ascending" />
                        <tr>
                            <xsl:call-template name="applyPassedOrFailed"/>
                            <xsl:call-template name="fullNameOnClick"/>

                            <td><xsl:value-of select="@name"/></td>
                            <td><xsl:value-of select="@result"/></td>
                            <td class="number"><xsl:value-of select="@test_cases_passed"/></td>
                            <td class="number"><xsl:value-of select="@test_cases_failed"/></td>
                            <td class="number"><xsl:value-of select="@test_cases_skipped"/></td>
                            <td class="number"><xsl:value-of select="@test_cases_aborted"/></td>
                        </tr>
                    </xsl:for-each>
                </tbody>
            </table>
        </xsl:if>

        <xsl:if test="TestCase">
            <div style="clear: both">
                <h4>Test cases</h4>
                <table class="summary">
                    <thead>
                        <th>Test name</th>
                        <th>Status</th>
                        <th>Assertions passed</th>
                        <th>Assertions failed</th>
                        <th>Warnings</th>
                        <th>Expected failures</th>
                    </thead>
                    <tbody>
                        <xsl:for-each select="TestCase">
                            <xsl:sort select="@result"/>
                            <xsl:sort select="translate(@name, 'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ')" order="ascending" />
                                <tr>
                                    <xsl:call-template name="applyPassedOrFailed"/>
                                    <xsl:call-template name="fullNameOnClick"/>

                                    <td><xsl:value-of select="@name"/></td>
                                    <td><xsl:value-of select="@result"/></td>
                                    <td class="number"><xsl:value-of select="@assertions_passed"/></td>
                                    <td class="number"><xsl:value-of select="@assertions_failed"/></td>
                                    <td class="number"><xsl:value-of select="@warnings_failed"/></td>
                                    <td class="number"><xsl:value-of select="@expected_failures"/></td>
                                </tr>
                        </xsl:for-each>
                    </tbody>
                </table>
            </div>
        </xsl:if>
    </div>
</xsl:template>

<xsl:template name="applyPassedOrFailed">
    <xsl:choose>
        <xsl:when test="@result = 'passed'">
            <xsl:attribute name="class">passed passed-link</xsl:attribute>
        </xsl:when>
        <xsl:otherwise>
            <xsl:attribute name="class">failed failed-link</xsl:attribute>
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<xsl:template name="fullNameOnClick">
    <xsl:attribute name="onclick">
            <xsl:text>window.document.location='#</xsl:text>
            <xsl:apply-templates select="current()" mode="fullName"/>
            <xsl:text>'</xsl:text>
    </xsl:attribute>
</xsl:template>

<xsl:template match="TestCase[ancestor::TestLog]">
    <xsl:if test="Error|Warning|Message">
        <div class="test_summary">
            <xsl:attribute name="id"><xsl:apply-templates select="current()" mode="fullName"/></xsl:attribute>
            <h2 class="summary_header">Test case log summary</h2>

            <table class="summary section">
                <tr>
                    <td>Test case name</td>
                    <td><xsl:value-of select="@name"/></td>
                </tr>
                <tr>
                    <td>Full name</td>
                    <td><xsl:apply-templates select="current()" mode="fullName"/></td>
                </tr>
                <tr>
                    <td>Parent suite</td>
                    <td>
                        <a class="suite_anchor">
                            <xsl:attribute name="href">#<xsl:apply-templates select=".." mode="fullName"/></xsl:attribute>
                            <xsl:apply-templates select=".." mode="fullName"/>
                        </a>
                    </td>
                </tr>
            </table>

            <table class="summary section log_table">
                <thead>
                    <th>File</th>
                    <th>Line</th>
                    <th>Type</th>
                    <th>Log message</th>
                </thead>
                <tbody>
                    <xsl:for-each select="Error|Warning|Message">
                        <tr>
                            <td><xsl:value-of select="@file"/></td>
                            <td><xsl:value-of select="@line"/></td>
                            <td><xsl:value-of select="local-name(.)"/></td>
                            <td><xsl:value-of select="."/></td>
                        </tr>
                    </xsl:for-each>
                </tbody>
            </table>
        </div>
    </xsl:if>
</xsl:template>

</xsl:stylesheet>
