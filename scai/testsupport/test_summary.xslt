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
            padding-right: 0.5em;
            padding-left: 0.5em;
        }

        .summary td {
            border-top: 1px solid white;
        }

        td.number {
            text-align: right;
        }

        .test_suite {
            padding: 20px;
            margin-top: 20px;
            margin-bottom: 20px;
            border: 1px solid black;
        }

        .test_case {
            padding: 20px;
            margin-top: 20px;
            margin-bottom: 20px;
            border: 1px solid black;
        }

        .suite_anchor, .suite_anchor:visited {
            text-decoration: none;
            color: white;
        }

        tr.passed {
            background-color: #5cb85c;
            color: white;
        }

        tr.passed-link:hover {
            background-color: #92D092;
            cursor: pointer;
        }

        tr.failed {
            background-color: #d9534f;
            color: white;
        }

        tr.failed-link:hover {
            background-color: #E68D78;
            cursor: pointer;
        }

    </style>
    </head>
    <body>
        <h1>Test summary</h1>
        <h4>Top-level test suite results</h4>
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
                    <xsl:choose>
                        <xsl:when test="@result = 'passed'">
                            <xsl:attribute name="class">passed passed-link</xsl:attribute>
                        </xsl:when>
                        <xsl:otherwise>
                            <xsl:attribute name="class">failed failed-link</xsl:attribute>
                        </xsl:otherwise>
                    </xsl:choose>
                    <xsl:attribute name="onclick">
                        <xsl:text>window.document.location='#</xsl:text>
                        <xsl:apply-templates select="current()" mode="fullName"/>
                        <xsl:text>'</xsl:text>
                    </xsl:attribute>

                    <td><xsl:value-of select="@name"/></td>
                    <td><xsl:value-of select="@result"/></td>
                    <td class="number"><xsl:value-of select="@test_cases_passed"/></td>
                    <td class="number"><xsl:value-of select="@test_cases_failed"/></td>
                </tr>
            </xsl:for-each>
            </tbody>
        </table>

        <h1>Test suite details</h1>
        <xsl:apply-templates select="//TestSuite[ancestor::TestResult]"/>

        <h1>Test suite logs</h1>
        <xsl:apply-templates select="//TestCase[ancestor::TestLog]">
            <!-- Make sure most relevant messages (i.e. errors) come first -->
            <xsl:sort select="Error" order="descending"/>
            <xsl:sort select="Warning" order="descending"/>
            <xsl:sort select="Message" order="descending"/>
        </xsl:apply-templates>
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
    <div class="test_suite" style="clear: both;">
        <xsl:attribute name="id"><xsl:apply-templates select="current()" mode="fullName"/></xsl:attribute>
        <table class="summary" style="clear: both;">
            <tr>
                <td>Suite name</td>
                <td><xsl:value-of select="@name"/></td>
            </tr>
            <tr>
                <td>Full name</td>
                <td><xsl:apply-templates select="current()" mode="fullName"/></td>
            </tr>
            <tr>
                <td>Status</td>
                <td><xsl:value-of select="@result"/></td>
            </tr>
        </table>
        <div>
            <h4>Summary</h4>
            <table class="summary">
                <tr>
                    <td>Passed tests</td>
                    <td class="number"><xsl:value-of select="@test_cases_passed"/></td>
                </tr>
                <tr>
                    <td>Failed tests</td>
                    <td class="number"><xsl:value-of select="@test_cases_failed"/></td>
                </tr>
                <tr>
                    <td>Assertions failed</td>
                    <td class="number"><xsl:value-of select="@assertions_failed"/></td>
                </tr>
                <tr>
                    <td>Warnings failed</td>
                    <td class="number"><xsl:value-of select="@warnings_failed"/></td>
                </tr>
                <tr>
                    <td>Expected failures</td>
                    <td class="number"><xsl:value-of select="@expected_failures"/></td>
                </tr>
                <tr>
                    <td>Passed with warning</td>
                    <td class="number"><xsl:value-of select="@test_cases_passed_with_warnings"/></td>
                </tr>
                <tr>
                    <td>Skipped tests</td>
                    <td class="number"><xsl:value-of select="@test_cases_skipped"/></td>
                </tr>
                <tr>
                    <td>Aborted tests</td>
                    <td class="number"><xsl:value-of select="@test_cases_aborted"/></td>
                </tr>
            </table>
        </div>

        <xsl:if test="TestSuite">
            <div>
                <h4>Test suites</h4>
                <table class="summary">
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
                                <xsl:choose>
                                    <xsl:when test="@result = 'passed'">
                                        <xsl:attribute name="class">passed passed-link</xsl:attribute>
                                    </xsl:when>
                                    <xsl:otherwise>
                                        <xsl:attribute name="class">failed failed-link</xsl:attribute>
                                    </xsl:otherwise>
                                </xsl:choose>
                                <xsl:attribute name="onclick">
                                    <xsl:text>window.document.location='#</xsl:text>
                                    <xsl:apply-templates select="current()" mode="fullName"/>
                                    <xsl:text>'</xsl:text>
                                </xsl:attribute>

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
            </div>
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
                                    <xsl:choose>
                                        <xsl:when test="@result = 'passed'">
                                            <xsl:attribute name="class">passed passed-link</xsl:attribute>
                                        </xsl:when>
                                        <xsl:otherwise>
                                            <xsl:attribute name="class">failed failed-link</xsl:attribute>
                                        </xsl:otherwise>
                                    </xsl:choose>
                                    <xsl:attribute name="onclick">
                                        <xsl:text>window.document.location='#</xsl:text>
                                        <xsl:apply-templates select="current()" mode="fullName"/>
                                        <xsl:text>'</xsl:text>
                                    </xsl:attribute>

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

<xsl:template match="TestCase[ancestor::TestLog]">
    <div class="test_case">
        <xsl:attribute name="id"><xsl:apply-templates select="current()" mode="fullName"/></xsl:attribute>

        <h4><xsl:apply-templates select="." mode="fullName"/></h4>

        <table class="summary">
            <thead>

            </thead>
                <th>File</th>
                <th>Line</th>
                <th>Type</th>
                <th>Log message</th>
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
</xsl:template>

</xsl:stylesheet>
