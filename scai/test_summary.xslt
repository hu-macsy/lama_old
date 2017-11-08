<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="html" indent="yes" />

<xsl:template match="TestOutput">
<!--     <xsl:for-each select="TestResult"> -->
    <!--<xsl:copy>
        <xsl:apply-templates select="@*|node()" />
    </xsl:copy>-->
<!--     </xsl:for-each> -->
    <html>
    <head>
    <style>
        table.summary {
            text-align: left;
        }

        .summary td, .summary th {
            padding: 0.2em;
            padding-right: 0.5em;
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
                <xsl:sort select="@name"/>
                <tr>
                    <td><a>
                            <xsl:attribute name="href">#<xsl:apply-templates select="current()" mode="fullName"/></xsl:attribute>
                            <xsl:value-of select="@name"/>
                        </a>
                    </td>
                    <td><xsl:value-of select="@result"/></td>
                    <td class="number"><xsl:value-of select="@test_cases_passed"/></td>
                    <td class="number"><xsl:value-of select="@test_cases_failed"/></td>
                </tr>
            </xsl:for-each>
            </tbody>
        </table>

        <h1>Test suite details</h1>
        <xsl:apply-templates select="//TestSuite"/>
    </body>
    </html>
</xsl:template>

<xsl:template match="TestSuite" mode="fullName">
    <xsl:if test="not(local-name(..) = 'TestResult')">
        <xsl:apply-templates select=".." mode="fullName"/> ::
    </xsl:if>
    <xsl:value-of select="@name"/>
</xsl:template>

<xsl:template match="TestSuite">
    <div class="test_suite" style="float: left; clear: both; width: 90%">
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
        <div style="float: left">
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
            <div style="float: left">
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
                            <tr>
                                <td>
                                    <a>
                                        <xsl:attribute name="href">#<xsl:apply-templates select="current()" mode="fullName"/></xsl:attribute>
                                        <xsl:value-of select="@name"/>
                                    </a>
                                </td>
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
                        <th>Assertions passed</th>
                        <th>Assertions failed</th>
                        <th>Warnings</th>
                        <th>Expected failures</th>
                    </thead>
                    <tbody>
                        <xsl:for-each select="TestCase">
                            <xsl:sort select="@result"/>
                                <tr>
                                    <td><xsl:value-of select="@name"/></td>
                                    <td><xsl:value-of select="@result"/></td>
                                    <td class="number"><xsl:value-of select="@assertions_passed"/></td>
                                    <td class="number"><xsl:value-of select="@assertions_failed"/></td>
                                    <td class="number"><xsl:value-of select="@warnings_failed"/></td>
                                    <td class="number"><xsl:value-of select="@expcted_failures"/></td>
                                </tr>
                        </xsl:for-each>
                    </tbody>
                </table>
            </div>
        </xsl:if>
    </div>
</xsl:template>

<!--
<xsl:template match="*//TestSuite">
    <div style="clear: both;">
        <h3><xsl:value-of select="@name"/></h3>
        <div style="padding: 20px 20px 20px 20px;" id="{@name}">
            <h3>Child suites</h3>
            <div style="float: left">
                <xsl:apply-templates select="TestSuite" mode="suite_summary" />
            </div>
        <div style="float: left">
            <xsl:apply-templates select="TestCase"/>
        </div>
        </div>
    </div>
</xsl:template>
-->

<!--<xsl:template match="TestCase">
    <div><xsl:value-of select="@name"/></div>
</xsl:template>-->



</xsl:stylesheet>
