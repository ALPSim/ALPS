<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
	<xsl:template match="/">
		<html>
			<head/>
			<body>
				<xsl:for-each select="SIMULATION">
          <br/>
          <H1>ALPS Simulation</H1>
            <H2>Parameters</H2>
						<table border="1">
							<thead><tr><td>
											<B>Parameter</B>
										</td>
										<td>
											<B>Value</B>
										</td>
									</tr>
								</thead>
					  <tbody>
						<xsl:for-each select="PARAMETERS">
						  <xsl:for-each select="PARAMETER">
							  <tr>
								  <td>
									  <xsl:value-of select="@name"/>
								  </td>
								  <td>
									  <xsl:apply-templates/>
								  </td>
							  </tr>
						  </xsl:for-each>
						</xsl:for-each>
						</tbody></table>
            <H2>Execution information</H2>
						<xsl:for-each select="MCRUN">
            	<H3>Run <xsl:value-of select="position()"/></H3>
							<table border="1"><thead><tr><td><b>Executed from</b></td><td><b>to</b></td>
              	<td><b>on</b></td></tr></thead><tbody>
								<xsl:for-each select="EXECUTED">
								<tr>
									<td><xsl:for-each select="FROM"><xsl:apply-templates/></xsl:for-each></td>
									<td><xsl:for-each select="TO"><xsl:apply-templates/></xsl:for-each></td>
									<td><xsl:for-each select="MACHINE">
												<xsl:for-each select="NAME"><xsl:apply-templates/></xsl:for-each>
											</xsl:for-each></td>
                </tr>
            		</xsl:for-each>
							</tbody></table>
            </xsl:for-each>
					<H2>Averages</H2>
            <H3>Total</H3>
              <UL>
							<xsl:for-each select="AVERAGES/SCALAR_AVERAGE">
                <LI><A HREF="#scalar{position()}"><xsl:value-of select="@name"/></A></LI>
							</xsl:for-each>
							<xsl:for-each select="AVERAGES/VECTOR_AVERAGE">
                <LI><A HREF="#vector{position()}"><xsl:value-of select="@name"/></A></LI>
							</xsl:for-each>
              </UL>
				      <TABLE BORDER="1" COLS="5" WIDTH="800">
							
							<THEAD><TR><TD><B>Name</B></TD><TD><B>Count</B></TD><TD><B>Mean</B></TD>
							           <TD><B>Error</B></TD><TD><B>Tau</B></TD>
												 <TD><B>Method</B></TD></TR></THEAD>
							<TBODY>
							<xsl:for-each select="AVERAGES/SCALAR_AVERAGE">
              <A NAME="scalar{position()}"></A>
							<TR>
							  <TD><B><xsl:value-of select="@name"/></B></TD>
							  <TD><xsl:value-of select="COUNT"/></TD>
							  <TD><xsl:value-of select="MEAN"/></TD>
                              <xsl:choose>
                              <xsl:when test="ERROR/@converged = 'maybe'">
								<TD bgcolor="#ffff00"><xsl:value-of select="ERROR"/><BR/>check convergence</TD>
                              </xsl:when>
                              <xsl:when test="ERROR/@converged = 'no'">
								<TD bgcolor="#ff0000"><blink><xsl:value-of select="ERROR"/><BR/>not converged</blink></TD>
                              </xsl:when>
                              <xsl:otherwise>
								<TD ><xsl:value-of select="ERROR"/></TD>
                              </xsl:otherwise>
                              </xsl:choose>
                              <TD><xsl:value-of select="AUTOCORR"/></TD>
                              <TD><xsl:value-of select="ERROR/@method"/></TD>
							</TR>
							</xsl:for-each>
							<xsl:for-each select="AVERAGES/VECTOR_AVERAGE">
              <A NAME="vector{position()}"></A>
							<xsl:for-each select="SCALAR_AVERAGE">
							<TR>
							  <TD><B><xsl:value-of select="../@name"/></B>[<xsl:value-of 
							  select="@indexvalue"/>]</TD>
							  <TD><xsl:value-of select="COUNT"/></TD>
							  <TD><xsl:value-of select="MEAN"/></TD>
                              <xsl:choose>
                              <xsl:when test="ERROR/@converged = 'maybe'">
								<TD bgcolor="#ffff00"><xsl:value-of select="ERROR"/><BR/>check convergence</TD>
                              </xsl:when>
                              <xsl:when test="ERROR/@converged = 'no'">
								<TD bgcolor="#ff0000"><blink><xsl:value-of select="ERROR"/><BR/>not converged</blink></TD>
                              </xsl:when>
                              <xsl:otherwise>
								<TD ><xsl:value-of select="ERROR"/></TD>
                              </xsl:otherwise>
                              </xsl:choose>
								<TD><xsl:value-of select="AUTOCORR"/></TD>
								<TD><xsl:value-of select="ERROR/@method"/></TD>
							</TR>
							</xsl:for-each>
							</xsl:for-each>
							</TBODY>
							</TABLE>
						<xsl:for-each select="MCRUN">
            	<H3>Run <xsl:value-of select="position()"/></H3>
				      <TABLE BORDER="1" COLS="5" WIDTH="800">
							
							<THEAD><TR><TD><B>Name</B></TD><TD><B>Count</B></TD><TD><B>Mean</B></TD>
							           <TD><B>Error</B></TD><TD><B>Tau</B></TD>
												 <TD><B>Method</B></TD></TR></THEAD>
							<TBODY>
							<xsl:for-each select="AVERAGES/SCALAR_AVERAGE">
							<TR>
							  <TD><B><xsl:value-of select="@name"/></B></TD>
							  <TD><xsl:value-of select="COUNT"/></TD>
							  <TD><xsl:value-of select="MEAN"/></TD>
                              <xsl:choose>
                              <xsl:when test="ERROR/@converged = 'maybe'">
								<TD bgcolor="#ffff00"><xsl:value-of select="ERROR"/><BR/>check convergence</TD>
                              </xsl:when>
                              <xsl:when test="ERROR/@converged = 'no'">
								<TD bgcolor="#ff0000"><blink><xsl:value-of select="ERROR"/><BR/>not converged</blink></TD>
                              </xsl:when>
                              <xsl:otherwise>
								<TD ><xsl:value-of select="ERROR"/></TD>
                              </xsl:otherwise>
                              </xsl:choose>
								<TD><xsl:value-of select="AUTOCORR"/></TD>
								<TD><xsl:value-of select="ERROR/@method"/></TD>
							</TR>
							<xsl:for-each select="BINNED">
							<TR>
							  <TD><B><xsl:value-of select="@name"/></B></TD>
							  <TD><xsl:value-of select="COUNT"/> bins</TD>
							  <TD><xsl:value-of select="MEAN"/></TD>
                              <TD ><xsl:value-of select="ERROR"/></TD>
                <TD/>
								<TD/>
							</TR>
							</xsl:for-each>
							</xsl:for-each>
							<xsl:for-each select="AVERAGES/VECTOR_AVERAGE">
							<xsl:for-each select="SCALAR_AVERAGE">
							<TR>
							  <TD><B><xsl:value-of select="../@name"/></B>[<xsl:value-of 
							  select="@indexvalue"/>]</TD>
							  <TD><xsl:value-of select="COUNT"/></TD>
							  <TD><xsl:value-of select="MEAN"/></TD>
                              <xsl:choose>
                              <xsl:when test="ERROR/@converged = 'maybe'">
								<TD bgcolor="#ffff00"><xsl:value-of select="ERROR"/><BR/>check convergence</TD>
                              </xsl:when>
                              <xsl:when test="ERROR/@converged = 'no'">
								<TD bgcolor="#ff0000"><blink><xsl:value-of select="ERROR"/><BR/>not converged</blink></TD>
                              </xsl:when>
                              <xsl:otherwise>
								<TD ><xsl:value-of select="ERROR"/></TD>
                              </xsl:otherwise>
                              </xsl:choose>
								<TD><xsl:value-of select="AUTOCORR"/></TD>
								<TD><xsl:value-of select="ERROR/@method"/></TD>
							</TR>
							<xsl:for-each select="BINNED">
							<TR>
							  <TD><B><xsl:value-of select="@name"/></B></TD>
							  <TD><xsl:value-of select="COUNT"/> bins</TD>
							  <TD><xsl:value-of select="MEAN"/></TD>
                              <TD ><xsl:value-of select="ERROR"/></TD>
                <TD/>
								<TD/>
							</TR>
							</xsl:for-each>
							</xsl:for-each>
							</xsl:for-each>
							</TBODY>
							</TABLE>
</xsl:for-each>
          <H2>Eigenvalues</H2>
          <xsl:for-each select="EIGENVALUES">
            <xsl:for-each select="QUANTUMNUMBER">
              <B> Quantum number <xsl:value-of select="@name"/>=<xsl:value-of select="@value"/></B><BR/>
            </xsl:for-each>
            <PRE><xsl:value-of select="."/></PRE>
          </xsl:for-each>
          <H2>Eigen states</H2>
          <xsl:for-each select="EIGENSTATES">
            <H3>Sector</H3>
            <xsl:for-each select="QUANTUMNUMBER">
              <B> Quantum number <xsl:value-of select="@name"/>=<xsl:value-of select="@value"/></B><BR/>
            </xsl:for-each>
          <xsl:for-each select="EIGENSTATE">
            <H4>Eigenstate #<xsl:value-of select="position()"/></H4>
            <TABLE BORDER="1" COLS="5" WIDTH="800">
							<THEAD><TR><TD><B>Name</B></TD><TD><B>Expectation value</B></TD></TR></THEAD>
							<TBODY>
							  <xsl:for-each select="SCALAR_AVERAGE">
							  <TR>
							    <TD><B><xsl:value-of select="@name"/></B></TD>
							    <TD><xsl:value-of select="MEAN"/></TD>
							  </TR>
							  </xsl:for-each>
							  <xsl:for-each select="VECTOR_AVERAGE">
							    <xsl:for-each select="SCALAR_AVERAGE">
							    <TR>
							      <TD><B><xsl:value-of select="../@name"/></B>[<xsl:value-of select="@indexvalue"/>]</TD>
							      <TD><xsl:value-of select="MEAN"/></TD>
							    </TR>
                </xsl:for-each>
              </xsl:for-each>
            </TBODY>
          </TABLE>
				</xsl:for-each>
				</xsl:for-each>
				</xsl:for-each>
			</body>
		</html>
	</xsl:template>
</xsl:stylesheet>
