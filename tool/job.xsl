<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <xsl:template match="/">
    <html>
      <head/>
      <body>
        <xsl:for-each select="JOB">
        <H1>ALPS Job Description File</H1>
        <xsl:for-each select="OUTPUT">
	<p>Output to new job file 
          <A href="{@file}">
            <xsl:value-of select="@file"/>
          </A></p>
	</xsl:for-each>
	<table border="1">
	  <thead><tr><td><b>Status</b></td><td><b>Input file</b></td>
                     <td><b>Output file</b></td></tr></thead>
	  <tbody>
	    <xsl:for-each select="TASK">
	      <tr>
	        <td><xsl:value-of select="@status"/></td>
	        <td>
		  <xsl:for-each select="INPUT">
                    <A href="{@file}">
		      <xsl:value-of select="@file"/>
                    </A>
		  </xsl:for-each>
	        </td>
	        <td>
		  <xsl:for-each select="OUTPUT">
                    <A href="{@file}">
		        <xsl:value-of select="@file"/>
                    </A>
		  </xsl:for-each>
	        </td>
	      </tr>
	    </xsl:for-each>
	  </tbody>
	</table>
	    </xsl:for-each>
      </body>
    </html>
  </xsl:template>
</xsl:stylesheet>
