<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <xsl:template match="/SIMULATION/PARAMETERS">
    <xsl:text>{</xsl:text>
    <xsl:for-each select="PARAMETER">
      <xsl:text>'</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>' = '</xsl:text>
      <xsl:apply-templates/>
      <xsl:text>', '</xsl:text>
    </xsl:for-each>
    <xsl:text>}</xsl:text>
  </xsl:template>
</xsl:stylesheet>
