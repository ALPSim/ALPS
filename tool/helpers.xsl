<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:variable name="newline">
<xsl:text>
</xsl:text>
</xsl:variable>

<xsl:template name="PrintXMLHeader">
   <xsl:text disable-output-escaping = "yes">
&lt;?xml-stylesheet type="text/xsl" href="http://xml.comp-phys.org/2003/4/plot2html.xsl"?>
   </xsl:text>   
</xsl:template>

<xsl:template name="PrintToken">
  <xsl:param name="token"/>
  <xsl:param name="close" select="'false'"/>
  <xsl:param name="error" select="'false'"/>
  
  <xsl:text disable-output-escaping = "yes">&lt;</xsl:text>
  <xsl:if test="$close = 'true'"><xsl:text>/</xsl:text></xsl:if>
  <xsl:if test="$error = 'true'"><xsl:text>d</xsl:text></xsl:if>
  <xsl:value-of select="$token"/>
  <xsl:text disable-output-escaping = "yes">&gt;</xsl:text>
</xsl:template>

<xsl:template name="PrintVariable">
  <xsl:param name="name"/>
  <xsl:param name="value"/>
  
  <xsl:text disable-output-escaping = "yes">&lt;</xsl:text>
</xsl:template>

</xsl:stylesheet>