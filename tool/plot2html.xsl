<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:include href="helpers.xsl"/>
<xsl:output method="html"/>

<xsl:template match="plot">
  <html>
   <head/>
   <body>
    <h1><xsl:value-of select="@name"/></h1>
    <xsl:value-of select="yaxis/@label"/><xsl:text> versus </xsl:text>
    <xsl:value-of select="xaxis/@label"/>

    <xsl:for-each select="set">
      <xsl:if test= "../legend/@show = 'true' and @label">
        <h2><xsl:value-of select="@label"/></h2>
        <xsl:value-of select="$newline"/>
      </xsl:if>

      <table border="1">
      <xsl:apply-templates select="point">
        <xsl:sort select="x" data-type="number"/>
      </xsl:apply-templates>
      </table>
      
    </xsl:for-each>

   </body>
  </html>
</xsl:template>

<xsl:template match="plot/set/point">
 <xsl:if test="position()=1">
 <tr>
   <td><b><xsl:value-of select="../../xaxis/@label"/></b></td>
   <xsl:if test= "dx">
      <td><b><xsl:value-of select="../../xaxis/@label"/>
          <xsl:text> Error</xsl:text></b></td>
   </xsl:if>
   <td><b><xsl:value-of select="../../yaxis/@label"/></b></td>
   <xsl:if test= "dy">
      <td><b><xsl:value-of select="../../yaxis/@label"/>
          <xsl:text> Error</xsl:text></b></td>
   </xsl:if>
 </tr>
 </xsl:if>
 <tr>
   <td><xsl:value-of select="x"/></td>
   <xsl:if test= "dx">
      <td><xsl:value-of select="dx"/></td>
   </xsl:if>
   <td><xsl:value-of select="y"/></td>
   <xsl:if test= "dy">
      <td><xsl:value-of select="dy"/></td>
    </xsl:if>
 </tr> 
</xsl:template>

</xsl:stylesheet>
