<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:include href="helpers.xsl"/>
<xsl:output method="text"/>

<xsl:template match="plot">
  <xsl:text disable-output-escaping = "yes">set title "</xsl:text><xsl:value-of select='/plot/@name'/><xsl:text disable-output-escaping = "yes">"</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text disable-output-escaping = "yes">set xlabel "</xsl:text><xsl:value-of select='/plot/xaxis/@label'/><xsl:text disable-output-escaping = "yes">"</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text disable-output-escaping = "yes">set ylabel "</xsl:text><xsl:value-of select='/plot/yaxis/@label'/><xsl:text disable-output-escaping = "yes">"</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text disable-output-escaping = "yes">plot </xsl:text>
  <xsl:for-each select="set">
    <xsl:if test="position() > 1">
      <xsl:text disable-output-escaping = "yes">, </xsl:text>
    </xsl:if>
    <xsl:text disable-output-escaping = "yes">'-' </xsl:text>
    <xsl:choose>
      <xsl:when test="point/dx">
        <xsl:choose>
          <xsl:when test="point/dy">
            <xsl:text>using 1:2:3:4 notitle with xyerrorlines</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>using 1:2:3 notitle with xerrorlines</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <xsl:choose>
          <xsl:when test="point/dy">
            <xsl:text>using 1:2:3 notitle with yerrorlines</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>using 1:2 notitle with errorlines</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>    
  <xsl:value-of select="$newline"/>
  <xsl:for-each select="set">
    <xsl:apply-templates select="point">
      <xsl:sort select="x" data-type="number"/>
    </xsl:apply-templates>
    <xsl:text>e</xsl:text><xsl:value-of select="$newline"/>
  </xsl:for-each>    
</xsl:template>

<xsl:template match="plot/set/point">
  <xsl:if test="x and y">
    <xsl:value-of select="x"/>
    <xsl:text>	</xsl:text>
    <xsl:if test= "dx">
      <xsl:value-of select="dx"/>
      <xsl:text>	</xsl:text>
    </xsl:if>
    <xsl:value-of select="y"/>
    <xsl:text>	</xsl:text>
    <xsl:if test= "dy">
      <xsl:value-of select="dy"/>
    </xsl:if>
    <xsl:value-of select="$newline"/> 
  </xsl:if>
</xsl:template>

</xsl:stylesheet>
