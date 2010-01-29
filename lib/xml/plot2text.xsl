<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!--
   Copyright (c) 2003-2010 Matthias Troyer <troyer@ethz.ch>
  
   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->
  
<xsl:include href="helpers.xsl"/>
<xsl:output method="text"/>

<xsl:template match="plot">
  <xsl:for-each select="set">
    <xsl:if test= "../legend/@show = 'true' and @label">
      <xsl:value-of select="@label"/>
      <xsl:value-of select="$newline"/>
    </xsl:if>

    <xsl:apply-templates select="point">
      <xsl:sort select="x" data-type="number"/>
    </xsl:apply-templates>
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
