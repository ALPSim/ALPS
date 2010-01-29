<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!--
   Copyright (c) 2003-2010 Matthias Troyer (troyer@ethz.ch)
    
   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->
  
<xsl:variable name="newline">
<xsl:text>
</xsl:text>
</xsl:variable>

<xsl:template name="PrintXMLHeader">
   <xsl:text disable-output-escaping = "yes">
&lt;?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>
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

<xsl:template name="NameOrLabel">
  <xsl:choose>
    <xsl:when test="@name"><xsl:value-of select="@name"/></xsl:when>
    <xsl:otherwise><xsl:value-of select="@label"/></xsl:otherwise>
  </xsl:choose>   
</xsl:template>

<xsl:template name="LabelOrName">
  <xsl:choose>
    <xsl:when test="@label"><xsl:value-of select="@label"/></xsl:when>
    <xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
  </xsl:choose>  
</xsl:template>

<xsl:template name="Print_XMGRACE_PlotHeader">
<xsl:text disable-output-escaping = "yes">
# Grace project file
#
@g0 on
@with g0
@    frame linewidth 2.0
@    title "</xsl:text>   
<xsl:value-of select='/plot/@name'/>
<xsl:text disable-output-escaping = "yes">"
@    xaxis  label "</xsl:text>
<xsl:value-of select='/plot/xaxis/@label'/>
<xsl:text disable-output-escaping = "yes">"
@    xaxis  label char size 1.500000
@    xaxis  ticklabel char size 1.250000
@    xaxis  tick minor ticks 4
@    yaxis  label "</xsl:text>
<xsl:value-of select='/plot/yaxis/@label'/>
<xsl:text disable-output-escaping = "yes">"
@    yaxis  label char size 1.500000
@    yaxis  ticklabel char size 1.250000
@    yaxis  tick minor ticks 4</xsl:text>
</xsl:template>

<xsl:template name="Print_XMGRACE_SetHeader">
<xsl:text disable-output-escaping = "yes">
@    s</xsl:text><xsl:number/><xsl:text> symbol 1
@    s</xsl:text><xsl:number/><xsl:text> symbol size 0.500000
@    s</xsl:text><xsl:number/><xsl:text> line type 1
@target G0.S</xsl:text><xsl:number/><xsl:text>
@type </xsl:text>   
<xsl:if test="./point/x">
  <xsl:text>x</xsl:text>
</xsl:if>
<xsl:if test="./point/dx">
  <xsl:text>dx</xsl:text>
</xsl:if>
<xsl:if test="./point/y">
  <xsl:text>y</xsl:text>
</xsl:if>
<xsl:if test="./point/dy">
  <xsl:text>dy</xsl:text>
</xsl:if>
<xsl:text>
</xsl:text>
</xsl:template>

</xsl:stylesheet>
