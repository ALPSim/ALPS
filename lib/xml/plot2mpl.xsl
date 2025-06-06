<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!--
   Copyright (c) 2003-2010 Matthias Troyer <troyer@ethz.ch>
  
   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the “Software”),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:
  
   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.
  
   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
  -->
  
<xsl:include href="helpers.xsl"/>
<xsl:output method="text"/>

<xsl:template match="plot">
  <xsl:text>from pylab import *</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text>from numpy import *</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text>subplot(111)</xsl:text><xsl:value-of select="$newline"/>
  
  <xsl:for-each select="set">
    <xsl:text>data = zip(*[</xsl:text><xsl:value-of select="$newline"/>
    <xsl:if test= "point">
      <xsl:apply-templates select="point">
        <xsl:sort select="x" data-type="number"/>
      </xsl:apply-templates>
    </xsl:if>
    <xsl:text>])</xsl:text><xsl:value-of select="$newline"/>
    <xsl:text>if len(data):</xsl:text><xsl:value-of select="$newline"/><xsl:text>     errorbar(data[0],data[1]</xsl:text>
    <xsl:if test= "point/dy">
      <xsl:text>,yerr=data[3]</xsl:text>
    </xsl:if>
    <xsl:if test= "point/dx">
      <xsl:text>,xerr=data[2]</xsl:text>
    </xsl:if>
    <xsl:text>)</xsl:text><xsl:value-of select="$newline"/>
    
  </xsl:for-each> 
  
 
  <xsl:text>title("</xsl:text><xsl:value-of select='/plot/@name'/><xsl:text>")</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text>xlabel("</xsl:text><xsl:value-of select='/plot/xaxis/@label'/><xsl:text>")</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text>ylabel("</xsl:text><xsl:value-of select='/plot/yaxis/@label'/><xsl:text>")</xsl:text><xsl:value-of select="$newline"/>
  <xsl:choose>
    <xsl:when test="plot/legend/@show = 'true'">
      <xsl:text>legend()</xsl:text><xsl:value-of select="$newline"/>
    </xsl:when>
  </xsl:choose>

</xsl:template>

<xsl:template match="plot/set/point">
  <xsl:text>[</xsl:text>
  <xsl:value-of select="x"/>
  <xsl:text>, </xsl:text>
  <xsl:value-of select="y"/>
  <xsl:text>,	</xsl:text>
  <xsl:choose>
    <xsl:when test="dx">
      <xsl:value-of select="dx"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>0. </xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>,	</xsl:text>
  <xsl:choose>
    <xsl:when test="dy">
      <xsl:value-of select="dy"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>0. </xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>], </xsl:text><xsl:value-of select="$newline"/>
</xsl:template>

</xsl:stylesheet>
