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
