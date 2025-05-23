<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!--
   Copyright (c) 2003-2010 Matthias Troyer (troyer@ethz.ch)
    
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
  <xsl:text disable-output-escaping = "yes">set title "</xsl:text><xsl:value-of select='/plot/@name'/><xsl:text disable-output-escaping = "yes">"</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text disable-output-escaping = "yes">set xlabel "</xsl:text><xsl:value-of select='/plot/xaxis/@label'/><xsl:text disable-output-escaping = "yes">"</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text disable-output-escaping = "yes">set ylabel "</xsl:text><xsl:value-of select='/plot/yaxis/@label'/><xsl:text disable-output-escaping = "yes">"</xsl:text><xsl:value-of select="$newline"/>
  <xsl:text disable-output-escaping = "yes">plot </xsl:text>
  <xsl:for-each select="set">
    <xsl:if test= "point">
      <xsl:if test="preceding-sibling::set/point">
        <xsl:text disable-output-escaping = "yes">, </xsl:text>
      </xsl:if>
      <xsl:text disable-output-escaping = "yes">'-' </xsl:text>
      <xsl:choose>
        <xsl:when test="point/dx">
          <xsl:choose>
            <xsl:when test="point/dy">
              <xsl:text>using 1:2:3:4 </xsl:text>
              <xsl:choose>
                <xsl:when test="../legend/@show = 'true' and @label">
                  <xsl:text disable-output-escaping = "yes">title "</xsl:text>
                  <xsl:value-of select="@label"/>
                  <xsl:text disable-output-escaping = "yes">" </xsl:text>
                </xsl:when>
                <xsl:otherwise>
                  <xsl:text>notitle </xsl:text>
                </xsl:otherwise>
              </xsl:choose>
              <xsl:text>with xyerrorlines</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>using 1:2:3 </xsl:text>
              <xsl:choose>
                <xsl:when test="../legend/@show = 'true' and @label">
                  <xsl:text disable-output-escaping = "yes">title "</xsl:text>
                  <xsl:value-of select="@label"/>
                  <xsl:text disable-output-escaping = "yes">" </xsl:text>
                </xsl:when>
                <xsl:otherwise>
                  <xsl:text>notitle </xsl:text>
                </xsl:otherwise>
              </xsl:choose>
              <xsl:text>with xerrorlines</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
        </xsl:when>
        <xsl:otherwise>
          <xsl:choose>
            <xsl:when test="point/dy">
              <xsl:text>using 1:2:3 </xsl:text>
              <xsl:choose>
                <xsl:when test="../legend/@show = 'true' and @label">
                  <xsl:text disable-output-escaping = "yes">title "</xsl:text>
                  <xsl:value-of select="@label"/>
                  <xsl:text disable-output-escaping = "yes">" </xsl:text>
                </xsl:when>
                <xsl:otherwise>
                  <xsl:text>notitle </xsl:text>
                </xsl:otherwise>
              </xsl:choose>
              <xsl:text> with yerrorlines</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>using 1:2 </xsl:text>
              <xsl:choose>
                <xsl:when test="../legend/@show = 'true' and @label">
                  <xsl:text disable-output-escaping = "yes">title "</xsl:text>
                  <xsl:value-of select="@label"/>
                  <xsl:text disable-output-escaping = "yes">" </xsl:text>
                </xsl:when>
                <xsl:otherwise>
                  <xsl:text>notitle </xsl:text>
                </xsl:otherwise>
              </xsl:choose>
              <xsl:text>with linespoints</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:if>
  </xsl:for-each>    
  <xsl:value-of select="$newline"/>
  <xsl:for-each select="set">
    <xsl:if test= "point">
      <xsl:apply-templates select="point">
        <xsl:sort select="x" data-type="number"/>
      </xsl:apply-templates>
      <xsl:text>e</xsl:text><xsl:value-of select="$newline"/>
    </xsl:if>
  </xsl:for-each>    
</xsl:template>

<xsl:template match="plot/set/point">
  <xsl:if test="x and y">
    <xsl:value-of select="x"/>
    <xsl:text>	</xsl:text>
    <xsl:if test="dx">
      <xsl:value-of select="dx"/>
      <xsl:text>	</xsl:text>
    </xsl:if>
    <xsl:value-of select="y"/>
    <xsl:text>	</xsl:text>
    <xsl:if test="dy">
      <xsl:value-of select="dy"/>
    </xsl:if>
    <xsl:value-of select="$newline"/> 
  </xsl:if>
</xsl:template>

</xsl:stylesheet>
