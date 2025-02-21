<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                              xmlns:xslout="http://xml.comp-phys.org/schema/plot">

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
  
<xsl:namespace-alias stylesheet-prefix="xslout" result-prefix="xsl"/>
<xsl:output method="xml" indent="yes"/>
<xsl:include href="helpers.xsl"/>

<xsl:template match="/">
  <xslout:stylesheet version="1.0">
    <xsl:value-of select="$newline"/>
    <xsl:value-of select="$newline"/>
  
    <xslout:include href="archive2plot.xsl"/>
    <xsl:value-of select="$newline"/>

    <!-- general plot description -->
    <xslout:variable name="plot_name">
      <xsl:value-of select='plot/@name'/>
    </xslout:variable>
    <xsl:value-of select="$newline"/>
    
    <xslout:variable name="legend_show">
      <xsl:value-of select='plot/legend/@show'/>
    </xslout:variable>
    <xsl:value-of select="$newline"/>
    
    <!-- xaxis description -->
    <xslout:variable name="xaxis_label">
      <xsl:value-of select='plot/xaxis/@label'/>
    </xslout:variable>
    <xsl:value-of select="$newline"/>
    
    <xslout:variable name="xaxis_type">
      <xsl:value-of select='plot/xaxis/@type'/>
    </xslout:variable>
    <xsl:value-of select="$newline"/>

    <xslout:variable name="xaxis_name">
      <xsl:choose>
        <xsl:when test="plot/xaxis/@name">
          <xsl:value-of select='plot/xaxis/@name'/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select='plot/xaxis/@label'/>
        </xsl:otherwise>  
      </xsl:choose>
    </xslout:variable>
    <xsl:value-of select="$newline"/>
    
    <xsl:choose>
      <xsl:when test="plot/xaxis/@index">
        <xslout:variable name="xaxis_index">
        <xsl:value-of select='plot/xaxis/@index'/>
        </xslout:variable>
      </xsl:when>
      <xsl:otherwise>
        <xslout:variable name="xaxis_index"/>  
      </xsl:otherwise>  
    </xsl:choose>
    <xsl:value-of select="$newline"/>

    <!-- yaxis description -->
    <xslout:variable name="yaxis_label">
      <xsl:value-of select='plot/yaxis/@label'/>
    </xslout:variable>
    <xsl:value-of select="$newline"/>
    
    <xslout:variable name="yaxis_type">
      <xsl:value-of select='plot/yaxis/@type'/>
    </xslout:variable>
    <xsl:value-of select="$newline"/>

    <xslout:variable name="yaxis_name">
      <xsl:choose>
        <xsl:when test="plot/yaxis/@name">
          <xsl:value-of select='plot/yaxis/@name'/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select='plot/yaxis/@label'/>
        </xsl:otherwise>  
      </xsl:choose>
    </xslout:variable>
    <xsl:value-of select="$newline"/>
    
    <xsl:choose>
      <xsl:when test="plot/yaxis/@index">
        <xslout:variable name="yaxis_index">
        <xsl:value-of select='plot/yaxis/@index'/>
        </xslout:variable>
      </xsl:when>
      <xsl:otherwise>
        <xslout:variable name="yaxis_index"/>  
      </xsl:otherwise>  
    </xsl:choose>
    <xsl:value-of select="$newline"/>
    <xsl:value-of select="$newline"/>
    <xsl:value-of select="$newline"/>
  
  
    <!-- constraints -->
    <xsl:for-each select="plot">

    <!-- global index constraints -->
    <xslout:template name="CheckIndexConstraints"><xsl:value-of select="$newline"/>
      <xsl:text disable-output-escaping = "yes">&lt;xslout:if test='true()</xsl:text>
      <xsl:for-each select="constraint">
        <xsl:if test="@type = 'INDEX'">
          <xsl:text> and @indexvalue</xsl:text>
          <xsl:value-of select="@condition"/>
        </xsl:if>    
      </xsl:for-each>
      <xsl:text disable-output-escaping = "yes">'&gt;</xsl:text><xsl:value-of select="$newline"/>
      <xslout:call-template name="ExtractSingleIndexAndValue"/><xsl:value-of select="$newline"/>
      <xsl:text disable-output-escaping = "yes">&lt;/xslout:if&gt;</xsl:text><xsl:value-of select="$newline"/>
    </xslout:template><xsl:value-of select="$newline"/>
    
    <xsl:value-of select="$newline"/>
    <xsl:value-of select="$newline"/> 

    <xslout:template name="CheckHistogramIndexConstraints"><xsl:value-of select="$newline"/>
    <xsl:text disable-output-escaping = "yes">&lt;xslout:if test='true()</xsl:text>
      <xsl:for-each select="constraint">
        <xsl:if test="@type = 'INDEX'">
          <xsl:text> and @indexvalue</xsl:text>
          <xsl:value-of select="@condition"/>
        </xsl:if>
      </xsl:for-each>
      <xsl:text disable-output-escaping = "yes">'&gt;</xsl:text><xsl:value-of select="$newline"/>
      <xslout:call-template name="ExtractHistogramIndexAndValue"/><xsl:value-of select="$newline"/>
      <xsl:text disable-output-escaping = "yes">&lt;/xslout:if&gt;</xsl:text><xsl:value-of select="$newline"/>
    </xslout:template><xsl:value-of select="$newline"/>

    <xsl:value-of select="$newline"/>
    <xsl:value-of select="$newline"/>
 
       
    <!-- global constraints -->
    <xslout:template name="ExtractDataFromArchive"><xsl:value-of select="$newline"/>
      <xslout:call-template name="DoLoop"/>
    </xslout:template> <!-- CheckConstraints -->
      
    <xsl:value-of select="$newline"/>
    <xsl:value-of select="$newline"/>

    <!-- global loop -->
    <xslout:template name="DoLoop"><xsl:value-of select="$newline"/>
      <xsl:text disable-output-escaping = "yes">&lt;!-- global loop --&gt;</xsl:text>
      <xsl:value-of select="$newline"/>

      <xslout:call-template name="OpenSet">
        <xslout:with-param name="empty_set" select="true()"/>
      </xslout:call-template>
        
      <!-- iterate all simulation data -->      
      <xslout:for-each select="//SIMULATION">
           
        <xsl:for-each select="for-each">
          <xsl:choose>
            <xsl:when test="@type = 'PARAMETER' or not(@type)">
              <xsl:text disable-output-escaping = "yes">&lt;xslout:sort select="PARAMETERS/PARAMETER[@name='</xsl:text>
              <xsl:call-template name="NameOrLabel"/>
              <xsl:text disable-output-escaping = "yes">']" data-type="number"/&gt;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <!-- no other types implemented yet -->
            </xsl:otherwise>
          </xsl:choose>
        </xsl:for-each>  
           
        <!-- new loop increment? -->
        <xsl:for-each select="for-each">
          <xsl:choose>
            <xsl:when test="@type = 'PARAMETER' or not(@type)">
              <xsl:text disable-output-escaping = "yes">&lt;xslout:variable name="last" select="PARAMETERS/PARAMETER[@name='</xsl:text>
              <xsl:call-template name="NameOrLabel"/>
              <xsl:text disable-output-escaping = "yes">']"/&gt;</xsl:text>                   
              <xsl:text disable-output-escaping = "yes">&lt;xslout:if test="not(preceding-sibling::SIMULATION[PARAMETERS/PARAMETER[@name='</xsl:text>
              <xsl:call-template name="NameOrLabel"/>
              <xsl:text disable-output-escaping = "yes">'] = $last])"&gt;</xsl:text>
              <xsl:text disable-output-escaping = "yes">&lt;xslout:call-template name="CloseSet"/&gt;</xsl:text>
              <xsl:text disable-output-escaping = "yes">&lt;xslout:call-template name="OpenSet"/&gt;</xsl:text>
              <xsl:text disable-output-escaping = "yes">&lt;/xslout:if&gt;</xsl:text>
            </xsl:when>
           </xsl:choose>
        </xsl:for-each>  <!-- loop -->
        
        <xsl:call-template name="CheckSetConstraints"/>
      </xslout:for-each>  <!-- SIMULATION -->
        
      <xslout:call-template name="CloseSet"/>
    </xslout:template> <!-- DoLoop -->
                      
      <xsl:value-of select="$newline"/>
      <xsl:value-of select="$newline"/>
      
      <xslout:template name="OpenSet"><xsl:value-of select="$newline"/>
        <xslout:param name="empty_set" select="false()"/>
      
        <xslout:choose>
        <xslout:when test="not($empty_set)">
        <!-- generate label from set constraint -->
        <xsl:value-of select="$newline"/>
        <xsl:text disable-output-escaping = "yes">&lt;!-- set constraints --&gt;</xsl:text>
        <xsl:value-of select="$newline"/>
        <xslout:text disable-output-escaping = "yes">&lt;set label="</xslout:text>
          <xsl:value-of select="@label"/>
          
            <!-- check loops -->
             <xsl:for-each select="//for-each">
              <xsl:call-template name="LabelOrName"/>
              <xsl:text>=</xsl:text>
              <xsl:choose>
                <xsl:when test="@type = 'PARAMETER' or not(@type)">
                  <xsl:text disable-output-escaping = "yes">&lt;xslout:value-of select="PARAMETERS/PARAMETER[@name='</xsl:text>
                  <xsl:call-template name="NameOrLabel"/>
                  <xsl:text disable-output-escaping = "yes">']"/&gt;</xsl:text>
                </xsl:when>
              </xsl:choose>

              <xsl:if test="//constraint"><xsl:text> and </xsl:text></xsl:if> 
            </xsl:for-each>
         
            <!-- check constraints -->
            <xsl:for-each select="//constraint">
              <xsl:call-template name="LabelOrName"/>
              <xsl:value-of select="@condition"/>
              <xsl:if test="not(position() = last())"><xsl:text> and </xsl:text></xsl:if> 
            </xsl:for-each>
            
        <xslout:text disable-output-escaping = "yes">"&gt;</xslout:text>
        <xsl:value-of select="$newline"/>
        </xslout:when>
        <xslout:otherwise>
          <xslout:text disable-output-escaping = "yes">&lt;set&gt;</xslout:text>
        </xslout:otherwise>
        </xslout:choose>
        
        <xslout:value-of select="$newline"/>
        <xsl:value-of select="$newline"/>
      </xslout:template><!-- OpenSet -->

      <xsl:value-of select="$newline"/>
      <xsl:value-of select="$newline"/>
 
      <xslout:template name="CloseSet"><xsl:value-of select="$newline"/>
          <xslout:text disable-output-escaping = "yes">&lt;/set&gt;</xslout:text><xsl:value-of select="$newline"/>
          <xslout:value-of select="$newline"/><xsl:value-of select="$newline"/>
      </xslout:template><!-- CloseSet -->

      <xsl:value-of select="$newline"/>
      <xsl:value-of select="$newline"/>
      
    </xsl:for-each> <!-- plot constraints -->
    
   <xsl:value-of select="$newline"/>
  </xslout:stylesheet>
  
</xsl:template>


<xsl:template name="CheckSetConstraints">
  <!-- implement set constraints -->
  <xsl:if test="//constraint">      
    <xsl:text disable-output-escaping = "yes">&lt;xslout:if test='true()</xsl:text>
  </xsl:if>  
  <xsl:for-each select="//constraint">
    <xsl:if test="@condition">
      <xsl:if test="@type = 'PARAMETER' or @type = 'SCALAR_AVERAGE' or not(@type)"><xsl:text> and </xsl:text></xsl:if>  
      <xsl:choose>
        <xsl:when test="@type = 'PARAMETER' or not(@type)">
          <xsl:text disable-output-escaping = "yes">PARAMETERS/PARAMETER[@name="</xsl:text>
          <xsl:call-template name="NameOrLabel"/>
          <xsl:text disable-output-escaping = "yes">"]</xsl:text>
          <xsl:value-of select="@condition"/>
        </xsl:when>
        <xsl:when test="@type = 'SCALAR_AVERAGE'">
          <xsl:text disable-output-escaping = "yes">AVERAGES/SCALAR_AVERAGE[@name="</xsl:text>
          <xsl:call-template name="NameOrLabel"/>
          <xsl:text disable-output-escaping = "yes">"]/MEAN</xsl:text>
          <xsl:value-of select="@condition"/>
        </xsl:when>
      </xsl:choose>
    </xsl:if> 
  </xsl:for-each>  <!-- constraint -->
  <xsl:if test="//constraint">      
    <xsl:text disable-output-escaping = "yes">'&gt;</xsl:text><xsl:value-of select="$newline"/>
  </xsl:if>  

  <!-- extract data -->
  <xslout:call-template name="ExtractValuesFromSimulation"/><xsl:value-of select="$newline"/>

  <xsl:if test="//constraint">      
    <xsl:text disable-output-escaping = "yes">&lt;/xslout:if&gt;</xsl:text><xsl:value-of select="$newline"/>
  </xsl:if>  
</xsl:template><!-- CheckSetConstraints -->
 

</xsl:stylesheet>
