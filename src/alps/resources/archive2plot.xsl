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

<xsl:output method="xml"/>

<xsl:template match="/">
  <xsl:call-template name = "PrintXMLHeader"/>

  <plot name="{$plot_name}">
    <!-- general plot description -->
    <legend show="{$legend_show}"/>
    <xaxis label="{$xaxis_label}"/>
    <yaxis label="{$yaxis_label}"/>
    
    <!-- extract data -->      
    <xsl:call-template name="ExtractDataFromArchive"/>
  </plot>
</xsl:template>

<xsl:template name="ExtractValuesFromSimulation">

  <xsl:choose>
    <xsl:when test="$xaxis_type = 'INDEX'">
      <!-- xaxis value=index, extract yaxis value -->
      <xsl:call-template name="ExctractIndexAndAxisValue">
        <xsl:with-param name="axis_type" select="$yaxis_type"/>
        <xsl:with-param name="axis_name" select="$yaxis_name"/>
      </xsl:call-template>
    </xsl:when>

    <xsl:otherwise>
     <point>
      <!-- extract xaxis value -->
      <xsl:call-template name="ExtractAxisValue">
        <xsl:with-param name="axis_token" select="'x'"/>
        <xsl:with-param name="axis_type" select="$xaxis_type"/>
        <xsl:with-param name="axis_name" select="$xaxis_name"/>
        <xsl:with-param name="axis_index" select="$xaxis_index"/>
      </xsl:call-template>
  
      <!-- extract yaxis value -->
      <xsl:call-template name="ExtractAxisValue">
        <xsl:with-param name="axis_token" select="'y'"/>
        <xsl:with-param name="axis_type" select="$yaxis_type"/>
        <xsl:with-param name="axis_name" select="$yaxis_name"/>
        <xsl:with-param name="axis_index" select="$yaxis_index"/>
      </xsl:call-template>
      </point>
      <xsl:value-of select="$newline"/>
    </xsl:otherwise>
  </xsl:choose>  
    
</xsl:template>

<!-- template which extracts the value of a given parameter/scalar average/vector average -->
<xsl:template name="ExtractAxisValue">
  <xsl:param name="axis_token"/>
  <xsl:param name="axis_type"/>
  <xsl:param name="axis_name"/>
  <xsl:param name="axis_index" select="0"/>

  <xsl:choose>
  
  <!-- parameter -->
  <xsl:when test="$axis_type = 'PARAMETER'">
    <xsl:for-each select="PARAMETERS/PARAMETER">
      <xsl:if test = "@name = $axis_name">
      
        <!-- print parameter with axis tokens -->
        <xsl:call-template name="PrintToken">
           <xsl:with-param name="token" select="$axis_token"/>
        </xsl:call-template>
        <xsl:apply-templates/>
        <xsl:call-template name="PrintToken">
           <xsl:with-param name="token" select="$axis_token"/>
           <xsl:with-param name="close" select="'true'"/>
        </xsl:call-template>
        
      </xsl:if>
    </xsl:for-each>
  </xsl:when>

  <!-- scalar average -->
  <xsl:when test="$axis_type = 'SCALAR_AVERAGE'">
    <xsl:for-each select="AVERAGES/SCALAR_AVERAGE">
      <xsl:if test = "@name = $axis_name">
      
         <!-- print mean with axis tokens -->
         <xsl:call-template name="PrintToken">
           <xsl:with-param name="token" select="$axis_token"/>
         </xsl:call-template>
         <xsl:apply-templates select="MEAN"/>
         <xsl:call-template name="PrintToken">
           <xsl:with-param name="token" select="$axis_token"/>
           <xsl:with-param name="close" select="'true'"/>
         </xsl:call-template>
         
         <!-- print error with axis tokens -->
         <xsl:call-template name="PrintToken">
           <xsl:with-param name="token" select="$axis_token"/>
           <xsl:with-param name="error" select="'true'"/>
         </xsl:call-template>      
	 <xsl:choose>
  	   <xsl:when test="ERROR/@converged = 'no'">   
	     <xsl:value-of select="ERROR * 1000"/>
           </xsl:when>
  	   <xsl:when test="ERROR/@converged = 'maybe'">   
	     <xsl:value-of select="ERROR"/>
           </xsl:when>
	   <xsl:otherwise>
	     <xsl:value-of select="ERROR"/>
           </xsl:otherwise>
         </xsl:choose>
         <xsl:call-template name="PrintToken">
           <xsl:with-param name="token" select="$axis_token"/>
           <xsl:with-param name="close" select="'true'"/>
           <xsl:with-param name="error" select="'true'"/>
         </xsl:call-template>
         
      </xsl:if>
    </xsl:for-each>
  </xsl:when>

  <!-- vector average -->
  <xsl:when test="$axis_type = 'VECTOR_AVERAGE'">
    <xsl:for-each select="AVERAGES/VECTOR_AVERAGE">
      <xsl:if test = "@name = $axis_name">
        <xsl:for-each select="SCALAR_AVERAGE">
          <xsl:if test = "@indexvalue = $axis_index">

             <!-- print mean with axis tokens -->
             <xsl:call-template name="PrintToken">
               <xsl:with-param name="token" select="$axis_token"/>
             </xsl:call-template>
             <xsl:apply-templates select="MEAN"/>
             <xsl:call-template name="PrintToken">
               <xsl:with-param name="token" select="$axis_token"/>
               <xsl:with-param name="close" select="'true'"/>
             </xsl:call-template>
         
             <!-- print error with axis tokens -->
             <xsl:call-template name="PrintToken">
               <xsl:with-param name="token" select="$axis_token"/>
               <xsl:with-param name="error" select="'true'"/>
             </xsl:call-template>         
   	     <xsl:choose>
  	       <xsl:when test="ERROR/@converged = 'no'">   
	         <xsl:value-of select="ERROR * 1000"/>
               </xsl:when>
  	       <xsl:when test="ERROR/@converged = 'maybe'">   
	         <xsl:value-of select="ERROR"/>
               </xsl:when>
	       <xsl:otherwise>
	         <xsl:value-of select="ERROR"/>
               </xsl:otherwise>
             </xsl:choose>
             <xsl:call-template name="PrintToken">
               <xsl:with-param name="token" select="$axis_token"/>
               <xsl:with-param name="close" select="'true'"/>
               <xsl:with-param name="error" select="'true'"/>
             </xsl:call-template>

          </xsl:if>
        </xsl:for-each>
      </xsl:if>
    </xsl:for-each>
  </xsl:when>
            
  <!-- unknown type -->
  <xsl:otherwise>
    <xsl:text>ERROR: You are using an unknown </xsl:text>
    <xsl:value-of select="$axis_token"/>
    <xsl:text>-axis type!</xsl:text>
  </xsl:otherwise>
  
  </xsl:choose>
</xsl:template>

<!-- template which extracts index and value of a given vector average or histogram -->
<xsl:template name="ExctractIndexAndAxisValue">
  <xsl:param name="axis_type"/>
  <xsl:param name="axis_name"/>
  
  <xsl:choose>

  <!-- vector average -->
  <xsl:when test="$axis_type = 'VECTOR_AVERAGE'">
    <xsl:for-each select="AVERAGES/VECTOR_AVERAGE">
      <xsl:if test = "@name = $axis_name">
        <xsl:for-each select="SCALAR_AVERAGE">
          <xsl:call-template name="CheckIndexConstraints"/>
        </xsl:for-each>
      </xsl:if>
    </xsl:for-each>
  </xsl:when>
            
  <!-- histogram -->
  <xsl:when test="$axis_type = 'HISTOGRAM'">
    <xsl:for-each select="AVERAGES/HISTOGRAM">
      <xsl:if test = "@name = $axis_name">
        <xsl:for-each select="ENTRY">
          <xsl:call-template name="CheckHistogramIndexConstraints"/>
        </xsl:for-each>
      </xsl:if>
    </xsl:for-each>
  </xsl:when>

  <!-- unknown type -->
  <xsl:otherwise>
    <xsl:text>ERROR: xaxis INDEX type can only be combined with yaxis VECTOR_AVERAGE or HISTOGRAM type.</xsl:text>
  </xsl:otherwise>
  
  </xsl:choose>

</xsl:template>

<xsl:template name="ExtractSingleIndexAndValue">    
  <point>
    <x><xsl:value-of select="@indexvalue"/></x>
    <y><xsl:apply-templates select="MEAN"/></y>
    <dy>
      <xsl:choose>
        <xsl:when test="ERROR/@converged = 'no'">   
	  <xsl:value-of select="ERROR * 1000"/>
        </xsl:when>
  	<xsl:when test="ERROR/@converged = 'maybe'">   
	  <xsl:value-of select="ERROR"/>
        </xsl:when>
	<xsl:otherwise>
	  <xsl:value-of select="ERROR"/>
        </xsl:otherwise>
      </xsl:choose>
    </dy>
  </point>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template name="ExtractHistogramIndexAndValue">     
  <point>
    <x><xsl:value-of select="@indexvalue"/></x>
    <y><xsl:apply-templates select="VALUE"/></y>
  </point>
  <xsl:value-of select="$newline"/>
</xsl:template>

</xsl:stylesheet>
