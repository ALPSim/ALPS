<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
	<xsl:template match="/">

<!--
   Copyright (c) 2003-2010 Matthias Troyer (troyer@ethz.ch)
    
   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->
  

<xsl:for-each select="JOB/TASK">Task <xsl:value-of select="position()"/>: <xsl:value-of select="@status"/>, <xsl:for-each select="INPUT">input: <xsl:value-of select="@file"/></xsl:for-each><xsl:for-each select="OUTPUT">, output: <xsl:value-of select="@file"/></xsl:for-each><xsl:text>
</xsl:text></xsl:for-each>

    <xsl:for-each select="SIMULATION">
      <xsl:text>
Monte Carlo Simulation
======================

Parameters:

</xsl:text>
      <xsl:for-each select="PARAMETERS/PARAMETER">
        <xsl:value-of select="@name"/> = <xsl:apply-templates/>
        <xsl:text>
</xsl:text>
      </xsl:for-each>
      <xsl:text>

Execution information
---------------------
      
</xsl:text>
      <xsl:for-each select="MCRUN">
        <xsl:text>
Run </xsl:text><xsl:value-of select="position()"/>:
        <xsl:for-each select="EXECUTED">
Executed from <xsl:for-each select="FROM"><xsl:apply-templates/></xsl:for-each> to <xsl:for-each select="TO"><xsl:apply-templates/></xsl:for-each> on <xsl:for-each select="MACHINE/NAME"><xsl:apply-templates/></xsl:for-each>
        </xsl:for-each>
        <xsl:text>
</xsl:text>
      </xsl:for-each>
      <xsl:text>


Averages
========

Total:
------

</xsl:text>
        <xsl:for-each select="AVERAGES/SCALAR_AVERAGE">
          <xsl:value-of select="@name"/>: <xsl:value-of select="MEAN"/> +/- <xsl:value-of select="ERROR"/>, count=<xsl:value-of select="COUNT"/>, tau=<xsl:value-of select="AUTOCORR"/>
          <xsl:if test= "ERROR/@converged = 'maybe'"> <xsl:text>
</xsl:text>WARNING: check error convergence</xsl:if>
          <xsl:if test= "ERROR/@converged = 'no'"> <xsl:text>
</xsl:text>WARNING: ERRORS NOT CONVERGED!!!</xsl:if><xsl:text>
    
</xsl:text>
        </xsl:for-each>
        <xsl:for-each select="AVERAGES/VECTOR_AVERAGE/SCALAR_AVERAGE">
					 <xsl:value-of select="../@name"/>[<xsl:value-of select="@indexvalue"/>]: <xsl:value-of select="MEAN"/> +/- <xsl:value-of select="ERROR"/>, count=<xsl:value-of select="COUNT"/>, tau=<xsl:value-of select="AUTOCORR"/>
                     <xsl:if test= "ERROR/@converged = 'maybe'"> <xsl:text>
</xsl:text>WARNING: check error convergence</xsl:if>
                     <xsl:if test= "ERROR/@converged = 'no'"> <xsl:text>
</xsl:text>WARNING: ERRORS NOT CONVERGED!!!</xsl:if>
                     <xsl:text>
           
</xsl:text>
        </xsl:for-each>
        <xsl:text>
        
</xsl:text>
        <xsl:for-each select="MCRUN/AVERAGES">
Run :<xsl:text>
------

</xsl:text>
        <xsl:for-each select="SCALAR_AVERAGE">
          <xsl:value-of select="@name"/>: <xsl:value-of select="MEAN"/> +/- <xsl:value-of select="ERROR"/>, count=<xsl:value-of select="COUNT"/>, tau=<xsl:value-of select="AUTOCORR"/>        
          <xsl:if test= "ERROR/@converged = 'maybe'"> <xsl:text>
</xsl:text>WARNING: check error convergence</xsl:if>
          <xsl:if test= "ERROR/@converged = 'no'"> <xsl:text>
</xsl:text>WARNING: ERRORS NOT CONVERGED!!!</xsl:if><xsl:text>    
</xsl:text>
          <xsl:for-each select="BINNED">
             <xsl:text>    </xsl:text><xsl:value-of select="COUNT"/> bins: <xsl:value-of select="MEAN"/> +/- <xsl:value-of select="ERROR"/><xsl:text>    
</xsl:text>
          </xsl:for-each>
          <xsl:text>    
</xsl:text>
        </xsl:for-each>
        <xsl:for-each select="VECTOR_AVERAGE/SCALAR_AVERAGE">
					 <xsl:value-of select="../@name"/>[<xsl:value-of select="@indexvalue"/>]: <xsl:value-of select="MEAN"/> +/- <xsl:value-of select="ERROR"/>, count=<xsl:value-of select="COUNT"/>, tau=<xsl:value-of select="AUTOCORR"/>        
                     <xsl:if test= "ERROR/@converged = 'maybe'"> <xsl:text>
</xsl:text>WARNING: check error convergence</xsl:if>
          <xsl:if test= "ERROR/@converged = 'no'"> <xsl:text>
</xsl:text>WARNING: ERRORS NOT CONVERGED!!!</xsl:if><xsl:text>
</xsl:text>
          <xsl:for-each select="BINNED">
             <xsl:text>    </xsl:text><xsl:value-of select="COUNT"/> bins: <xsl:value-of select="MEAN"/> +/- <xsl:value-of select="ERROR"/>
             <xsl:text>    
</xsl:text>
          </xsl:for-each>
          <xsl:text>    
</xsl:text>
        </xsl:for-each>
        <xsl:text>
        
</xsl:text>
							</xsl:for-each>
</xsl:for-each>
	</xsl:template>
</xsl:stylesheet>
