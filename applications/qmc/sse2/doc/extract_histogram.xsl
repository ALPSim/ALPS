<!--
   Copyright (c) 2003-2010 Matthias Troyer (troyer@ethz.ch)
    
   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->

<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE xsl:stylesheet [
 <!ENTITY nbsp "&#160;">
]>

    <xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"
    >
    <xsl:output method="text" encoding="utf-8" indent="no"
       />
      

<xsl:template match="SIMULATION" >
<xsl:for-each select="AVERAGES/VECTOR_AVERAGE[@name='Histogram']/SCALAR_AVERAGE">
<xsl:value-of select="@indexvalue"/><xsl:text> </xsl:text><xsl:value-of select="MEAN"/><xsl:text>
</xsl:text>
</xsl:for-each>
</xsl:template>



</xsl:stylesheet>
