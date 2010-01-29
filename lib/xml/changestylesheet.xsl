<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!--
   Copyright (c) 2003-2010 Matthias Troyer (troyer@ethz.ch)
    
   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->
  
  <xsl:output method="xml"/>

  <xsl:template match="/">
   <xsl:text disable-output-escaping = "yes">
&lt;?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>
</xsl:text>   
    <xsl:copy-of select="."/>
  </xsl:template>

</xsl:stylesheet>
