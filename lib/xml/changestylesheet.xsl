<?xml version="1.0" encoding="UTF-8"?> 
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="xml"/>

  <xsl:template match="/">
   <xsl:text disable-output-escaping = "yes">
&lt;?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>
</xsl:text>   
    <xsl:copy-of select="."/>
  </xsl:template>

</xsl:stylesheet>
