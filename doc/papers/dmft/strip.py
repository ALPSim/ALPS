#  Copyright Emanuel Gull 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


stripped=open("refs_shortened.bib",'w')
for line in open("refs.bib"):
  if "url" in line or "URL" in line or "Url" in line or "eprint" in line:
      continue
  else:
    line=line.replace("Physical Review Letters","Phys. Rev. Lett.")
    line=line.replace("Zeitschrift fur Physik B Condensed Matter","Z. Phys. B")
    line=line.replace("Zeitschrift f{\"u}r Physik B Condensed Matter","Z. Phys. B")
    line=line.replace("Journal of Applied Physics","J. Appl. Phys.")
    line=line.replace("Journal of Physics: Condensed Matter","J. Phys. Condens. Matter")
    line=line.replace("Journal of Physics A: Mathematical and General","J. Phys. A")
    line=line.replace("Journal of the Physical Society of Japan","J. Phys. Soc. Jpn.")
    line=line.replace("New Journal of Physics","New J. Phys.")
    line=line.replace("Journal of Magnetism and Magnetic Materials","J. Magn. Magn. Mater.")
    line=line.replace("Physical Review B (Condensed Matter and Materials Physics)", "Phys. Rev. B")
    line=line.replace("Reviews of Modern Physics", "Rev. Mod. Phys.")
    line=line.replace("EPL (Europhysics Letters)", "Europhys. Lett.")
    line=line.replace("EPL", "Europhys. Lett.")
    line=line.replace("Progress of Theoretical Physics", "Prog. Theor. Phys.")
    stripped.write(line)
