// Copyright (C) 2010 - 2012 by Andreas Hehn <hehn@phys.ethz.ch>.
// SPDX-License-Identifier: MIT
#ifndef ALPS_XML_MATRIX_HPP
#define ALPS_XML_MATRIX_HPP

#include <alps/numeric/matrix.hpp>
#include <alps/parser/xmlstream.h>
#include <sstream>

namespace alps::numeric {

// Preserve the historical MATRIX/ROW/ELEMENT representation and formatting.
template <typename T, typename MemoryBlock>
void matrix<T, MemoryBlock>::write_xml(oxstream& xml) const
{
    xml << start_tag("MATRIX") << attribute("cols", num_cols())
        << attribute("rows", num_rows());
    for (size_type i = 0; i < num_rows(); ++i) {
        xml << start_tag("ROW");
        for (size_type j = 0; j < num_cols(); ++j) {
            std::stringstream value;
            value << (*this)(i, j);
            xml << start_tag("ELEMENT") << alps::convert(value.str()) << end_tag("ELEMENT");
        }
        xml << end_tag("ROW");
    }
    xml << end_tag("MATRIX");
}

template <typename T, typename MemoryBlock>
oxstream& operator<<(oxstream& xml, matrix<T, MemoryBlock> const& value)
{
    value.write_xml(xml);
    return xml;
}

} // namespace alps::numeric
#endif
