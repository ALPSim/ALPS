// SPDX-License-Identifier: MIT
#include <alps/xml/matrix.hpp>
#include <algorithm>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <string>

template<class T>
std::string render(alps::numeric::matrix<T> const& matrix, bool member) {
    std::ostringstream stream;
    {
        alps::oxstream xml(stream);
        if (member) matrix.write_xml(xml);
        else xml << matrix;
    }
    auto text = stream.str();
    text.erase(std::remove_if(text.begin(), text.end(), [](char c) {
        return c == ' ' || c == '\n' || c == '\t' || c == '\r';
    }), text.end());
    return text;
}

int main() {
    alps::numeric::matrix<double> rectangle(2, 3);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 3; ++j) rectangle(i, j) = 3 * i + j;
    const std::string expected =
        "<MATRIXcols=\"3\"rows=\"2\"><ROW><ELEMENT>0</ELEMENT><ELEMENT>1</ELEMENT>"
        "<ELEMENT>2</ELEMENT></ROW><ROW><ELEMENT>3</ELEMENT><ELEMENT>4</ELEMENT>"
        "<ELEMENT>5</ELEMENT></ROW></MATRIX>";
    if (render(rectangle, true) != expected || render(rectangle, false) != expected)
        throw std::runtime_error("Historical rectangular matrix XML changed");
    alps::numeric::matrix<std::complex<double>> complex(1, 1);
    complex(0, 0) = {1., -2.};
    if (render(complex, false) !=
        "<MATRIXcols=\"1\"rows=\"1\"><ROW><ELEMENT>(1,-2)</ELEMENT></ROW></MATRIX>")
        throw std::runtime_error("Historical complex matrix XML changed");
    alps::numeric::matrix<std::string> escaped(1, 1);
    escaped(0, 0) = "<&>";
    if (render(escaped, true).find("<ELEMENT>&lt;&amp;&gt;</ELEMENT>") == std::string::npos)
        throw std::runtime_error("Matrix XML did not escape element text");
    alps::numeric::matrix<double> empty(0, 3);
    const auto empty_xml = render(empty, false);
    if (empty_xml != "<MATRIXcols=\"3\"rows=\"0\"/>"
        && empty_xml != "<MATRIXcols=\"3\"rows=\"0\"></MATRIX>")
        throw std::runtime_error("Empty matrix XML lost its shape");
}
